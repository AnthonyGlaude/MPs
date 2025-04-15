import os
import csv
import re
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster

# -----------------------------------------------------------
# Étape 1 : Extraction des meilleures lignes par variant pour chaque groupe
# -----------------------------------------------------------
def extract_best_entries(tsv_file, gene_dict):
    """
    Lit le fichier TSV d'OpenProt et, pour chaque gène défini dans gene_dict,
    ne conserve que les lignes dont le gene symbol commence par le deuxième terme (celui avec le "P").
    Pour chaque variant (ex. GAPDHP1, GAPDHP2, …), on sélectionne la ligne avec le meilleur MS score,
    et en cas d'égalité, celle avec la plus grande longueur protéique.
    
    Les résultats sont sauvegardés dans des dossiers distincts pour chaque groupe,
    avec un fichier nommé "<groupe>_transcripts.txt" contenant le header :
       Gene    Transcript_ID    PA
    (Ici, "Gene" correspond au nom du groupe, par exemple "GAPDH")
    
    :param tsv_file: chemin vers le fichier 'openprot_database.tsv'
    :param gene_dict: dictionnaire des gènes, par exemple :
                      {
                          "GAPDH": ["GAPDH", "GAPDHP"],
                          "PPIA": ["PPIA", "PPIAP"],
                          "PHB":  ["PHB",  "PHBP"],
                          "EEF1A1": ["EEF1A1", "EEF1A1P"],
                          "HMGB1":  ["HMGB1", "HMGB1P"]
                      }
    """
    best_entries = {}

    with open(tsv_file, 'r', newline='') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            gene_symbol = row["gene symbol"]

            for group, terms in gene_dict.items():
                if len(terms) < 2:
                    continue
                search_term = terms[1]  # par exemple "GAPDHP"
                pattern = f"^{re.escape(search_term)}(\\d+)$"
                match = re.match(pattern, gene_symbol)
                if match:
                    variant_number = match.group(1)
                    variant = f"P{variant_number}"
                    try:
                        ms_score = float(row["MS score"])
                    except ValueError:
                        ms_score = 0.0
                    try:
                        protein_length = int(row["protein length (a.a.)"])
                    except ValueError:
                        protein_length = 0

                    if group not in best_entries:
                        best_entries[group] = {}
                    if variant not in best_entries[group]:
                        best_entries[group][variant] = (row, ms_score, protein_length)
                    else:
                        _, best_ms, best_length = best_entries[group][variant]
                        if (ms_score > best_ms) or (ms_score == best_ms and protein_length > best_length):
                            best_entries[group][variant] = (row, ms_score, protein_length)
                    break

    for group, variant_dict in best_entries.items():
        os.makedirs(group, exist_ok=True)
        output_file = os.path.join(group, f"{group}_transcripts.txt")
        with open(output_file, 'w') as out_file:
            out_file.write("Gene\tTranscript_ID\tPA\n")
            for variant in sorted(variant_dict.keys(), key=lambda v: int(v[1:])):
                row, _, _ = variant_dict[variant]
                # Vous pouvez choisir d'écrire le nom du groupe ou la valeur extraite ;
                # ici nous utilisons la valeur extraite (par exemple "GAPDHP1")
                out_file.write(f"{row['gene symbol']}\t{row['transcript accession']}\t{row['protein accession']}\n")

# -----------------------------------------------------------
# Étape 1b : Fusionner avec le fichier biomart pour chaque gene
# -----------------------------------------------------------
def fuse_with_biomart(gene_dict):
    """
    Pour chaque gène, si un fichier {gene}_biomart.txt existe dans le dossier,
    fusionne son contenu avec le fichier {gene}_transcripts.txt en standardisant les colonnes.
    Le fichier biomart a pour header: "Gene name", "Transcript stable ID version", "Protein stable ID version"
    et sera renommé en: "Gene", "Transcript_ID", "PA".
    Le fichier final est sauvegardé sous le nom "{gene}_final_transcripts.txt".
    Les valeurs de la colonne Gene sont conservées telles qu'elles apparaissent.
    """
    import os
    import pandas as pd
    
    for gene in gene_dict.keys():
        folder = gene
        openprot_file = os.path.join(folder, f"{gene}_transcripts.txt")
        biomart_file = os.path.join(folder, f"{gene}_biomart.txt")
        dfs = []
        if os.path.exists(openprot_file):
            df_openprot = pd.read_csv(openprot_file, sep='\t')
            dfs.append(df_openprot)
        if os.path.exists(biomart_file):
            df_biomart = pd.read_csv(biomart_file, sep='\t')
            df_biomart = df_biomart.rename(columns={
                "Gene name": "Gene",
                "Transcript stable ID version": "Transcript_ID",
                "Protein stable ID version": "PA"
            })
            dfs.append(df_biomart)
        if dfs:
            df_final = pd.concat(dfs, ignore_index=True)
            final_file = os.path.join(folder, f"{gene}_final_transcripts.txt")
            df_final.to_csv(final_file, sep='\t', index=False)
        else:
            print(f"Aucun fichier trouvé pour le gène {gene}")

# -----------------------------------------------------------
# Étape 2 : Combinaison des fichiers finaux en un fichier unique pour la suite du pipeline
# -----------------------------------------------------------
def combine_gene_files(gene_dict, combined_filename="genes_list.tsv"):
    """
    Concatène les fichiers {gene}_final_transcripts.txt de chaque dossier en un seul fichier.
    
    :param gene_dict: dictionnaire des gènes utilisé précédemment.
    :param combined_filename: nom du fichier combiné qui sera utilisé pour la suite.
    """
    combined_rows = []
    header = "Gene\tTranscript_ID\tPA\n"
    for gene in gene_dict.keys():
        file_path = os.path.join(gene, f"{gene}_final_transcripts.txt")
        if os.path.exists(file_path):
            with open(file_path, "r") as infile:
                next(infile)  # ignorer le header
                for line in infile:
                    combined_rows.append(line)
    with open(combined_filename, "w") as outfile:
        outfile.write(header)
        for row in combined_rows:
            outfile.write(row)
    print("Fichier combiné créé:", combined_filename)
    return combined_filename

# -----------------------------------------------------------
# Étape 3 : Création du dictionnaire des gènes à partir du fichier combiné
# -----------------------------------------------------------
def creating_genes_dic_from_tsv(enriched_file):
    """
    Lit le fichier genes_list.tsv et retourne un dictionnaire {Gene: [list of Transcript_ID]}.
    """
    df = pd.read_csv(enriched_file, sep='\t')
    df.columns = df.columns.str.strip()
    gene_dic = df.groupby("Gene")["Transcript_ID"].apply(list).to_dict()
    return gene_dic

# -----------------------------------------------------------
# Étape 4 : Annotation et conversion (filtrage du fichier GTEx)
# -----------------------------------------------------------
def annotation_conversion(gene_dic, gtex_mapping_file, gtex_file):
    annotation_df = pd.read_csv(gtex_mapping_file, sep='\t', low_memory=False)
    gtex_to_tissue = dict(zip(annotation_df['SAMPID'], annotation_df['SMTS']))
    
    transcripts = []
    for tx_list in gene_dic.values():
        transcripts.extend(tx_list)
    transcripts_pattern = "|".join(t.replace(".", "\\.") for t in transcripts)
    print("Pattern de transcrits :", transcripts_pattern)
    
    os.makedirs("inter_output", exist_ok=True)
    awk_command = f"awk 'NR==1 || $1 ~ /^({transcripts_pattern})$/' {gtex_file} > inter_output/filtered_gtex_file.txt"
    subprocess.run(awk_command, shell=True)
    
    dat_gene = pd.read_csv("inter_output/filtered_gtex_file.txt", sep='\t')
    dat_gene.set_index('transcript_id', inplace=True)
    dat_gene_renamed = dat_gene.rename(columns=gtex_to_tissue)
    renamed_columns = set(gtex_to_tissue.values())
    dat_gene_tpm = dat_gene_renamed.loc[:, dat_gene_renamed.columns.isin(renamed_columns) | (dat_gene_renamed.columns == 'transcript_id')]
    dat_gene_tpm = dat_gene_renamed.loc[:, dat_gene_renamed.columns.isin(renamed_columns) | (dat_gene_renamed.columns == 'transcript_id')]
    
    non_numeric_cols = dat_gene_tpm.select_dtypes(include=['object']).columns
    numeric_cols = dat_gene_tpm.select_dtypes(include=['number']).columns
    
    df_merged = dat_gene_tpm[numeric_cols].T.groupby(level=0).mean().T
    for col in non_numeric_cols:
        series_col = dat_gene_tpm.loc[:, col]
        if isinstance(series_col, pd.DataFrame):
            series_col = series_col.iloc[:, 0]
        if col not in df_merged.columns:
            df_merged.insert(0, col, series_col)
    
    df_merged.to_csv("inter_output/filtered_gtex_file_merged.txt", sep='\t', index=False)
    return df_merged

# -----------------------------------------------------------
# Étape 5 : Clustering
# -----------------------------------------------------------
def cluster(merged_tissues, gene_key, directory):
    gene_key = gene_key.strip().split(':')[0]
    # Exclure la colonne 'transcript_id' pour garder uniquement les données numériques
    if 'transcript_id' in merged_tissues.columns:
        data_for_cluster = merged_tissues.drop(columns=['transcript_id'])
    else:
        data_for_cluster = merged_tissues.copy()
    
    # Transposer les données pour avoir les tissus en lignes
    transposed_data = data_for_cluster.T
    
    # Normalisation
    scaler = StandardScaler()
    normalized_data = scaler.fit_transform(transposed_data)
    
    # Clustering hiérarchique
    Z = linkage(normalized_data, method='ward')
    
    # --- Lignes pour générer le dendrogram et les figures (désactivées)
    # plt.figure(figsize=(10, 7))
    # dendrogram(Z, labels=transposed_data.index)
    # plt.title('Dendrogram for Hierarchical Clustering')
    # plt.xlabel('Tissues')
    # plt.ylabel('Distance')
    # plt.savefig(os.path.join(directory, 'dendrogram.png'))
    # plt.close()
    
    distances = Z[:, 2]
    max_d = np.percentile(distances, 75)
    clusters = fcluster(Z, max_d, criterion='distance')
    
    # Ajout de l'information cluster aux données transposées
    transposed_data['Cluster'] = clusters
    
    # --- Lignes pour générer la distribution des clusters (désactivées)
    # plt.figure(figsize=(10, 6))
    # sns.countplot(x='Cluster', data=transposed_data)
    # plt.title('Distribution des clusters')
    # plt.xlabel('Cluster')
    # plt.ylabel('Nombre de tissues')
    # clusters_distribution_path = os.path.join(directory, 'clusters_distribution.png')
    # plt.savefig(clusters_distribution_path)
    # plt.close()
    
    return transposed_data

# -----------------------------------------------------------
# Fonction principale
# -----------------------------------------------------------
import sys


def extract_transcript_ids(final_transcript_file, output_file="transcript_ids.txt"):
    """
    Extrait uniquement la colonne Transcript_ID du fichier final_transcripts.txt, sans inclure le header.
    
    :param final_transcript_file: Chemin vers le fichier final_transcripts.txt.
    :param output_file: Chemin du fichier de sortie où seront enregistrés les transcript_ID. Par défaut, "transcript_ids.txt".
    :return: Une liste contenant les transcript_ID.
    """
    transcript_ids = []
    with open(final_transcript_file, 'r') as f:
        # Ignorer le header
        next(f)
        for line in f:
            # Supposons que le fichier est au format "Gene\tTranscript_ID\tPA"
            columns = line.strip().split('\t')
            if len(columns) >= 2:
                transcript_ids.append(columns[1])
    
    # Optionnel : écrire les transcript_ID dans un fichier de sortie
    with open(output_file, 'w') as out:
        for tid in transcript_ids:
            out.write(tid + "\n")
    
    return transcript_ids



def main():
    # Fichiers de référence
    fasta_file = "references/human-openprot-2_1-refprots+altprots+isoforms-uniprot2022_06_01.fasta"
    tsv_file = "reference/openprot_database.tsv"
    gtex_mapping_file = "reference/annoted.txt"  # Mapping GTEx (SAMPID -> SMTS)
    gtex_file = "reference/GTEx_Analysis_2022-06-06_v10_RSEMv1.3.3_transcripts_tpm.txt"
    
    gene_dict = {
        "GAPDH": ["GAPDH", "GAPDHP"],
        "PPIA": ["PPIA", "PPIAP"],
        "PHB":  ["PHB1",  "PHBP"],
        "EEF1A1": ["EEF1A1", "EEF1A1P"],
        "HMGB1":  ["HMGB1", "HMGB1P"]
    }

    # Extraction des entrées et fusion avec biomart (si applicable)
    extract_best_entries(tsv_file, gene_dict)
    fuse_with_biomart(gene_dict)
    
    # Combinaison des fichiers finaux en un fichier global
    enriched_file = combine_gene_files(gene_dict, combined_filename="genes_list.tsv")
    
    try:
        gene_dic = creating_genes_dic_from_tsv(enriched_file)
        print("Dictionnaire des gènes:", gene_dic)
    except ValueError as e:
        print("❌ Erreur lors de la création du dictionnaire de gènes :", e)
        return
    
    # Si des arguments sont passés, on ne traite que ces gènes
    if len(sys.argv) > 1:
        selected_genes = sys.argv[1:]
    else:
        selected_genes = gene_dic.keys()
    
    combined_cluster_file = "clusterdata_combined.txt"
    if os.path.exists(combined_cluster_file):
        os.remove(combined_cluster_file)
    
    for gene_key in selected_genes:
        if gene_key not in gene_dic:
            print(f"Le gène {gene_key} n'est pas présent dans le dictionnaire.")
            continue

        transcripts = gene_dic[gene_key]
        gene_key_clean = gene_key.strip('>:')
        directory = f"./{gene_key_clean}/"
        os.makedirs(directory, exist_ok=True)
        
        merged_tissues = annotation_conversion({gene_key: transcripts}, gtex_mapping_file, gtex_file)
        transposed_data = cluster(merged_tissues, gene_key, directory)
        if transposed_data is not None:
            transposed_data.reset_index(inplace=True)
            transposed_data.insert(0, 'Gene', gene_key_clean)
            value_vars = [col for col in transposed_data.columns if col not in ['Gene', 'index', 'Cluster', 'transcript_id']]
            tidy_data = transposed_data.melt(id_vars=['Gene', 'index', 'Cluster'],
                                             value_vars=value_vars,
                                             var_name='Transcript',
                                             value_name='TPM')
            tidy_data.rename(columns={'index': 'Tissue'}, inplace=True)
            tidy_data.to_csv(combined_cluster_file, sep='\t', mode='a',
                               header=not os.path.exists(combined_cluster_file), index=False)
            gene_cluster_file = os.path.join(directory, f"{gene_key_clean}_clusterdata.txt")
            tidy_data.to_csv(gene_cluster_file, sep='\t', index=False)
    # Exemple d'appel de la fonction
    transcript_ids = extract_transcript_ids("GAPDH/GAPDH_final_transcripts.txt", output_file="GAPDH/transcript_ids.txt")
    print("Les Transcript_ID extraits sont :", transcript_ids)

if __name__ == "__main__":
    main()
