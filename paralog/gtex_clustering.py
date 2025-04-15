import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
import numpy as np
import os

from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster

# -----------------------------------------------------------
# Étape 1 : Tu obtiens le fichier initial BioMart du gène d'intérêt
# (Ce fichier contient les colonnes suivantes :
#  Gene stable ID version, Transcript stable ID version, 
#  Human paralogue gene stable ID, Human paralogue associated gene name, 
#  Human paralogue protein or transcript ID, Paralogue %id. target Human gene identical to query gene,
#  Gene name, Protein stable ID version)
# -----------------------------------------------------------

# -----------------------------------------------------------
# Étape 2 : Extraction des ENSG des paralogues
# Cette fonction lit le fichier BioMart (fichier 1) et crée un fichier texte (fichier 2)
# contenant la liste unique des "Human paralogue gene stable ID"
# -----------------------------------------------------------
def extract_paralogues_ENSG(biomart_file, output_file="paralogues_ENSG.txt"):
    """
    Extrait les ENSG des paralogues depuis le fichier BioMart (fichier 1)
    et les enregistre dans output_file (fichier 2).
    """
    df = pd.read_csv(biomart_file, sep="\t")
    df.columns = df.columns.str.strip()
    
    col = "Human paralogue gene stable ID"
    if col not in df.columns:
        print(f"⚠️ Colonne '{col}' introuvable dans le fichier {biomart_file}")
        return

    paralogue_ensg = df[col].dropna().unique()
    if len(paralogue_ensg) == 0:
        print("⚠️ Aucun paralogue trouvé.")
        return

    with open(output_file, "w") as f:
        for gene_id in paralogue_ensg:
            f.write(f"{gene_id}\n")
    print(f"✅ {len(paralogue_ensg)} paralogue(s) ENSG extrait(s) vers {output_file}")


# -----------------------------------------------------------
# Étape 3 (manuelle) : 
# Tu copies le contenu du fichier 'paralogues_ENSG.txt' dans BioMart afin d'obtenir 
# un fichier (fichier 3) contenant tous les transcrits (ENST) pour chaque ENSG.
# Ce fichier doit contenir au minimum les colonnes "Gene name" et "Transcript stable ID version".
# -----------------------------------------------------------

# -----------------------------------------------------------
# Étape 4 : Fusionner et dédupliquer les transcrits issus du fichier 1 et du fichier 3.
# On crée un fichier unique "genes_list.tsv" (fichier 4) avec les colonnes : 
# "Gene name" et "Transcript_ID" (qui sera au format "ENSTxxxx.y")
# -----------------------------------------------------------
def combine_transcripts(file1, file3, output_file="reference/genes_list.tsv"):
    """
    Combine les transcrits issus du fichier 3 et du fichier 1 pour obtenir une liste unique.
    
    Paramètres :
      - file3 : chemin vers le fichier BioMart contenant les transcrits du gène d'intérêt 
                (header attendu : Gene stable ID, Gene stable ID version, Transcript stable ID version, 
                               Gene name, Protein stable ID, Protein stable ID version)
      - file1 : chemin vers le fichier BioMart des paralogues 
                (header attendu : Gene stable ID version, Transcript stable ID version, Human paralogue gene stable ID,
                               Human paralogue associated gene name, Human paralogue protein or transcript ID,
                               Paralogue %id. target Human gene identical to query gene, Gene name, Protein stable ID version)
      - output_file : chemin de sortie pour le fichier combiné (fichier 4)
    
    Le fichier de sortie aura deux colonnes :
      - Name : le nom du gène (tel quel)
      - Transcript_ID : le Transcript stable ID version (ex. ENST00000278302.5)
    """
    import pandas as pd
    
    # Lecture des deux fichiers
    df3 = pd.read_csv(file3, sep="\t")
    df1 = pd.read_csv(file1, sep="\t")
    
    # Nettoyage des noms de colonnes
    df3.columns = df3.columns.str.strip()
    df1.columns = df1.columns.str.strip()
    
    # Sélectionner les colonnes "Gene name" et "Transcript stable ID version"
    try:
        df3_sub = df3[["Gene name", "Transcript stable ID version"]]
    except KeyError as e:
        raise ValueError(f"Le fichier 3 doit contenir les colonnes 'Gene name' et 'Transcript stable ID version'. Erreur: {e}")
    
    try:
        df1_sub = df1[["Gene name", "Transcript stable ID version"]]
    except KeyError as e:
        raise ValueError(f"Le fichier 1 doit contenir les colonnes 'Gene name' et 'Transcript stable ID version'. Erreur: {e}")
    
    # Concaténer et supprimer les doublons
    df_combined = pd.concat([df3_sub, df1_sub], ignore_index=True).drop_duplicates()
    
    # On ne modifie pas le "Gene name", on le conserve tel quel
    df_out = df_combined.copy()
    df_out.columns = ["Gene", "Transcript_ID"]
    
    # Sauvegarder le résultat
    df_out.to_csv(output_file, sep="\t", index=False)
    print(f"✅ Fichier combiné généré : {output_file} ({df_out.shape[0]} transcrits uniques)")
    return output_file


# Exemple d'utilisation :
# combine_transcripts("path/to/fichier1.tsv", "path/to/fichier3.tsv", output_file="reference/genes_list.tsv")



# -----------------------------------------------------------
# Étape 5 : Utilisation de genes_list.tsv pour filtrer le fichier GTEx et produire les TPM
# (Les fonctions suivantes restent inchangées pour filtrer GTEx, générer la heatmap et le clustering)
# -----------------------------------------------------------
def creating_genes_dic_from_tsv(enriched_file):
    """
    Lit le fichier genes_list.tsv et retourne un dictionnaire {gene: [list of transcript IDs]}.
    """
    df = pd.read_csv(enriched_file, sep='\t')
    df.columns = df.columns.str.strip()
    
    # Ici, on se base uniquement sur les colonnes "Gene" et "Transcript_ID" du fichier préparé
    gene_dic = df.groupby("Gene")["Transcript_ID"].apply(list).to_dict()
    return gene_dic

def annotation_conversion(gene_dic, gtex_mapping_file, gtex_file):
    annotation_df = pd.read_csv(gtex_mapping_file, sep='\t', low_memory=False)
    # Création du dictionnaire de mapping : SAMPID -> SMTS (nom du tissu)
    gtex_to_tissue = dict(zip(annotation_df['SAMPID'], annotation_df['SMTS']))
    
    # Récupérer la liste de tous les identifiants de transcrits
    transcripts = []
    for tx_list in gene_dic.values():
        transcripts.extend(tx_list)
    transcripts_pattern = "|".join(t.replace(".", "\\.") for t in transcripts)
    print(transcripts_pattern)
    # S'assurer que le dossier inter_output existe
    os.makedirs("inter_output", exist_ok=True)
    
    # Filtrer le fichier GTEx pour ne garder que les lignes correspondant aux transcrits d'intérêt
    awk_command = f"awk 'NR==1 || $1 ~ /^({transcripts_pattern})$/' {gtex_file} > inter_output/filtered_gtex_file.txt"
    subprocess.run(awk_command, shell=True)
    
    dat_gene = pd.read_csv("inter_output/filtered_gtex_file.txt", sep='\t')
    dat_gene.set_index('transcript_id', inplace=True)

    dat_gene_renamed = dat_gene.rename(columns=gtex_to_tissue)
    renamed_columns = set(gtex_to_tissue.values())
    dat_gene_tpm = dat_gene_renamed.loc[:, dat_gene_renamed.columns.isin(renamed_columns) | (dat_gene_renamed.columns == 'transcript_id')]
   
    
    # Sélectionner les colonnes numériques (TPM) et la colonne 'transcript_id'
    dat_gene_tpm = dat_gene_renamed.loc[:, dat_gene_renamed.columns.isin(renamed_columns) | (dat_gene_renamed.columns == 'transcript_id')]
    
    # Séparer les colonnes non numériques et numériques
    non_numeric_cols = dat_gene_tpm.select_dtypes(include=['object']).columns
    numeric_cols = dat_gene_tpm.select_dtypes(include=['number']).columns
    
    # Agréger les données numériques (moyenne par tissu)

    df_merged = dat_gene_tpm[numeric_cols].T.groupby(level=0).mean().T
    
    # Pour chaque colonne non numérique, on insère seulement si elle n'existe pas déjà dans df_merged.
    for col in non_numeric_cols:
        series_col = dat_gene_tpm.loc[:, col]
        # Si plusieurs colonnes portent le même nom, prendre la première
        if isinstance(series_col, pd.DataFrame):
            series_col = series_col.iloc[:, 0]
        if col not in df_merged.columns:
            df_merged.insert(0, col, series_col)
        else:
            # Si la colonne existe déjà, on peut l'ignorer
            pass
    
    df_merged.to_csv("inter_output/filtered_gtex_file_merged.txt", sep='\t', index=False)
    return df_merged

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
    
    plt.figure(figsize=(10, 7))
    dendrogram(Z, labels=transposed_data.index)
    plt.title('Dendrogram for Hierarchical Clustering')
    plt.xlabel('Tissues')
    plt.ylabel('Distance')
    plt.savefig(os.path.join(directory, 'dendrogram.png'))
    plt.close()
    
    distances = Z[:, 2]
    max_d = np.percentile(distances, 75)
    clusters = fcluster(Z, max_d, criterion='distance')
    
    # Ajout de l'information cluster aux données transposées
    transposed_data['Cluster'] = clusters
    
    plt.figure(figsize=(10, 6))
    sns.countplot(x='Cluster', data=transposed_data)
    plt.title('Distribution des clusters')
    plt.xlabel('Cluster')
    plt.ylabel('Nombre de tissues')
    clusters_distribution_path = os.path.join(directory, 'clusters_distribution.png')
    plt.savefig(clusters_distribution_path)
    plt.close()
    
    return transposed_data



def main():
    # Fichier initial (fichier 1) obtenu depuis BioMart
    biomart_file = "reference/trim27.tsv"
    ENST_paralogs = "reference/mart_export.txt"
    # Étape 2 : Extraction des ENSG des paralogues (fichier 2)
    extract_paralogues_ENSG(biomart_file, "reference/paralogues_ENSG.txt")
    print("📌 Copiez le contenu de 'reference/paralogues_ENSG.txt' dans BioMart pour obtenir les transcrits (fichier 3).")
    
    # Étape 4 : Préparation du fichier standardisé (fichier 4)
    combine_transcripts(biomart_file, ENST_paralogs, output_file="reference/genes_list.tsv")
    
    # Maintenant, on utilise ce fichier pour la suite du pipeline GTEx (étape 5)
    enriched_file = "reference/genes_list.tsv"
    gtex_mapping_file = 'reference/annoted.txt'  # Mapping GTEx (SAMPID -> SMTS)
    gtex_file = "reference/GTEx_Analysis_2022-06-06_v10_RSEMv1.3.3_transcripts_tpm.txt"

    try:
        gene_dic = creating_genes_dic_from_tsv(enriched_file)
    except ValueError as e:
        print("❌ Erreur lors de la création du dictionnaire de gènes :", e)
        return
    
    all_transposed_data = pd.DataFrame()
    if os.path.exists("clusterdata_combined.txt"):
        os.remove("clusterdata_combined.txt")
    for gene_key, transcripts in gene_dic.items():
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
                tidy_data.to_csv("clusterdata_combined.txt", sep='\t', mode='a',
                                header=not os.path.exists("clusterdata_combined.txt"), index=False)
                all_transposed_data = pd.concat([all_transposed_data, tidy_data], ignore_index=True)

if __name__ == "__main__":
    main()