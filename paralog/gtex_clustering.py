import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
import numpy as np
import os

from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster

# =============================================================================
#                        SECTION 1: INITIAL SETUP
# =============================================================================
# This script processes BioMart data, extracts paralogue identifiers,
# combines transcript lists, filters the GTEx file to obtain TPM values,
# performs hierarchical clustering, and generates plots.
#
# Input file (file 1) obtained from BioMart contains the following columns:
#   - Gene stable ID version
#   - Transcript stable ID version
#   - Human paralogue gene stable ID
#   - Human paralogue associated gene name
#   - Human paralogue protein or transcript ID
#   - Paralogue %id. target Human gene identical to query gene
#   - Gene name
#   - Protein stable ID version
# =============================================================================

# =============================================================================
# STEP 2: EXTRACT ENSG OF PARALOGUES
# =============================================================================
def extract_paralogues_ENSG(biomart_file, output_file="paralogues_ENSG.txt"):
    """
    Extracts the ENSG identifiers of paralogues from the BioMart file (file 1)
    and writes them to output_file (file 2).
    """
    df = pd.read_csv(biomart_file, sep="\t")
    df.columns = df.columns.str.strip()
    
    col = "Human paralogue gene stable ID"
    if col not in df.columns:
        print(f" Column '{col}' not found in file {biomart_file}")
        return

    paralogue_ensg = df[col].dropna().unique()
    if len(paralogue_ensg) == 0:
        print(" No paralogue found.")
        return

    with open(output_file, "w") as f:
        for gene_id in paralogue_ensg:
            f.write(f"{gene_id}\n")
    print(f" {len(paralogue_ensg)} paralogue(s) ENSG extracted to {output_file}")


# =============================================================================
# STEP 3 (MANUAL STEP):
# Copy the content of 'paralogues_ENSG.txt' into BioMart to obtain a file
# (file 3) containing all the transcripts (ENST) for each ENSG.
# This file must include at least the columns "Gene name" and 
# "Transcript stable ID version".
# =============================================================================

# =============================================================================
# STEP 4: COMBINE AND DEDUPLICATE TRANSCRIPTS FROM FILE 1 AND FILE 3
# =============================================================================
def combine_transcripts(file1, file3, output_file="reference/genes_list.tsv"):
    """
    Combines transcripts from file 3 and file 1 to obtain a unique list.
    
    Parameters:
      - file3: path to the BioMart file containing transcripts for the gene of interest 
               (expected header: Gene stable ID, Gene stable ID version, Transcript stable ID version, 
                                Gene name, Protein stable ID, Protein stable ID version)
      - file1: path to the BioMart paralogues file 
               (expected header: Gene stable ID version, Transcript stable ID version, Human paralogue gene stable ID,
                                Human paralogue associated gene name, Human paralogue protein or transcript ID,
                                Paralogue %id. target Human gene identical to query gene, Gene name, Protein stable ID version)
      - output_file: path for the combined output file (file 4)
    
    The output file will have two columns:
      - Gene: the gene name (unchanged)
      - Transcript_ID: the Transcript stable ID version (e.g., ENST00000278302.5)
    """
    # Read both files
    df3 = pd.read_csv(file3, sep="\t")
    df1 = pd.read_csv(file1, sep="\t")
    
    # Clean column names
    df3.columns = df3.columns.str.strip()
    df1.columns = df1.columns.str.strip()
    
    # Select the columns "Gene name" and "Transcript stable ID version"
    try:
        df3_sub = df3[["Gene name", "Transcript stable ID version"]]
    except KeyError as e:
        raise ValueError(f"File 3 must contain columns 'Gene name' and 'Transcript stable ID version'. Error: {e}")
    
    try:
        df1_sub = df1[["Gene name", "Transcript stable ID version"]]
    except KeyError as e:
        raise ValueError(f"File 1 must contain columns 'Gene name' and 'Transcript stable ID version'. Error: {e}")
    
    # Concatenate and remove duplicates
    df_combined = pd.concat([df3_sub, df1_sub], ignore_index=True).drop_duplicates()
    
    # Do not modify "Gene name"; keep it as is.
    df_out = df_combined.copy()
    df_out.columns = ["Gene", "Transcript_ID"]
    
    # Save the result
    df_out.to_csv(output_file, sep="\t", index=False)
    print(f" Combined file generated: {output_file} ({df_out.shape[0]} unique transcripts)")
    return output_file


# =============================================================================
# STEP 5: USE genes_list.tsv TO FILTER THE GTEx FILE AND GENERATE TPM VALUES
# =============================================================================
def creating_genes_dic_from_tsv(enriched_file):
    """
    Reads the genes_list.tsv file and returns a dictionary {gene: [list of transcript IDs]}.
    """
    df = pd.read_csv(enriched_file, sep='\t')
    df.columns = df.columns.str.strip()
    
    # Build dictionary based solely on "Gene" and "Transcript_ID" columns.
    gene_dic = df.groupby("Gene")["Transcript_ID"].apply(list).to_dict()
    return gene_dic

def annotation_conversion(gene_dic, gtex_mapping_file, gtex_file):
    """
    Filters the GTEx file for transcripts of interest and aggregates TPM values.
    
    Parameters:
      - gene_dic: dictionary {gene: [list of transcript IDs]}
      - gtex_mapping_file: GTEx mapping file (SAMPID -> SMTS)
      - gtex_file: GTEx transcripts TPM file
      
    Returns:
      - A DataFrame with merged TPM data.
    """
    annotation_df = pd.read_csv(gtex_mapping_file, sep='\t', low_memory=False)
    # Create mapping dictionary: SAMPID -> SMTS (tissue name)
    gtex_to_tissue = dict(zip(annotation_df['SAMPID'], annotation_df['SMTS']))
    
    # Retrieve a list of all transcript IDs
    transcripts = []
    for tx_list in gene_dic.values():
        transcripts.extend(tx_list)
    transcripts_pattern = "|".join(t.replace(".", "\\.") for t in transcripts)
    print(transcripts_pattern)
    
    # Ensure the 'inter_output' directory exists
    os.makedirs("inter_output", exist_ok=True)
    
    # Filter the GTEx file: keep rows with headers or matching transcripts
    awk_command = f"awk 'NR==1 || $1 ~ /^({transcripts_pattern})$/' {gtex_file} > inter_output/filtered_gtex_file.txt"
    subprocess.run(awk_command, shell=True)
    
    dat_gene = pd.read_csv("inter_output/filtered_gtex_file.txt", sep='\t')
    dat_gene.set_index('transcript_id', inplace=True)
    
    # Rename columns based on the mapping (SAMPID -> SMTS)
    dat_gene_renamed = dat_gene.rename(columns=gtex_to_tissue)
    renamed_columns = set(gtex_to_tissue.values())
    dat_gene_tpm = dat_gene_renamed.loc[:, dat_gene_renamed.columns.isin(renamed_columns) | (dat_gene_renamed.columns == 'transcript_id')]
    
    # Select TPM (numeric) columns and the 'transcript_id' column
    dat_gene_tpm = dat_gene_renamed.loc[:, dat_gene_renamed.columns.isin(renamed_columns) | (dat_gene_renamed.columns == 'transcript_id')]
    
    # Separate non-numeric and numeric columns
    non_numeric_cols = dat_gene_tpm.select_dtypes(include=['object']).columns
    numeric_cols = dat_gene_tpm.select_dtypes(include=['number']).columns
    
    # Aggregate numeric data (compute mean per tissue)
    df_merged = dat_gene_tpm[numeric_cols].T.groupby(level=0).mean().T
    
    # For each non-numeric column, add it if not already present in df_merged.
    for col in non_numeric_cols:
        series_col = dat_gene_tpm.loc[:, col]
        if isinstance(series_col, pd.DataFrame):
            series_col = series_col.iloc[:, 0]
        if col not in df_merged.columns:
            df_merged.insert(0, col, series_col)
    
    df_merged.to_csv("inter_output/filtered_gtex_file_merged.txt", sep='\t', index=False)
    return df_merged

def cluster(merged_tissues, gene_key, directory):
    """
    Performs hierarchical clustering on TPM data per tissue.
    
    Parameters:
      - merged_tissues: DataFrame with merged TPM data.
      - gene_key: gene name string.
      - directory: path to save output plots.
      
    Returns:
      - A transposed DataFrame with an added 'Cluster' column.
    """
    gene_key = gene_key.strip().split(':')[0]
    # Exclude the 'transcript_id' column for numeric data
    if 'transcript_id' in merged_tissues.columns:
        data_for_cluster = merged_tissues.drop(columns=['transcript_id'])
    else:
        data_for_cluster = merged_tissues.copy()
    
    # Transpose data so that tissues become rows
    transposed_data = data_for_cluster.T
    
    # Normalize the data
    scaler = StandardScaler()
    normalized_data = scaler.fit_transform(transposed_data)
    
    # Perform hierarchical clustering using Ward's method
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
    
    # Add cluster information to the transposed data
    transposed_data['Cluster'] = clusters
    
    plt.figure(figsize=(10, 6))
    sns.countplot(x='Cluster', data=transposed_data)
    plt.title('Cluster Distribution')
    plt.xlabel('Cluster')
    plt.ylabel('Number of tissues')
    clusters_distribution_path = os.path.join(directory, 'clusters_distribution.png')
    plt.savefig(clusters_distribution_path)
    plt.close()
    
    return transposed_data

# =============================================================================
# STEP 6: MAIN PIPELINE EXECUTION
# =============================================================================
def main():
    # File paths (adjust as necessary)
    biomart_file = "reference/trim27.tsv"         # File 1: Initial BioMart file
    ENST_paralogs = "reference/mart_export.txt"     # File with paralog transcripts (file from BioMart)
    
    # STEP 2: Extract paralogue ENSG IDs (file 2)
    extract_paralogues_ENSG(biomart_file, "reference/paralogues_ENSG.txt")
    print(" Copy the content of 'reference/paralogues_ENSG.txt' into BioMart to obtain transcripts (file 3).")
    
    # STEP 4: Prepare the standardized transcript file (file 4)
    combine_transcripts(biomart_file, ENST_paralogs, output_file="reference/genes_list.tsv")
    
    # Use the prepared file for the downstream GTEx pipeline (STEP 5)
    enriched_file = "reference/genes_list.tsv"
    gtex_mapping_file = 'reference/annoted.txt'   # GTEx mapping file (SAMPID -> SMTS)
    gtex_file = "reference/GTEx_Analysis_2022-06-06_v10_RSEMv1.3.3_transcripts_tpm.txt"
    
    try:
        gene_dic = creating_genes_dic_from_tsv(enriched_file)
    except ValueError as e:
        print(" Error creating gene dictionary:", e)
        return

    all_transposed_data = pd.DataFrame()
    # Remove old combined file if it exists
    if os.path.exists("clusterdata_combined.txt"):
        os.remove("clusterdata_combined.txt")
    
    # Process each gene in the dictionary
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
