# =============================================================================
#                           LOAD LIBRARIES
# =============================================================================
library(dplyr)
library(tidyr)
library(flextable)
library(ggplot2)
library(pagedown)
library(patchwork)
library(stringr)
library(parallel)
library(rlang)
library(officer)
library(readr)
library(png)
library(grid)

# =============================================================================
#                SECTION 1: READ AND PREPARE EXPRESSION DATA
# =============================================================================

# Input file: clusterdata_combined.txt (columns: Gene, Tissue, Cluster, Transcript, TPM) (from the python script output)

# 1. Read the input files
results_list <- read_tsv("PHB/clusterdata_combined.txt", show_col_types = FALSE)
pseudogenes   <- read_tsv("PHB/PHB_pseudogen.txt", show_col_types = FALSE)
info          <- read_tsv("PHB/PHB_information.txt", show_col_types = FALSE)

# 2. Merge info and pseudogenes to obtain MANE select data.
# A full join keeps canonical information from 'info' and also includes pseudogenes (with extra columns as NA)
mane_select <- full_join(info, pseudogenes,
                         by = c("Gene name", "Gene stable ID version", "Transcript stable ID version", "Protein stable ID version"))

# 3. Left join results with MANE select using the transcript column.
results_list_enriched <- results_list %>%
  left_join(
    mane_select %>% select(`Transcript stable ID version`, 
                           `Protein stable ID version`, 
                           `RefSeq match transcript (MANE Select)`, 
                           `UniProtKB isoform ID`,
                           `UniProtKB/TrEMBL ID`),
    by = c("Transcript" = "Transcript stable ID version")
  ) %>%
  # Rename new columns for clarity
  rename(
    ProteinID = `Protein stable ID version`,
    MANE_Transcript = `RefSeq match transcript (MANE Select)`,
    UniProt_ID_isoform = `UniProtKB isoform ID`,
    UniProt_ID_trembl = `UniProtKB/TrEMBL ID`
  ) %>%
  mutate(
    UniProt_ID = coalesce(UniProt_ID_isoform, UniProt_ID_trembl)
  ) %>%
  select(-UniProt_ID_isoform, -UniProt_ID_trembl) %>%
  mutate(is_canonical = !is.na(MANE_Transcript) & MANE_Transcript != "")

# Visualize the enriched result
print(results_list_enriched)

# =============================================================================
#                SECTION 2: CALCULATE EXPRESSION RATIOS
# =============================================================================

results_list_enriched <- results_list_enriched %>%
  group_by(Gene, Tissue, Cluster) %>%
  mutate(ref_tpm = TPM[is_canonical == TRUE][1],
         ratio = TPM / ref_tpm) %>%
  ungroup()

agg_df <- results_list_enriched %>%
  group_by(Gene, Transcript, Cluster, is_canonical) %>%
  summarise(
    count      = n(),
    ratios     = if (!first(is_canonical)) paste(round(ratio, 3), collapse = "; ") else NA_character_,
    mean_ratio = if (!first(is_canonical)) round(mean(ratio), 3) else NA_real_,
    .groups    = "drop"
  ) %>%
  mutate(
    cell = ifelse(is_canonical,
                  paste0("n = ", count),
                  ifelse(!is.na(ratios), paste0(ratios, " (", mean_ratio, ")"), NA_character_)
    )
  )

# =============================================================================
#                SECTION 3: FORMAT SUMMARY TABLE FOR ISOFORMS
# =============================================================================
    
# Pivot to wide format (one column per cluster)
wide_df <- agg_df %>%
  mutate(cluster_col = paste0("cluster_", Cluster)) %>%
  select(Gene, Transcript, cluster_col, cell, is_canonical) %>%
  pivot_wider(names_from = cluster_col, values_from = cell) %>%
  arrange(Gene, desc(is_canonical))

# Append "(REF)" to the reference transcript(s)
wide_df <- wide_df %>%
  mutate(Transcript = ifelse(Transcript %in% unique(results_list_enriched$Transcript[results_list_enriched$is_canonical]),
                             paste0(Transcript, " (REF)"),
                             Transcript))

# Create a flextable, export as HTML and then convert to PDF
ft <- flextable(wide_df)
ft <- autofit(ft)
ft
save_as_html(ft, path = "expression_isoform.html")
pagedown::chrome_print(input = "expression_isoform.html", output = "expression_isoform.pdf")

# =============================================================================
#           SECTION 4: PREPARE FINAL TABLE (CLUSTER MEANS)
# =============================================================================
    
mean_df <- agg_df %>%
  select(Gene, Transcript, Cluster, mean_ratio) %>%
  pivot_wider(names_from = Cluster, values_from = mean_ratio, names_prefix = "cluster_") %>%
  arrange(Gene, desc(Transcript %in% unique(results_list_enriched$Transcript[results_list_enriched$is_canonical]))) %>%
  mutate(Transcript = ifelse(Transcript %in% unique(results_list_enriched$Transcript[results_list_enriched$is_canonical]),
                             paste0(Transcript, " (REF)"),
                             Transcript))

mean_df_long <- mean_df %>%
  pivot_longer(
    cols = starts_with("cluster_"),
    names_to = "Cluster",
    values_to = "mean_ratio"
  )

# Replace values lower than the 75th quantile for each (Gene, Cluster) group with NA
mean_df_long <- mean_df_long %>%
  group_by(Gene, Cluster) %>%
  mutate(mean_ratio = ifelse(mean_ratio > quantile(mean_ratio, probs = 0.75, na.rm = TRUE),
                             mean_ratio, NA)) %>%
  ungroup()

# Optionally, pivot back to wide format
mean_df_filtered_wide <- mean_df_long %>%
  pivot_wider(
    names_from = Cluster,
    values_from = mean_ratio
  ) %>%
  filter(str_detect(Transcript, fixed(" (REF)")) | !if_all(starts_with("cluster_"), is.na))

# =============================================================================
#        SECTION 5: LOCALIZATION ANALYSIS USING COSINE SIMILARITY
# =============================================================================
    
# Load additional libraries (if not already loaded)
library(dplyr)
library(tidyr)
library(stringr)

## Part A: Prepare mean_df_filtered_wide by extracting transcript IDs and flagging references.
mean_df_filtered_wide <- mean_df_filtered_wide %>%
  mutate(TranscriptID = str_remove(Transcript, fixed(" (REF)")),
         is_ref = str_detect(Transcript, fixed(" (REF)")))

## Part B: Build the localization matrix.
# Read and parse the prediction file.
pred_lines <- readLines("isoforms_localisation.txt")
pred_raw <- tibble(raw_line = pred_lines) %>%
  filter(raw_line != "") %>%
  mutate(
    Transcript = sub("^(\\S+).*", "\\1", raw_line),
    details = ifelse(str_detect(raw_line, "details"),
                     sub(".*details\\s*", "", raw_line),
                     NA)
  )

# Split the "localization: score" pairs.
pred_details <- pred_raw %>%
  filter(!is.na(details)) %>%
  separate_rows(details, sep = ",\\s*") %>%
  separate(details, into = c("Localization", "Score"), sep = ":\\s*", convert = TRUE)

# Aggregate scores by Transcript and Localization.
agg_predictions <- pred_details %>%
  mutate(Score = as.numeric(trimws(Score))) %>%
  group_by(Transcript, Localization) %>%
  summarise(TotalScore = sum(Score, na.rm = TRUE), .groups = "drop") %>%
  group_by(Transcript) %>%
  mutate(MaxScore = max(TotalScore, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(FilteredScore = TotalScore) %>%  
  group_by(Transcript) %>%
  mutate(TotalFiltered = sum(FilteredScore, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Percentage = ifelse(!is.na(FilteredScore) & TotalFiltered > 0,
                             FilteredScore / TotalFiltered * 100, NA))

# Build the localization matrix (one vector per Transcript).
localization_matrix <- agg_predictions %>%
  select(Transcript, Localization, Percentage) %>%
  pivot_wider(names_from = Localization, values_from = Percentage, 
              values_fill = list(Percentage = 0)) %>%
  mutate(TranscriptID = Transcript)

## Part C: Calculate cosine similarity.
cosine_similarity <- function(vec1, vec2) {
  dot_product <- sum(vec1 * vec2)
  norm1 <- sqrt(sum(vec1^2))
  norm2 <- sqrt(sum(vec2^2))
  return(dot_product / (norm1 * norm2))
}

# Determine the localization columns to compare.
localization_cols <- setdiff(names(localization_matrix), c("Transcript", "TranscriptID"))

# Replace any NA with 0.
localization_matrix <- localization_matrix %>%
  mutate(across(all_of(localization_cols), ~ replace_na(.x, 0)))

# Join gene and reference information from mean_df_filtered_wide.
localization_matrix <- localization_matrix %>%
  left_join(mean_df_filtered_wide %>% select(Gene, TranscriptID, is_ref),
            by = "TranscriptID")

# For each gene, calculate cosine similarity by comparing each transcript to the reference.
similarity_by_gene <- localization_matrix %>%
  group_by(Gene) %>%
  mutate(ref_vector = list({
    ref_row <- filter(cur_data(), is_ref)
    if(nrow(ref_row) > 0) {
      as.numeric(ref_row[1, localization_cols])
    } else {
      rep(NA, length(localization_cols))
    }
  })) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(CosineSimilarity = if(all(is.na(unlist(ref_vector)))) NA else
           cosine_similarity(c_across(all_of(localization_cols)), unlist(ref_vector))) %>%
  ungroup()

## Part D: Integrate the cosine similarity into mean_df_filtered_wide.
final_df <- left_join(mean_df_filtered_wide, 
                      similarity_by_gene %>% select(TranscriptID, CosineSimilarity),
                      by = "TranscriptID")

View(final_df)

# =============================================================================
#                SECTION 6: TRANSCRIPT FUNCTION ANALYSIS
#                     (Processing the protnlm.txt file)
# =============================================================================
         
# Read the text file using an alternate name.
protnlm_lines <- readLines("protnlm.txt")

# Identify the transcript block start lines (lines starting with ">").
transcript_indices <- grep("^>", protnlm_lines)

result_list <- list()
reference_functions <- list()

for (i in seq_along(transcript_indices)) {
  start_idx <- transcript_indices[i]
  end_idx <- if (i < length(transcript_indices)) transcript_indices[i + 1] - 1 else length(protnlm_lines)
  block <- protnlm_lines[start_idx:end_idx]
  transcript_id <- sub("^>", "", block[1])
  gene_name <- sub("^ENST[0-9]+\\.[0-9]+\\|", "", block[1])
  
  func_scores <- list()
  if (length(block) > 1) {
    for (j in seq(2, length(block), by = 2)) {
      if (j + 1 <= length(block)) {
        func_name <- block[j]
        score <- as.numeric(block[j + 1])
        func_scores[[func_name]] <- score
      }
    }
  }
  
  # Consider the first transcript for each gene as the reference.
  if (!(gene_name %in% names(reference_functions))) {
    if (length(names(func_scores)) > 0) {
      ref_func <- if (any(grepl("TRIM5", names(func_scores)))) {
        names(func_scores)[grepl("TRIM5", names(func_scores))][1]
      } else {
        head(names(func_scores), 1)
      }
    } else {
      ref_func <- NA
    }
    reference_functions[[gene_name]] <- list(
      Fonction1 = ref_func,
      Score1 = if (!is.na(ref_func)) func_scores[[ref_func]] else NA
    )
  }
  
  # Compare the current function with the reference.
  if (length(names(func_scores)) > 0) {
    current_func <- head(names(func_scores), 1)
    is_match <- identical(current_func, reference_functions[[gene_name]]$Fonction1)
    current_score <- func_scores[[current_func]]
  } else {
    current_func <- NA
    is_match <- NA
    current_score <- NA
  }
  
  result_list[[i]] <- data.frame(
    GeneName        = gene_name,
    Transcript_id   = transcript_id,
    Fonction1       = current_func,
    Score1          = current_score,
    MatchReference  = is_match,
    stringsAsFactors = FALSE
  )
}

final_results_df <- do.call(rbind, result_list)
print(final_results_df)

# Extract only the relevant results for each gene.
final_results <- list()
unique_genes <- unique(final_results_df$GeneName)
for (gene in unique_genes) {
  gene_data <- final_results_df[final_results_df$GeneName == gene, ]
  gene_data <- gene_data[!is.na(gene_data$Fonction1), ]
  final_results[[gene]] <- gene_data
}
final_results_df <- do.call(rbind, final_results)
print(final_results_df)

# Clean the TranscriptID column.
final_results_df$TranscriptID <- sapply(strsplit(final_results_df$Transcript_id, "\\|"), `[`, 1)

# Merge with final_df but only select available columns.
merged_df <- merge(final_df, 
                   final_results_df[, c("TranscriptID", "Fonction1", "Score1")], 
                   by = "TranscriptID", 
                   all.x = TRUE)

# Remove unused columns.
merged_df <- subset(merged_df, select = -c(is_ref, TranscriptID)) 

# Create a flextable from the merged results.
merged_df <- merged_df %>%
  arrange(Gene, desc(grepl("\\(REF\\)", Transcript)))
ft <- flextable(merged_df)
ft <- autofit(ft)
save_as_html(ft, path = "sum_isoform.html")
html_content <- readLines("sum_isoform.html")
css <- "<style>@page { size: landscape; }</style>"
html_content <- sub("(<head.*?>)", paste0("\\1\n", css), html_content)
writeLines(html_content, "sum_isoform.html")
pagedown::chrome_print(input = "sum_isoform.html", output = "sum_isoform.pdf")

# =============================================================================
#                SECTION 7: GENERATE HEATMAP
# =============================================================================
      
# 1) Mark pseudogenes (adapt according to your columns)
results_list_enriched <- results_list_enriched %>%
  mutate(is_pseudogene = if_else(grepl("P\\d+$", Gene), TRUE, FALSE))

# 2) Identify canonical UniProt entries
canonical_uniprots <- results_list_enriched %>%
  filter(is_canonical, !is.na(UniProt_ID)) %>%
  pull(UniProt_ID) %>%
  unique()

# 3) Create a group label for each transcript/gene.
results_list_enriched <- results_list_enriched %>%
  mutate(group_label = case_when(
    is_pseudogene ~ "Pseudogenes",
    !is.na(UniProt_ID) & UniProt_ID %in% canonical_uniprots ~ paste0("Canonical: ", UniProt_ID),
    !is.na(UniProt_ID) ~ paste0("Protein: ", UniProt_ID),
    TRUE ~ "Non-coding transcripts"
  ))

# 4) List groups (canonical, protein, non-coding, pseudogenes)
protein_groups <- results_list_enriched %>%
  group_by(group_label) %>%
  summarise(transcripts = list(unique(Transcript)), .groups = "drop")

plots_list <- list()

# 5) Loop through each group label.
for (i in seq_len(nrow(protein_groups))) {
  label <- protein_groups$group_label[i]
  transcripts <- protein_groups$transcripts[[i]]
  
  # Subset data for this group.
  grp_data <- results_list_enriched %>%
    filter(Transcript %in% transcripts)
  
  # A) Aggregation: group by Gene for pseudogenes, else by Transcript.
  if (label == "Pseudogenes") {
    heatmap_df <- grp_data %>%
      group_by(Gene, Tissue) %>%
      summarise(avg_TPM = mean(TPM, na.rm = TRUE), .groups = "drop")
    y_var <- "Gene"
    y_label <- "Gene"
  } else {
    heatmap_df <- grp_data %>%
      group_by(Transcript, Tissue) %>%
      summarise(avg_TPM = mean(TPM, na.rm = TRUE), .groups = "drop")
    y_var <- "Transcript"
    y_label <- "Transcript"
  }
  
  # B) Min-max normalization.
  global_min <- min(heatmap_df$avg_TPM, na.rm = TRUE)
  global_max <- max(heatmap_df$avg_TPM, na.rm = TRUE)
  heatmap_df <- if (global_max - global_min == 0) {
    heatmap_df %>% mutate(norm_TPM = 0)
  } else {
    heatmap_df %>% mutate(norm_TPM = (avg_TPM - global_min) / (global_max - global_min))
  }
  
  # C) Reorder the y-axis: for non-pseudogenes, put canonical transcripts first.
  if (label != "Pseudogenes") {
    heatmap_df <- heatmap_df %>%
      mutate(is_can_transcript = Transcript %in% results_list_enriched$Transcript[results_list_enriched$is_canonical])
    can_levels <- unique(heatmap_df[[y_var]][heatmap_df$is_can_transcript])
    noncan_levels <- unique(heatmap_df[[y_var]][!heatmap_df$is_can_transcript])
    heatmap_df[[y_var]] <- factor(heatmap_df[[y_var]], levels = c(can_levels, noncan_levels))
  } else {
    heatmap_df <- heatmap_df %>%
      mutate(!!sym(y_var) := factor(!!sym(y_var), levels = sort(unique(!!sym(y_var)))))
  }
  
  # D) Create the main heatmap.
  heatmap_plot <- ggplot(heatmap_df, aes(x = Tissue, y = !!sym(y_var))) +
    geom_tile(aes(fill = norm_TPM), color = "white") +
    scale_fill_viridis_c(
      option = "magma",
      na.value = "grey",
      direction = -1,
      limits = c(0, 1)
    ) +
    labs(
      title = label,
      x = NULL,
      y = y_label,
      fill = "TPM (norm)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size = 10, face = "bold"),
      axis.title.y = element_text(size = 10)
    )
  
  # E) Create the mini summary heatmap (mean TPM per Tissue for this group).
  avg_df <- heatmap_df %>%
    group_by(Tissue) %>%
    summarise(mean_TPM = mean(avg_TPM, na.rm = TRUE), .groups = "drop")
  local_min <- min(avg_df$mean_TPM, na.rm = TRUE)
  local_max <- max(avg_df$mean_TPM, na.rm = TRUE)
  avg_df <- if (local_max - local_min == 0) {
    avg_df %>% mutate(norm_mean_TPM = 0)
  } else {
    avg_df %>% mutate(norm_mean_TPM = (mean_TPM - local_min) / (local_max - local_min))
  }
  avg_heatmap <- ggplot(avg_df, aes(x = Tissue, y = 1)) +
    geom_tile(aes(fill = norm_mean_TPM), color = "white") +
    scale_fill_viridis_c(
      option = "magma",
      na.value = "grey",
      direction = -1,
      limits = c(0, 1),
      guide = "none"
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank()
    )
  
  # Combine the main heatmap and the mini heatmap.
  combined_plot <- heatmap_plot / avg_heatmap +
    plot_layout(heights = c(1, 0.3))
  
  # Store the plot in the list.
  plots_list[[label]] <- combined_plot
}

# Reorder plots and combine into one vertical plot.
plots_list <- plots_list[order(
  grepl("Canonical", names(plots_list)),
  grepl("Protein", names(plots_list)),
  grepl("Non-coding", names(plots_list)),
  grepl("Pseudogenes", names(plots_list)),
  decreasing = TRUE
)]
if (length(plots_list) > 0) {
  final_plot <- wrap_plots(plots_list, ncol = 1) +
    plot_annotation(title = "Unique Heatmap (Canonical, Protein, Non-coding, Pseudogenes)")
  
  # Display on screen and save to file.
  print(final_plot)
  ggsave("unique_heatmap.png", final_plot, width = 12, height = 60, limitsize = FALSE)
}

# =============================================================================
#         SECTION 8: WOLF PSORT & INTERPROSCAN ANALYSIS
# =============================================================================
      
# Read the Wolf PSORT file which includes the Gene information.
psort_data <- read.delim("wolf_psort_table_trim27.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Split the combined ID column into Gene, TranscriptID, and ProteinID.
split_ids <- strsplit(psort_data$Gene.TranscriptID.ProteinID, "\\|")
psort_data$Gene <- sapply(split_ids, `[`, 1)
psort_data$TranscriptID <- sapply(split_ids, `[`, 2)
psort_data$ProteinID <- sapply(split_ids, `[`, 3)
psort_data$Gene.TranscriptID.ProteinID <- NULL

# Reorder columns (Gene will be removed later).
psort_data <- psort_data[, c("Gene", "TranscriptID", "ProteinID", "Localization", setdiff(names(psort_data), c("Gene", "TranscriptID", "ProteinID", "Localization")))]

# Read the InterProScan file.
interpro <- read.delim("iprscan5-R20250324-034509-0946-34713522-p1m.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
colnames(interpro)[1:13] <- c(
  "ProteinID", "MD5", "Length", "Analysis", "SignatureAccession", "SignatureDesc",
  "Start", "End", "Score", "Status", "Date", "InterPro_ID", "InterPro_Name"
)

# Define strict nuclear-related keywords.
keywords_nuclear <- c(
  "dna-binding", "transcription factor", "helix-turn-helix", "homeobox", "nuclear localization",
  "nucleoprotein", "histone", "zinc finger protein", "zinc finger domain", "classic zinc finger",
  "c2h2-type zinc finger"
)

# Clean the SignatureDesc column.
interpro$SignatureDesc <- ifelse(is.na(interpro$SignatureDesc), "", interpro$SignatureDesc)
interpro$SignatureDesc_clean <- tolower(trimws(interpro$SignatureDesc))

# Keep only rows with nuclear-related domains.
pattern <- paste(keywords_nuclear, collapse = "|")
interpro_filtered <- interpro[grepl(pattern, interpro$SignatureDesc_clean), ]

# Summarize the unique domain descriptions per ProteinID.
interpro_summary <- aggregate(
  SignatureDesc ~ ProteinID,
  data = interpro_filtered,
  FUN = function(x) paste(unique(x), collapse = "; ")
)
colnames(interpro_summary)[2] <- "InterPro_support"

# Split the ProteinID into TranscriptID and ProteinID.
id_split <- strsplit(interpro_summary$ProteinID, "\\|")
interpro_summary$TranscriptID <- sapply(id_split, `[`, 1)
interpro_summary$ProteinID <- sapply(id_split, `[`, 2)
interpro_summary <- interpro_summary[, c("TranscriptID", "ProteinID", "InterPro_support")]

# Merge the InterPro summary with the Wolf PSORT data.
psort_data <- merge(psort_data, interpro_summary, by = c("TranscriptID", "ProteinID"), all.x = TRUE)

# Generate a flextable for each gene and export to HTML and PDF.
unique_genes <- unique(psort_data$Gene)
for (gene in unique_genes) {
  dir.create(gene, showWarnings = FALSE)
  df_gene <- psort_data %>%
    filter(Gene == gene) %>%
    select(-Gene)  # Remove the Gene column
  
  ft <- flextable(df_gene)
  # Center numeric columns.
  numeric_cols <- names(df_gene)[sapply(df_gene, is.numeric)]
  ft <- align(ft, j = numeric_cols, align = "center", part = "all")
  ft <- align(ft, part = "header", align = "center")
  
  # Add a title in English.
  titre <- paste("Predicted subcellular localization for", gene)
  ft <- add_header_lines(ft, values = titre)
  ft <- align(ft, i = 1, align = "center", part = "header")
  
  # Save the flextable as HTML and then convert to PDF.
  html_file <- file.path(gene, paste0("wolf_psort_table_", gene, ".html"))
  pdf_file <- file.path(gene, paste0("wolf_psort_table_", gene, ".pdf"))
  save_as_html(ft, path = html_file)
  
  paper_height <- max(11, nrow(df_gene) * 0.22 + 1.5)  # Add a margin for the title
  pagedown::chrome_print(
    input = html_file,
    output = pdf_file,
    options = list(
      paperWidth = 11,
      paperHeight = paper_height,
      marginTop = 0,
      marginBottom = 0,
      marginLeft = 0,
      marginRight = 0
    )
  )
}
print("done")

# =============================================================================
#                SECTION 9: METADATA TABLE FOR TRIM27 & MID2
# =============================================================================

# Complete metadata for TRIM27 and MID2 transcripts.
tsl_data <- tibble(
  Transcript = c(
    # TRIM27
    "ENST00000377199.4", "ENST00000377194.7", "ENST00000414543.5",
    "ENST00000498117.1", "ENST00000496091.1", "ENST00000467742.1", "ENST00000481474.5",
    # MID2
    "ENST00000262843.11", "ENST00000443968.2", "ENST00000451923.1", "ENST00000474517.1"
  ),
  Name = c(
    # TRIM27
    "TRIM27-202", "TRIM27-201", "TRIM27-203",
    "TRIM27-207", "TRIM27-206", "TRIM27-204", "TRIM27-205",
    # MID2
    "MID2-201", "MID2-202", "MID2-203", "MID2-204"
  ),
  bp = c(2963, 2686, 1813, 728, 687, 541, 7006,
         7333, 2251, 758, 689),
  Protein = c("513aa", "358aa", "248aa", "No protein", "No protein", "No protein", "No protein",
              "735aa", "705aa", "218aa", "No protein"),
  Biotype = c("Protein coding", "Protein coding", "Protein coding",
              "Protein coding CDS not defined", "Protein coding CDS not defined",
              "Protein coding CDS not defined", "Retained intron",
              "Protein coding", "Protein coding", "Protein coding", "Retained intron"),
  UniProt_ID = c("P14373-1", "P14373-2", "H0Y551", NA, NA, NA, NA,
                 "Q9UJV3-1", "Q9UJV3-2", "A6PVI4", NA),
  TSL = c("TSL1", "TSL1", "TSL3", "TSL3", "TSL3", "TSL2", "TSL1",
          "TSL1", "TSL1", "TSL3", "TSL3")
)

# Filter for the two genes.
combined_table <- results_list_enriched %>%
  filter(Gene %in% c("TRIM27", "MID2")) %>%
  group_by(Transcript) %>%
  mutate(max_tpm = max(TPM, na.rm = TRUE)) %>%
  filter(TPM >= 0.85 * max_tpm) %>%
  arrange(Gene, Transcript, desc(is_canonical), desc(TPM)) %>%
  mutate(Tissue_TPM = paste0(Tissue, ": ", round(TPM, 2))) %>%
  summarise(
    Top_Expressed_Tissues = paste(Tissue_TPM, collapse = ", "),
    is_canonical = first(is_canonical),
    Gene = first(Gene),
    .groups = "drop"
  )

# Merge with metadata.
final_table <- combined_table %>%
  left_join(tsl_data, by = "Transcript") %>%
  arrange(Gene, desc(is_canonical)) %>%
  select(Gene, Transcript, Name, bp, Protein, Biotype, UniProt_ID,
         Top_Expressed_Tissues, TSL)

# Display the final table.
print(final_table)

# =============================================================================
#         SECTION 10: SUMMARY TABLE FOR ALL GENES WITH TSL & TISSUE INFO
# =============================================================================

# Step 1: Create the summary table for all genes.
all_genes_table <- results_list_enriched %>%
  group_by(Gene, Transcript) %>%
  mutate(max_tpm = max(TPM, na.rm = TRUE)) %>%
  filter(TPM >= 0.85 * max_tpm) %>%
  arrange(Gene, Transcript, desc(is_canonical), desc(TPM)) %>%
  mutate(Tissue_TPM = paste0(Tissue, ": ", round(TPM, 2))) %>%
  summarise(
    Top_Expressed_Tissues = paste(Tissue_TPM, collapse = ", "),
    is_canonical = first(is_canonical),
    .groups = "drop"
  ) %>%
  arrange(Gene, desc(is_canonical))

# Step 2: Read tsl_lvl.txt (3 columns) and rename columns.
tsl_info <- read_tsv("tsl_lvl.txt", show_col_types = FALSE)
colnames(tsl_info) <- c("Transcript", "Raw_TSL", "GENCODE_raw")

# Extract clean TSL and GENCODE Basic information.
tsl_info <- tsl_info %>%
  mutate(
    TSL = toupper(str_extract(Raw_TSL, "tsl[1-5]")),
    GENCODE_Basic = ifelse(str_detect(GENCODE_raw, "GENCODE basic"), TRUE, FALSE)
  ) %>%
  select(Transcript, TSL, GENCODE_Basic)

# Step 3: Merge with the main table.
all_genes_table <- all_genes_table %>%
  left_join(tsl_info, by = "Transcript") %>%
  mutate(GENCODE_Basic = ifelse(is.na(GENCODE_Basic), FALSE, GENCODE_Basic))

# Step 4: Add annotation status.
all_genes_table <- all_genes_table %>%
  mutate(
    Annotation_Status = case_when(
      TSL %in% c("TSL1", "TSL2") ~ "High confidence",
      TSL %in% c("TSL3", "TSL4", "TSL5") ~ "Low support",
      is.na(TSL) & GENCODE_Basic ~ "No TSL, but GENCODE Basic",
      is.na(TSL) & !GENCODE_Basic ~ "Likely poorly annotated",
      TRUE ~ "Unknown"
    )
  )

# Step 5: Create a flextable for the main table.
ft_main <- flextable(all_genes_table)
ft_main <- autofit(ft_main)

# Step 6: Create tissue summary table.
tissue_summary <- all_genes_table %>%
  select(Gene, Top_Expressed_Tissues) %>%
  separate_rows(Top_Expressed_Tissues, sep = ",\\s*") %>%
  mutate(Tissue = str_extract(Top_Expressed_Tissues, "^[^:]+")) %>%
  count(Gene, Tissue, sort = TRUE) %>%
  arrange(Gene, desc(n))

# Step 7: Create a flextable for the tissue summary.
ft_summary <- flextable(tissue_summary)
ft_summary <- autofit(ft_summary)
ft_summary <- add_header_lines(ft_summary, values = "Tissue summary by gene")

# Step 8: Combine both tables into a single HTML document.
save_as_html(ft_main, path = "main_table.html")
save_as_html(ft_summary, path = "summary_table.html")
main_content <- readLines("main_table.html")
summary_content <- readLines("summary_table.html")

# Extract the <body> content from both HTML files.
main_body <- main_content[which(grepl("<body>", main_content)):which(grepl("</body>", main_content))]
summary_body <- summary_content[which(grepl("<body>", summary_content)):which(grepl("</body>", summary_content))]

html_full <- c(
  '<html><head><style>@page { size: landscape; }</style></head>',
  '<body>',
  main_body[-c(1, length(main_body))],  # Remove <body> and </body> tags.
  '<hr><h2>Summary of tissue counts per gene</h2>',
  summary_body[-c(1, length(summary_body))],
  '</body></html>'
)

# Save the merged HTML and generate the PDF.
writeLines(html_full, "sum_isoform.html")
pagedown::chrome_print("sum_isoform.html", output = "sum_isoform.pdf")
