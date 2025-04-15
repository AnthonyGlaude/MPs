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


#les inputs sont : clusterdata_combined.txt (Gene, Tissue, Cluster, Transcript, TPM)
#                : genes_list.tsv 

### PARTIE 1 : Lecture et préparation des données d'expression ###



# 1. Lecture des deux jeux de données (si ce n'est pas déjà fait)
results_list <- read_tsv("clusterdata_combined.txt", show_col_types = FALSE)
mane_select <- read_tsv("reference/mane_select_trim27.tsv", show_col_types = FALSE)

# 2. Jointure gauche de results_list avec mane_select sur la colonne des transcrits
results_list_enriched <- results_list %>%
  left_join(
    mane_select %>% select(`Transcript stable ID version`, 
                           `Protein stable ID version`, 
                           `RefSeq match transcript (MANE Select)`, 
                           `UniProtKB isoform ID`,
                           `UniProtKB/TrEMBL ID`),  
    by = c("Transcript" = "Transcript stable ID version")
  ) %>%
  
  # 3. Renommage des nouvelles colonnes pour plus de clarté
  rename(
    ProteinID = `Protein stable ID version`,
    MANE_Transcript = `RefSeq match transcript (MANE Select)`,
    UniProt_ID_isoform = `UniProtKB isoform ID`,
    UniProt_ID_trembl = `UniProtKB/TrEMBL ID`  
  ) %>%
  mutate(
    UniProt_ID = coalesce(UniProt_ID_isoform, UniProt_ID_trembl)  
  ) %>%
  select(-UniProt_ID_isoform, -UniProt_ID_trembl) 


results_list_enriched <- results_list_enriched %>%
  mutate(is_canonical = !is.na(MANE_Transcript) & MANE_Transcript != "")



### PARTIE 4 : Calcul des ratios d'expression ###

results_list_enriched <- results_list_enriched %>%
  group_by(Gene, Tissue, Cluster) %>%
  mutate(ref_tpm = TPM[is_canonical == TRUE][1],
         ratio = TPM / ref_tpm) %>%
  ungroup()

agg_df <- results_list_enriched %>%
  group_by(Gene, Transcript, Cluster, is_canonical) %>%
  summarise(
    count = n(),
    ratios = if (!first(is_canonical)) paste(round(ratio, 3), collapse = "; ") else NA_character_,
    mean_ratio = if (!first(is_canonical)) round(mean(ratio), 3) else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    cell = ifelse(is_canonical,
                  paste0("n = ", count),
                  ifelse(!is.na(ratios), paste0(ratios, " (", mean_ratio, ")"), NA_character_)
    )
  )


### PARTIE 5 : Mise en forme du tableau récapitulatif ###

# Passage en format large (une colonne par cluster)
wide_df <- agg_df %>%
  mutate(cluster_col = paste0("cluster_", Cluster)) %>%
  select(Gene, Transcript, cluster_col, cell, is_canonical) %>%
  pivot_wider(names_from = cluster_col, values_from = cell)

wide_df <- wide_df %>%
  arrange(Gene, desc(is_canonical))

# Ajouter la mention (REF) pour le transcript de référence
wide_df <- wide_df %>%
  mutate(Transcript = ifelse(Transcript %in% unique(results_list_enriched$Transcript[results_list_enriched$is_canonical]),
                             paste0(Transcript, " (REF)"),
                             Transcript))

# Création d'un flextable et export en HTML puis PDF
ft <- flextable(wide_df)
ft <- autofit(ft)
ft
save_as_html(ft, path = "expression_isoform.html")
pagedown::chrome_print(input = "expression_isoform.html", output = "expression_isoform.pdf")

### PARTIE 6 : Préparation du tableau final (moyenne des clusters) ###

mean_df <- agg_df %>%
  select(Gene, Transcript, Cluster, mean_ratio) %>%
  pivot_wider(names_from = Cluster, values_from = mean_ratio, names_prefix = "cluster_")

mean_df <- mean_df %>%
  arrange(Gene, desc(Transcript %in% unique(results_list_enriched$Transcript[results_list_enriched$is_canonical])))

mean_df <- mean_df %>%
  mutate(Transcript = ifelse(Transcript %in% unique(results_list_enriched$Transcript[results_list_enriched$is_canonical]),
                             paste0(Transcript, " (REF)"),
                             Transcript))

mean_df_long <- mean_df %>%
  pivot_longer(
    cols = starts_with("cluster_"),
    names_to = "Cluster",
    values_to = "mean_ratio"
  )

# Remplacement par NA des valeurs inférieures au quantile 75 pour chaque groupe (Gene, Cluster)
mean_df_long <- mean_df_long %>%
  group_by(Gene, Cluster) %>%
  mutate(mean_ratio = ifelse(mean_ratio > quantile(mean_ratio, probs = 0.75, na.rm = TRUE),
                             mean_ratio, NA)) %>%
  ungroup()

# Reconversion en format large (optionnel)
mean_df_filtered_wide <- mean_df_long %>%
  pivot_wider(
    names_from = Cluster,
    values_from = mean_ratio
  )

mean_df_filtered_wide <- mean_df_filtered_wide %>%
  filter(str_detect(Transcript, fixed(" (REF)")) | !if_all(starts_with("cluster_"), is.na))


















#############ADD_localisation############################
#############Utilisation similarité cosinus #############
library(dplyr)
library(tidyr)
library(stringr)

#### Partie A : Préparation de mean_df_filtered_wide ####
mean_df_filtered_wide <- mean_df_filtered_wide %>%
  mutate(TranscriptID = str_remove(Transcript, fixed(" (REF)")),
         is_ref = str_detect(Transcript, fixed(" (REF)")))

#### Partie B : Construction de la matrice de localisations ####
# Lire et parser le fichier de prédictions

pred_lines <- readLines("isoforms_localisation.txt")
pred_raw <- tibble(raw_line = pred_lines) %>%
  filter(raw_line != "") %>%
  mutate(
    Transcript = sub("^(\\S+).*", "\\1", raw_line),
    details = ifelse(str_detect(raw_line, "details"),
                     sub(".*details\\s*", "", raw_line),
                     NA)
  )

# Séparer les paires "localisation: score"
pred_details <- pred_raw %>%
  filter(!is.na(details)) %>%
  separate_rows(details, sep = ",\\s*") %>%
  separate(details, into = c("Localization", "Score"), sep = ":\\s*", convert = TRUE)

# Agréger les scores par Transcript et Localization
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



# Construire la matrice de localisation : un vecteur par Transcript
localization_matrix <- agg_predictions %>%
  select(Transcript, Localization, Percentage) %>%
  pivot_wider(names_from = Localization, values_from = Percentage, 
              values_fill = list(Percentage = 0)) %>%
  mutate(TranscriptID = Transcript)

#### Partie C : Calcul de la similarité cosinus ####
cosine_similarity <- function(vec1, vec2) {
  dot_product <- sum(vec1 * vec2)
  norm1 <- sqrt(sum(vec1^2))
  norm2 <- sqrt(sum(vec2^2))
  return(dot_product / (norm1 * norm2))
}

# Déterminer les colonnes de localisation (score) à comparer
localization_cols <- setdiff(names(localization_matrix), c("Transcript", "TranscriptID"))

# Remplacer d'éventuels NA par 0 (en principe, pivot_wider les a remplacé)
localization_matrix <- localization_matrix %>%
  mutate(across(all_of(localization_cols), ~ replace_na(.x, 0)))

# Pour chaque gène, joindre l'information de Gene et de référence à partir de mean_df_filtered_wide
localization_matrix <- localization_matrix %>%
  left_join(mean_df_filtered_wide %>% select(Gene, TranscriptID, is_ref),
            by = "TranscriptID")

# Pour chaque gène, calculer la similarité cosinus en comparant chaque transcript à celui de référence
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

#### Partie D : Intégration dans mean_df_filtered_wide ####
final_df <- left_join(mean_df_filtered_wide, 
                      similarity_by_gene %>% select(TranscriptID, CosineSimilarity),
                      by = "TranscriptID")

View(final_df)

#########################################################################


# Lecture du fichier texte en utilisant un autre nom
protnlm_lines <- readLines("protnlm.txt")

# Identifier les lignes de début de bloc transcript (lignes commençant par ">")
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
  
  # Si c'est le premier transcript pour le gène, on le considère comme référence
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
  
  # Comparaison de la fonction courante avec la référence
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
    GeneName = gene_name,
    Transcript_id = transcript_id,
    Fonction1 = current_func,
    Score1 = current_score,
    MatchReference = is_match,
    stringsAsFactors = FALSE
  )
}

final_results_df <- do.call(rbind, result_list)
print(final_results_df)




# Pour extraire uniquement les résultats pertinents pour chaque gène
final_results <- list()
unique_genes <- unique(final_results_df$GeneName)

for (gene in unique_genes) {
  gene_data <- final_results_df[final_results_df$GeneName == gene, ]
  gene_data <- gene_data[!is.na(gene_data$Fonction1), ]
  final_results[[gene]] <- gene_data
}

final_results_df <- do.call(rbind, final_results)
print(final_results_df)
final_results_df$TranscriptID <- sub("^>", "", final_results_df$Transcript_id)
# Renommer la colonne "Transcript_id" en "TranscriptID"
final_results_df$TranscriptID <- sapply(strsplit(final_results_df$Transcript_id, "\\|"), `[`, 1)

# Puis, lors du merge, ne sélectionner que les colonnes disponibles
merged_df <- merge(final_df, 
                   final_results_df[, c("TranscriptID", "Fonction1", "Score1")], 
                   by = "TranscriptID", 
                   all.x = TRUE)


# dataframe "final" (faudrait changer de nom)
merged_df <- subset(merged_df, select = -c(is_ref, TranscriptID)) 


####### Tableau #######
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


































# Pour chaque gène, on crée un champ "canonical_status" basé sur is_canonical
results_list_enriched <- results_list_enriched %>%
  mutate(canonical_status = if_else(is_canonical, "Canonical", "Non-canonical"))

# Extraire la liste des gènes à traiter
genes <- unique(results_list_enriched$Gene)

for (g in genes) {
  
  # Filtrer les données pour le gène courant
  gene_data <- results_list_enriched %>% filter(Gene == g)
  
  # Nettoyer UniProt_ID
  gene_data <- gene_data %>%
    mutate(UniProt_ID = na_if(trimws(UniProt_ID), "")) %>%
    mutate(UniProt_ID = na_if(UniProt_ID, "NA")) %>%
    mutate(UniProt_ID = na_if(UniProt_ID, "NaN"))
  
  # Identifier les UniProt_ID des transcrits canoniques
  canonical_uniprots <- gene_data %>%
    filter(is_canonical, !is.na(UniProt_ID)) %>%
    pull(UniProt_ID) %>%
    unique()
  
  gene_data <- gene_data %>%
    mutate(group_label = case_when(
      !is.na(UniProt_ID) & UniProt_ID %in% canonical_uniprots ~ paste0("Canonical: ", UniProt_ID),
      !is.na(UniProt_ID) ~ paste0("Protein: ", UniProt_ID),
      TRUE ~ "Non-coding transcripts"
    ))
  
  # Agrégation uniquement par group_label
  protein_groups <- gene_data %>%
    group_by(group_label) %>%
    summarise(transcripts = list(unique(Transcript)), .groups = "drop")
  
  plots_list <- list()
  
  for (i in seq_len(nrow(protein_groups))) {
    
    label <- protein_groups$group_label[i]
    transcripts <- protein_groups$transcripts[[i]]
    
    # Extraire les données pour ce sous-groupe
    grp_data <- gene_data %>% filter(Transcript %in% transcripts)
    
    # Moyenne TPM par transcript et tissu
    heatmap_df <- grp_data %>%
      group_by(Transcript, Tissue) %>%
      summarise(avg_TPM = mean(TPM, na.rm = TRUE), .groups = "drop")
    
    # Normalisation
    global_min <- min(heatmap_df$avg_TPM, na.rm = TRUE)
    global_max <- max(heatmap_df$avg_TPM, na.rm = TRUE)
    
    heatmap_df <- if (global_max - global_min == 0) {
      heatmap_df %>% mutate(norm_TPM = 0)
    } else {
      heatmap_df %>% mutate(norm_TPM = (avg_TPM - global_min) / (global_max - global_min))
    }
    
    # Garder ordre logique : canonique d'abord
    heatmap_df <- heatmap_df %>%
      mutate(
        is_can_transcript = Transcript %in% gene_data$Transcript[gene_data$is_canonical],
        Transcript = factor(
          Transcript,
          levels = c(unique(Transcript[is_can_transcript]),
                     unique(Transcript[!is_can_transcript]))
        )
      )
    
    # Heatmap principal
    angle_x <- 60
    heatmap_plot <- ggplot(heatmap_df, aes(x = Tissue, y = Transcript)) +
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
        y = "Transcript",
        fill = "TPM (norm)"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 10, face = "bold")
      )
    
    # Mini heatmap résumé par tissu
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
      labs(
        title = NULL,
        x = "Tissue",
        y = NULL
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = angle_x, vjust = 1, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()
      )
    
    combined_plot <- heatmap_plot / avg_heatmap + 
      plot_layout(heights = c(1, 0.3))
    
    plots_list[[label]] <- combined_plot
  }
  
  # Trier les groupes (canonique > protein > non-coding)
  plots_list <- plots_list[order(
    grepl("Canonical", names(plots_list)),
    grepl("Protein", names(plots_list)),
    decreasing = TRUE
  )]
  
  # Sauvegarde PNG par gène
  if (length(plots_list) > 0) {
    combined_gene_plot <- wrap_plots(plots_list, ncol = 1) +
      plot_annotation(title = paste("Heatmap for gene:", g))
    
    dir.create(g, showWarnings = FALSE)
    ggsave(
      filename = paste0(g, "/heatmap_", g, ".png"),
      plot = combined_gene_plot,
      width = 10,
      height = 9 + length(plots_list) * 2,
      limitsize = FALSE  
    )
  }
}

print("bob")

# === PDF résumé à partir des PNG déjà générés ===

pdf("summary_heatmaps_from_pngs.pdf", width = 10, height = 11)

for (g in genes) {
  image_path <- file.path(g, paste0("heatmap_", g, ".png"))
  
  if (file.exists(image_path)) {
    img <- readPNG(image_path)
    grid::grid.raster(img)
    grid::grid.text(g, x = 0.5, y = 0.98, gp = gpar(fontsize = 14, fontface = "bold"))
    grid::grid.newpage()
  } else {
    warning(paste("Image manquante pour le gène :", g))
  }
}

dev.off()


print("bob")


















# === Lire les lignes ===
pred_lines <- readLines("isoforms_localisation.txt")

# === Table de correspondance Wolf → HPA ===
localization_mapping <- tribble(
  ~Wolf,         ~HPA,
  "nucl",        "Nucleus",
  "cyto",        "Cytosol",
  "cyto_nucl",   "Nucleus and Cytosol",  # modifié : localisation entière pour cyto_nucl
  "nucleop",     "Nucleoplasm",
  "extr",        "Extracellular",
  "plas",        "Plasma membrane",
  "mito",        "Mitochondria",
  "er",          "Endoplasmic reticulum",
  "golgi",       "Golgi apparatus",
  "lys",         "Lysosomes",
  "pero",        "Peroxisome"
)

# === Extraire ENSP + prédictions ===
pred_df <- tibble(raw = pred_lines) %>%
  filter(str_detect(raw, "\\S")) %>%
  mutate(
    ProteinID = str_extract(raw, "ENSP[0-9]+\\.[0-9]+"),
    details = str_extract(raw, "details\\s+.*")
  ) %>%
  filter(!is.na(ProteinID), !is.na(details)) %>%
  mutate(details = str_remove(details, "^details\\s+"))

# === Séparer localisation:score ===
loc_long <- pred_df %>%
  separate_rows(details, sep = ",\\s*") %>%
  separate(details, into = c("Localization", "Score"), sep = ":\\s*", convert = TRUE) %>%
  mutate(Score = as.numeric(Score)) %>%
  filter(!is.na(Score))

# Note : plus de traitement spécifique pour "cyto_nucl", elle sera gardée telle quelle

# === Appliquer le mapping Wolf → HPA ===
loc_standardized <- loc_long %>%
  left_join(localization_mapping, by = c("Localization" = "Wolf")) %>%
  mutate(Localization = coalesce(HPA, Localization)) %>%
  select(-HPA)

# === Calcul des % par localisation ===
localization_percent <- loc_standardized %>%
  group_by(ProteinID, Localization) %>%
  summarise(TotalScore = sum(Score, na.rm = TRUE), .groups = "drop") %>%
  group_by(ProteinID) %>%
  mutate(
    Total = sum(TotalScore, na.rm = TRUE),
    Percentage = ifelse(Total > 0, TotalScore / Total * 100, 0)
  ) %>%
  ungroup()

# === Ajouter infos Gene + UniProt ===
localization_percent <- localization_percent %>%
  left_join(
    results_list_enriched %>%
      select(Gene, ProteinID, UniProt_ID) %>%
      distinct(),
    by = "ProteinID"
  ) %>%
  mutate(ProteinLabel = ifelse(!is.na(UniProt_ID), UniProt_ID, ProteinID))

# === Lire fichier Human Protein Atlas ===
hpa_validation <- read.delim("hpa_localization.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# === Fusionner avec HPA ===
loc_long_ready <- localization_percent %>%
  left_join(hpa_validation, by = c("Gene", "Localization")) %>%
  select(Gene, ProteinLabel, Localization, Percentage, Status)

# === Générer les heatmaps ===
genes <- unique(loc_long_ready$Gene)

for (gene in genes) {
  cat("Traitement du gène:", gene, "\n")
  gene_loc_data <- loc_long_ready %>% filter(Gene == gene)
  if (nrow(gene_loc_data) == 0) next
  
  p_loc <- ggplot(gene_loc_data, aes(x = Localization, y = ProteinLabel)) +
    # Cercle externe : validation HPA (plus gros)
    geom_point(
      data = gene_loc_data %>% filter(Status %in% c("approved", "enhanced")),
      aes(color = Status),
      shape = 1, stroke = 2, size = 7
    ) +
    # Point principal : taille fixe
    geom_point(
      aes(fill = ifelse(Percentage == 0, NA, Percentage)),
      shape = 21, size = 6.5, color = "black"
    ) +
    scale_fill_viridis_c(option = "magma", na.value = "transparent", direction = -1) +
    scale_color_manual(
      values = c(
        "approved" = "green3",
        "enhanced" = "gold",
        "supported" = "dodgerblue",
        "uncertain" = "gray50"
      ),
      name = "Experimental validation (HPA)"
    ) +
    labs(
      title = paste("Gene subcellular localization prediction –", gene),
      x = "Subcellular localization (Wolf PSORT)",
      y = "Protein (UniProt)",
      fill = "Prediction (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA)
    )
  
  dir.create(gene, showWarnings = FALSE)
  ggsave(
    filename = file.path(gene, paste0("cellular_localization_", gene, ".png")),
    plot = p_loc,
    width = 10,
    height = 8
  )
}



















############### Tableau again ################
############### Script complet ################

library(flextable)
library(pagedown)
library(dplyr)

# === Étape 1 : lire le fichier Wolf PSORT avec Gene inclus
psort_data <- read.delim("wolf_psort_table_trim27.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Séparer Gene, TranscriptID et ProteinID
split_ids <- strsplit(psort_data$Gene.TranscriptID.ProteinID, "\\|")
psort_data$Gene <- sapply(split_ids, `[`, 1)
psort_data$TranscriptID <- sapply(split_ids, `[`, 2)
psort_data$ProteinID <- sapply(split_ids, `[`, 3)
psort_data$Gene.TranscriptID.ProteinID <- NULL

# Réorganiser les colonnes (Gene sera retiré plus tard)
psort_data <- psort_data[, c("Gene", "TranscriptID", "ProteinID", "Localization", setdiff(names(psort_data), c("Gene", "TranscriptID", "ProteinID", "Localization")))]

# === Étape 2 : lire le fichier InterProScan
interpro <- read.delim("iprscan5-R20250324-034509-0946-34713522-p1m.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
colnames(interpro)[1:13] <- c(
  "ProteinID", "MD5", "Length", "Analysis", "SignatureAccession", "SignatureDesc",
  "Start", "End", "Score", "Status", "Date", "InterPro_ID", "InterPro_Name"
)

# === Étape 3 : définir les mots-clés NUCLEAIRES (filtrage strict)
keywords_nuclear <- c(
  "dna-binding", "transcription factor", "helix-turn-helix", "homeobox", "nuclear localization",
  "nucleoprotein", "histone", "zinc finger protein", "zinc finger domain", "classic zinc finger",
  "c2h2-type zinc finger"
)

# Nettoyage de la colonne SignatureDesc
interpro$SignatureDesc <- ifelse(is.na(interpro$SignatureDesc), "", interpro$SignatureDesc)
interpro$SignatureDesc_clean <- tolower(trimws(interpro$SignatureDesc))

# === Étape 4 : ne garder que les domaines contenant les mots-clés nucléaires
pattern <- paste(keywords_nuclear, collapse = "|")
interpro_filtered <- interpro[grepl(pattern, interpro$SignatureDesc_clean), ]

# === Étape 5 : Résumer uniquement les descriptions pertinentes par ProteinID (complet)
interpro_summary <- aggregate(
  SignatureDesc ~ ProteinID,
  data = interpro_filtered,
  FUN = function(x) paste(unique(x), collapse = "; ")
)
colnames(interpro_summary)[2] <- "InterPro_support"

# === Étape 6 : séparer TranscriptID et ProteinID dans interpro_summary
id_split <- strsplit(interpro_summary$ProteinID, "\\|")
interpro_summary$TranscriptID <- sapply(id_split, `[`, 1)
interpro_summary$ProteinID <- sapply(id_split, `[`, 2)

# Réorganiser les colonnes
interpro_summary <- interpro_summary[, c("TranscriptID", "ProteinID", "InterPro_support")]

# === Joindre au tableau Wolf PSORT
psort_data <- merge(psort_data, interpro_summary, by = c("TranscriptID", "ProteinID"), all.x = TRUE)

unique_genes <- unique(psort_data$Gene)

for (gene in unique_genes) {
  dir.create(gene, showWarnings = FALSE)
  df_gene <- psort_data %>%
    filter(Gene == gene) %>%
    select(-Gene)  # Retirer la colonne Gene du tableau
  
  ft <- flextable(df_gene)
  
  # Centrer les colonnes numériques
  numeric_cols <- names(df_gene)[sapply(df_gene, is.numeric)]
  ft <- align(ft, j = numeric_cols, align = "center", part = "all")
  ft <- align(ft, part = "header", align = "center")
  
  # Ajouter un titre en français
  titre <- paste("Predicted subcellular localization for", gene)
  ft <- add_header_lines(ft, values = titre)
  ft <- align(ft, i = 1, align = "center", part = "header")
  
  # Sauvegarde
  html_file <- file.path(gene, paste0("wolf_psort_table_", gene, ".html"))
  pdf_file <- file.path(gene, paste0("wolf_psort_table_", gene, ".pdf"))
  
  save_as_html(ft, path = html_file)
  
  paper_height <- max(11, nrow(df_gene) * 0.22 + 1.5)  # Légère marge pour le titre
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













library(dplyr)
library(tibble)
library(stringr)

# Metadonnées complètes pour TRIM27 + MID2
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

# Filtrer pour les deux gènes
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

# Fusion avec métadonnées
final_table <- combined_table %>%
  left_join(tsl_data, by = "Transcript") %>%
  arrange(Gene, desc(is_canonical)) %>%
  select(Gene, Transcript, Name, bp, Protein, Biotype, UniProt_ID,
         Top_Expressed_Tissues, TSL)

# Affichage
print(final_table)







library(dplyr)
library(flextable)
library(pagedown)
library(readr)
library(stringr)
library(tidyr)

# Step 1: Create the summary table for all genes
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

# Step 2: Read tsl_lvl.txt with 3 columns
tsl_info <- read_tsv("tsl_lvl.txt", show_col_types = FALSE)
colnames(tsl_info) <- c("Transcript", "Raw_TSL", "GENCODE_raw")

# Extract clean TSL and GENCODE Basic info
tsl_info <- tsl_info %>%
  mutate(
    TSL = str_extract(Raw_TSL, "tsl[1-5]"),
    TSL = toupper(TSL),
    GENCODE_Basic = ifelse(str_detect(GENCODE_raw, "GENCODE basic"), TRUE, FALSE)
  ) %>%
  select(Transcript, TSL, GENCODE_Basic)

# Step 3: Merge with main table
all_genes_table <- all_genes_table %>%
  left_join(tsl_info, by = "Transcript") %>%
  mutate(GENCODE_Basic = ifelse(is.na(GENCODE_Basic), FALSE, GENCODE_Basic))

# Step 4: Add annotation status
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

# Step 5: Create flextable for main table
ft_main <- flextable(all_genes_table)
ft_main <- autofit(ft_main)

# Step 6: Create tissue summary table
tissue_summary <- all_genes_table %>%
  select(Gene, Top_Expressed_Tissues) %>%
  separate_rows(Top_Expressed_Tissues, sep = ",\\s*") %>%
  mutate(Tissue = str_extract(Top_Expressed_Tissues, "^[^:]+")) %>%
  count(Gene, Tissue, sort = TRUE) %>%
  arrange(Gene, desc(n))

# Step 7: Create flextable for summary
ft_summary <- flextable(tissue_summary)
ft_summary <- autofit(ft_summary)
ft_summary <- add_header_lines(ft_summary, values = "Tissue summary by gene")

# Step 8: Combine both tables into one HTML
# Create HTML for main table
save_as_html(ft_main, path = "main_table.html")

# Create HTML for summary table
save_as_html(ft_summary, path = "summary_table.html")

# Lire les deux HTML
main_content <- readLines("main_table.html")
summary_content <- readLines("summary_table.html")

# Extraire <body> des deux
main_body <- main_content[which(grepl("<body>", main_content)):which(grepl("</body>", main_content))]
summary_body <- summary_content[which(grepl("<body>", summary_content)):which(grepl("</body>", summary_content))]

# Fusionner dans un seul document
html_full <- c(
  '<html><head><style>@page { size: landscape; }</style></head>',
  '<body>',
  main_body[-c(1, length(main_body))],  # Retirer <body> et </body>
  '<hr><h2>Summary of tissue counts per gene</h2>',
  summary_body[-c(1, length(summary_body))],
  '</body></html>'
)

# Sauvegarder le fichier final HTML
writeLines(html_full, "sum_isoform.html")

# Générer le PDF
pagedown::chrome_print("sum_isoform.html", output = "sum_isoform.pdf")
