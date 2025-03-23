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

#les inputs sont : clusterdata_combined.txt (Gene, Tissue, Cluster, Transcript, TPM)
#                : genes_list.tsv 

### PARTIE 1 : Lecture et préparation des données d'expression ###
results_list <- read.delim("clusterdata_combined.txt", 
                           header = TRUE, sep = "\t", 
                           stringsAsFactors = FALSE)
colnames(results_list) <- c("Gene", "Tissue", "Cluster", "Transcript", "TPM")
results_list$TPM <- as.numeric(results_list$TPM)

### PARTIE 2 : Lecture et intégration du fichier MANE select ###
mane_select <- read.delim("reference/mane_select_trim27.tsv", 
                          header = TRUE, sep = "\t", 
                          stringsAsFactors = FALSE)
# Vérifiez les noms de colonnes
print(colnames(mane_select))

# Extraire le transcript canonique : on ne garde que les lignes où la colonne de MANE est renseignée
canonical_transcripts <- mane_select %>%
  filter(RefSeq.match.transcript..MANE.Select. != "") %>%
  select(Gene.name, Transcript.stable.ID.version) %>%
  rename(Gene = Gene.name, ref_transcript = Transcript.stable.ID.version)

# Fusionner avec results_list et créer la colonne booléenne is_canonical
results_list <- results_list %>%
  left_join(canonical_transcripts, by = "Gene") %>%
  mutate(is_canonical = (Transcript == ref_transcript))

# Vérification rapide
print(table(results_list$is_canonical))


### PARTIE 4 : Calcul des ratios d'expression ###

results_list <- results_list %>%
  group_by(Gene, Tissue, Cluster) %>%
  mutate(ref_tpm = TPM[is_canonical == TRUE][1],
         ratio = TPM / ref_tpm) %>%
  ungroup()

agg_df <- results_list %>%
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
  mutate(Transcript = ifelse(Transcript %in% unique(results_list$Transcript[results_list$is_canonical]),
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
  arrange(Gene, desc(Transcript %in% unique(results_list$Transcript[results_list$is_canonical])))

mean_df <- mean_df %>%
  mutate(Transcript = ifelse(Transcript %in% unique(results_list$Transcript[results_list$is_canonical]),
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











































library(dplyr)
library(stringr)

results_list <- results_list %>%
  # Nettoyer les espaces superflus dans RefSeq_Match et UniProt_Match
  mutate(
    RefSeq_Match = str_trim(RefSeq_Match),
    UniProt_Match = str_trim(UniProt_Match),
    UniProt_Match = if_else(UniProt_Match == "-", "no protein", UniProt_Match)
  ) %>%
  group_by(Gene) %>%
  # Réordonner pour que les transcrits avec une valeur non vide dans RefSeq_Match soient en tête
  arrange(desc(RefSeq_Match != "")) %>%
  mutate(canonical_uniprot = if_else(
    any(RefSeq_Match != ""),
    first(UniProt_Match),
    NA_character_
  )) %>%
  ungroup() %>%
  mutate(is_canonical = if_else(
    (RefSeq_Match != "") | (UniProt_Match == canonical_uniprot),
    TRUE, FALSE
  ))



library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(stringr)

# Avant la boucle, vous avez déjà votre data frame results_list
# ... (bloc de nettoyage et définition de canonical_uniprot, is_canonical) ...

# --- Boucle principale sur les gènes ---
for(g in genes) {
  
  # Filtrer et nettoyer les données pour le gène g
  gene_data <- results_list %>% 
    filter(Gene == g) %>%
    mutate(UniProt_Match = str_trim(UniProt_Match)) %>%  
    mutate(UniProt_Match = if_else(UniProt_Match == "-", "no protein", UniProt_Match))
  
  # Récupérer les groupes uniques de UniProt_Match
  uniprot_groups <- unique(gene_data$UniProt_Match)
  uniprot_groups <- uniprot_groups[!is.na(uniprot_groups) & uniprot_groups != ""]
  
  # Déterminer min_val et max_val globaux pour le gène (si vous le souhaitez,
  # mais ici on fera surtout un min/max local à chaque groupe)
  min_val <- min(gene_data$TPM, na.rm = TRUE)
  max_val <- max(gene_data$TPM, na.rm = TRUE)
  
  # Ordonner les groupes pour que le canonical apparaisse en premier et "no protein" en dernier
  can_up <- unique(gene_data$canonical_uniprot)
  if(length(can_up) > 1) { can_up <- can_up[1] }
  
  rank_uniprot <- function(up) {
    if(up == can_up) return(0)
    if(up == "no protein") return(2)
    return(1)
  }
  uniprot_groups <- uniprot_groups[order(sapply(uniprot_groups, rank_uniprot))]
  
  plots_list <- list()
  
  for(up in uniprot_groups) {
    
    grp_data <- gene_data %>% filter(UniProt_Match == up)
    
    # 1) Calcul de avg_TPM globalement (tous transcripts/tissues du groupe)
    heatmap_df <- grp_data %>%
      group_by(Transcript, Tissue) %>%
      summarise(avg_TPM = mean(TPM, na.rm = TRUE), .groups = "drop")
    
    # 1) Normalisation GLOBALE pour le heatmap principal
    global_min <- min(heatmap_df$avg_TPM, na.rm = TRUE)
    global_max <- max(heatmap_df$avg_TPM, na.rm = TRUE)
    
    if(global_max - global_min == 0) {
      heatmap_df <- heatmap_df %>% mutate(norm_TPM = 0)
    } else {
      heatmap_df <- heatmap_df %>% 
        mutate(norm_TPM = (avg_TPM - global_min) / (global_max - global_min))
    }
    
    # Réordonner les transcripts pour que les canoniques apparaissent en premier
    heatmap_df <- heatmap_df %>%
      mutate(
        is_can_transcript = Transcript %in% grp_data$Transcript[grp_data$is_canonical],
        Transcript = factor(
          Transcript,
          levels = c(unique(Transcript[is_can_transcript]),
                     unique(Transcript[!is_can_transcript]))
        )
      )
    
    # Heatmap principal (norm_TPM)
    angle_x <- 60  
    heatmap_plot <- ggplot(heatmap_df, aes(x = Tissue, y = Transcript)) +
      geom_tile(aes(fill = norm_TPM), color = "white") +
      scale_fill_viridis_c(
        option = "magma",
        na.value = "grey",
        direction = -1,
        limits = c(0, 1)  # car on a normalisé 0–1
      ) +
      labs(
        title = paste("UniProt:", up),
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
    
    # 2) Mini heatmap : moyenne TPM par Tissue, puis normalisation sur la ligne
    avg_df <- heatmap_df %>%
      group_by(Tissue) %>%
      summarise(mean_TPM = mean(avg_TPM, na.rm = TRUE), .groups = "drop")
    
    local_min <- min(avg_df$mean_TPM, na.rm = TRUE)
    local_max <- max(avg_df$mean_TPM, na.rm = TRUE)
    
    if(local_max - local_min == 0) {
      avg_df <- avg_df %>% mutate(norm_mean_TPM = 0)
    } else {
      avg_df <- avg_df %>% mutate(norm_mean_TPM = (mean_TPM - local_min) / (local_max - local_min))
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
    
    # Combinaison
    combined_grp_plot <- heatmap_plot / avg_heatmap + 
      plot_layout(heights = c(1, 0.3))
    
    plots_list[[up]] <- combined_grp_plot
  }
  
  if(length(plots_list) > 0) {
    combined_gene_plot <- wrap_plots(plots_list, ncol = 1) +
      plot_annotation(title = paste("Heatmap for gene:", g))
    
    dir.create(g, showWarnings = FALSE)
    ggsave(
      filename = paste0(g, "/heatmap_", g, ".png"),
      plot = combined_gene_plot,
      width = 10,
      height = 9
    )
  }
}












#####################FIGURE LOCALISATION CELL (WOLF)#######################



# Transformation des données de localisation en format long
localization_cols <- setdiff(names(localization_matrix), c("Transcript", "TranscriptID", "Gene", "is_ref"))
loc_long <- localization_matrix %>%
  pivot_longer(cols = all_of(localization_cols), 
               names_to = "Localization", 
               values_to = "Percentage")

# Pour chaque gène, créer un nuage de points avec fond blanc
genes <- unique(loc_long$Gene)
for (gene in genes) {
  gene_loc_data <- loc_long %>% filter(Gene == gene)
  if(nrow(gene_loc_data) == 0) next
  
  p_loc <- ggplot(gene_loc_data, aes(x = Transcript, y = Localization)) +
    geom_point(aes(color = ifelse(Percentage == 0, NA, Percentage),
                   size = Percentage)) +
    scale_color_viridis_c(option = "magma", na.value = "transparent", direction = -1) +
    scale_size_continuous(range = c(3, 10)) +
    labs(title = paste("Localisation cellulaire pour le gène :", gene),
         x = "Transcripts", y = "Localisations", 
         color = "Pourcentage", size = "Pourcentage") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA)
    )
  
  
  # Créer le dossier du gène s'il n'existe pas déjà
  gene_folder <- gene
  if (!dir.exists(gene_folder)) {
    dir.create(gene_folder)
  }
  
  # Sauvegarder la figure dans le dossier du gène
  filename_loc <- paste0(gene_folder, "/cellular_localization_", gene, ".png")
  ggsave(filename = filename_loc, plot = p_loc, width = 10, height = 8)
}
