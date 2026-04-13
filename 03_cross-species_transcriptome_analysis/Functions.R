# Functions.R
# Core functions for PCA, UMAP, and GO Enrichment used across the analytical scripts.

# ==========================================
# 0. Global Palettes
# ==========================================

species_palette <- c(
  "Anatherum_virginicum"     = "#1F77B4", "Andropogon_gerardii"      = "#FF7F0E",
  "Bothriochloa_laguroides"  = "#2CA02C", "Chrysopogon_zizanioides"  = "#D62728",
  "Elionorus_tripsacoides"   = "#9467BD", "Miscanthus_giganteus"     = "#8C564B",
  "Panicum_virgatum"         = "#98DF8A", "Schizachyrium_scoparium"  = "#E377C2",
  "Sorghastrum_nutans"       = "#7F7F7F", "Sorghum_bicolor"          = "#FF9896",
  "Themeda_triandra"         = "#E6AB02", "Tripsacum_dactyloides"    = "#17BECF",
  "Zea_diploperennis"        = "#BCBD22", "Zea_mays"                 = "#36eb11"
)

accession_palette <- c(
  "Anatherum_virginicum"     = "#1F77B4", "Andropogon_gerardii"      = "#FF7F0E",
  "Ag_Niagara"               = "#A6CEE3", "Bothriochloa_laguroides"  = "#2CA02C",
  "Chrysopogon_zizanioides"  = "#D62728", "Elionorus_tripsacoides"   = "#9467BD",
  "Miscanthus_giganteus"     = "#8C564B", "Panicum_virgatum_A1"      = "#98DF8A",
  "Panicum_virgatum_Kanlow"  = "#2CA02C", "Pv_Liberty"               = "#f69400",
  "Schizachyrium_scoparium"  = "#E377C2", "Sorghastrum_nutans"       = "#7F7F7F",
  "Sn_WI"                    = "#993404", "Sorghum_bicolor"          = "#FF9896",
  "Sb_BTx623"                = "#67000D", "Ss_Ahring"                = "#B3DE69",
  "Td_floridanum"            = "#B2DF8A", "Td_IL"                    = "#66C2A5",
  "Td_KS_B6"                 = "#1F78B4", "Td_KS_C9"                 = "#6A3D9A",
  "Td_McKain"                = "#A6CEE3", "Td_MO"                    = "#B2E2E2",
  "Td_OK"                    = "#63A3CB", "Td_Pete"                  = "#0ad2e4",
  "Td_TX"                    = "#084081", "Td_Iuka"                  = "#0f93e4",
  "Td_Texas_Sue"             = "#0cc688", "Themeda_triandra"         = "#E6AB02",
  "Zea_diploperennis"        = "#BCBD22", "Zd_momo"                  = "#006D2C",
  "Zm_lh244"                 = "#36eb11", "Zm_lh287Xlh244"           = "#77e204"
)

# ==========================================
# 1. Dimension Reduction Functions
# ==========================================

plot_pca_ggplot2 <- function(count_data, metadata, color_by, species_palette = NULL, accession_palette = NULL, x_pc = 1, y_pc = 2) {
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  
  # Ensure count data columns match metadata sample IDs
  metadata <- metadata %>% filter(sample_id %in% colnames(count_data))
  count_data <- count_data[, metadata$sample_id, drop = FALSE]
  
  # Check for dimension mismatch
  if (ncol(count_data) != nrow(metadata)) {
    stop("Mismatch between count_data columns and metadata rows after filtering.")
  }
  
  # Ensure metadata order matches count_data columns
  metadata <- metadata[match(colnames(count_data), metadata$sample_id), ]

  # Create DESeqDataSet and normalize
  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = metadata, design = ~ 1)
  dds <- estimateSizeFactors(dds)

  # Filter out low-expressed genes
  dds <- dds[rowSums(counts(dds, normalized = TRUE)) > 0, ] 

  # Perform variance stabilizing transformation
  vst_data <- vst(dds, blind = TRUE)
  vst_mat <- assay(vst_data)

  n_genes_used <- nrow(vst_mat)
  message("Number of genes used in PCA analysis: ", n_genes_used)

  # Perform PCA
  pca_res <- prcomp(t(vst_mat))
  
  # Prepare PCA data frame with metadata
  pca_df <- as.data.frame(pca_res$x)
  pca_df$sample_id <- rownames(pca_df)
  pca_df <- left_join(pca_df, metadata, by = "sample_id")
  
  # Validate x_pc and y_pc inputs
  num_pcs <- ncol(pca_df) - 1
  if (x_pc > num_pcs | y_pc > num_pcs) {
    stop("Requested PCA components exceed available components.")
  }
  
  # Get variance explained
  variance_explained <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 2)
  x_pc_var <- variance_explained[x_pc]
  y_pc_var <- variance_explained[y_pc]
  
  x_pc_col <- paste0("PC", x_pc)
  y_pc_col <- paste0("PC", y_pc)
  
  # Color scaling
  if (is.numeric(pca_df[[color_by]])) {
    color_scale <- scale_color_gradientn(colors = c("#ff991c", "#b200d5","#490057"))
  } else {
    unique_values <- unique(pca_df[[color_by]])
    if (color_by == "species" && !is.null(species_palette)) {
      color_scale <- scale_color_manual(values = species_palette)
    } else if (color_by == "accession" && !is.null(accession_palette)) {
      color_scale <- scale_color_manual(values = accession_palette)
    } else {
      color_palette <- RColorBrewer::brewer.pal(min(length(unique_values), 8), "Set2")
      color_scale <- scale_color_manual(values = setNames(color_palette, unique_values))
    }
  }
  
  # Plot
  pca_plot <- ggplot(pca_df, aes_string(x = x_pc_col, y = y_pc_col, color = color_by)) +
    geom_point(size = 0.8, alpha = 0.8) +
    labs(
      title = paste("PCA of VST Transformed Counts colored by", color_by),
      x = paste(x_pc_col, "-", x_pc_var, "% variance"),
      y = paste(y_pc_col, "-", y_pc_var, "% variance")
    ) +
    color_scale +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "right",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )

  attr(pca_plot, "n_genes_used") <- n_genes_used
  
  return(list(
    pca_plot = pca_plot,
    vst_matrix = vst_mat,
    pca_df = pca_df,
    variance_explained = variance_explained,
    n_genes_used = n_genes_used
  ))
}

run_umap_seurat_style <- function(count_data,
                                  metadata,
                                  color_by,
                                  species_palette = NULL,
                                  accession_palette = NULL,
                                  n_top_genes = 2000,
                                  scale_max = 10,
                                  n_pcs = 10,
                                  umap_n_neighbors = 50,
                                  umap_min_dist = 0.4,
                                  umap_metric = "cosine",
                                  spread = 0.5,
                                  seed = 14) {
  library(DESeq2)
  library(matrixStats)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(uwot)

  set.seed(seed)

  metadata <- metadata %>% filter(sample_id %in% colnames(count_data))
  count_data <- count_data[, metadata$sample_id, drop = FALSE]
  metadata <- metadata[match(colnames(count_data), metadata$sample_id), ]

  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = metadata, design = ~ 1)
  dds <- estimateSizeFactors(dds)
  dds <- dds[rowSums(counts(dds, normalized = TRUE)) > 0, ]
  vst_data <- vst(dds, blind = TRUE)
  vst_mat <- assay(vst_data)

  gene_vars <- rowVars(vst_mat)
  gene_order <- sort(gene_vars, decreasing = TRUE)
  top_genes <- names(gene_order)[1:n_top_genes]
  vst_top <- vst_mat[top_genes, ]

  gene_elbow_df <- data.frame(Rank = seq_along(gene_order), Variance = gene_order)
  gene_elbow_plot <- ggplot(gene_elbow_df, aes(x = Rank, y = Variance)) +
    geom_line() +
    geom_vline(xintercept = n_top_genes, color = "red", linetype = "dashed") +
    labs(title = "Elbow Plot of Gene Variance", x = "Gene Rank", y = "Variance") +
    theme_minimal()

  scaled_top <- t(scale(t(vst_top)))
  scaled_top[scaled_top > scale_max] <- scale_max
  scaled_top[scaled_top < -scale_max] <- -scale_max

  pca_res <- prcomp(t(scaled_top))
  pcs_mat <- pca_res$x[, 1:n_pcs]
  pc_var <- pca_res$sdev^2
  pc_var_explained <- pc_var / sum(pc_var) * 100

  max_pc_plot <- min(length(pc_var), 50)
  pc_elbow_df <- data.frame(PC = 1:max_pc_plot, Variance = pc_var_explained[1:max_pc_plot])
  pc_elbow_plot <- ggplot(pc_elbow_df, aes(x = PC, y = Variance)) +
    geom_line() +
    geom_vline(xintercept = n_pcs, color = "red", linetype = "dashed") +
    labs(title = "Elbow Plot of Principal Components", x = "PC", y = "Variance Explained (%)") +
    theme_minimal()

  umap_res <- umap(pcs_mat, n_neighbors = umap_n_neighbors, min_dist = umap_min_dist,
                   metric = umap_metric, spread = spread)
  umap_df <- as.data.frame(umap_res)
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$sample_id <- rownames(pcs_mat)
  umap_df <- left_join(umap_df, metadata, by = "sample_id")

  if (is.numeric(umap_df[[color_by]])) {
    color_scale <- scale_color_gradientn(colors = c("#ff991c", "#b200d5", "#490057"))
  } else {
    unique_values <- unique(umap_df[[color_by]])
    if (color_by == "species" && !is.null(species_palette)) {
      color_scale <- scale_color_manual(values = species_palette)
    } else if (color_by == "accession" && !is.null(accession_palette)) {
      color_scale <- scale_color_manual(values = accession_palette)
    } else {
      color_palette <- RColorBrewer::brewer.pal(min(length(unique_values), 8), "Set2")
      color_scale <- scale_color_manual(values = setNames(color_palette, unique_values))
    }
  }

  umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = .data[[color_by]])) +
    geom_point(size = 0.025, alpha = 0.8) +
    color_scale +
    coord_fixed() +
    theme_classic(base_size = 10) +
    labs(
      title = paste("Seurat-style UMAP (Top", n_top_genes, "genes, PC", n_pcs, ")"),
      x = "UMAP1", y = "UMAP2"
    )

  return(list(
    umap_plot = umap_plot,
    gene_elbow_plot = gene_elbow_plot,
    pc_elbow_plot = pc_elbow_plot,
    vst_matrix = vst_mat,
    vst_top = vst_top,
    scaled_matrix = scaled_top,
    pca_matrix = pcs_mat,
    umap_df = umap_df
  ))
}

# ==========================================
# 2. Pseudotime & Marker Gene Functions
# ==========================================

calculate_pseudotime <- function(count_data, metadata, covariate, elbo_tol = 0.00001) {
  library(DESeq2)
  library(phenopath)
  
  # Ensure sample IDs match
  count_data <- count_data[, metadata$sample_id, drop = FALSE]
  
  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = metadata, design = ~ 1)
  dds <- estimateSizeFactors(dds)
  
  min_libs <- ceiling(0.10 * ncol(dds))
  keep <- rowSums(counts(dds, normalized = TRUE) > 0) >= min_libs
  dds <- dds[keep, ]
  cat("Keeping", sum(keep), "genes expressed in at least", min_libs, "libraries\n")
  
  vst_data <- vst(dds, blind = TRUE)
  vst_mat <- assay(vst_data)
  
  vst_mat_transposed <- t(vst_mat)
  gene_variances <- apply(vst_mat_transposed, 2, var)
  zero_var_genes <- colnames(vst_mat_transposed)[gene_variances == 0]
  if(length(zero_var_genes) > 0) {
    cat("Removing", length(zero_var_genes), "genes with zero variance.\n")
    vst_mat_transposed <- vst_mat_transposed[, gene_variances > 0]
  }
  
  metadata[[covariate]] <- as.factor(metadata[[covariate]])
  
  fit <- phenopath(vst_mat_transposed, metadata[[covariate]], z_init = 1, elbo_tol = elbo_tol)
  pseudotime <- trajectory(fit)
  pseudotime_df <- data.frame(SampleID = rownames(vst_mat_transposed), Pseudotime = pseudotime)
  
  return(list(
    pseudotime_df = pseudotime_df,
    fit = fit,
    zero_var_genes = zero_var_genes,
    vst_matrix = rbind(vst_mat_transposed, matrix(NA, nrow=0, ncol=ncol(vst_mat_transposed))) # Not used but available
  ))
}

plot_pca_colored_by_gene <- function(gene_id, count_matrix, metadata, x_pc = 1, y_pc = 2) {
  library(DESeq2)
  library(ggplot2)
  library(dplyr)

  metadata <- metadata %>% filter(sample_id %in% colnames(count_matrix))
  count_matrix <- count_matrix[, metadata$sample_id, drop = FALSE]
  metadata <- metadata[match(colnames(count_matrix), metadata$sample_id), ]

  dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = metadata, design = ~1)
  dds <- estimateSizeFactors(dds)
  dds <- dds[rowSums(counts(dds, normalized = TRUE)) > 0, ]

  vst_data <- vst(dds, blind = TRUE)
  vst_mat <- assay(vst_data)

  if (!gene_id %in% rownames(vst_mat)) {
    stop(paste("Gene", gene_id, "not found in vst matrix"))
  }

  pca_res <- prcomp(t(vst_mat))
  pca_df <- as.data.frame(pca_res$x)
  pca_df$sample_id <- rownames(pca_df)
  pca_df <- left_join(pca_df, metadata, by = "sample_id")

  expr_vals <- vst_mat[gene_id, pca_df$sample_id]
  scaled_expr <- scale(expr_vals)
  scaled_expr[scaled_expr > 3] <- 3
  scaled_expr[scaled_expr < -3] <- -3
  pca_df$expression <- as.numeric(scaled_expr)

  var_exp <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 2)
  x_lab <- paste0("PC", x_pc, " (", var_exp[x_pc], "% variance)")
  y_lab <- paste0("PC", y_pc, " (", var_exp[y_pc], "% variance)")

  custom_colors <- rev(c("#7e203a", "#ca1875", "#ef2c69", "#f06848", "#fda25e", "#fceba5", "#f2f8d5", "#e4f2f8"))
  color_scale <- scale_color_gradientn(colors = custom_colors, limits = c(-3, 3))

  ggplot(pca_df, aes_string(x = paste0("PC", x_pc), y = paste0("PC", y_pc), color = "expression")) +
    geom_point(size = 0.8, alpha = 0.9) +
    color_scale +
    theme_classic(base_size = 10) +
    labs(title = paste("PCA colored by", gene_id, "expression (z-scored)"), x = x_lab, y = y_lab, color = "Z-score")
}

export_pca_plot <- function(pca_object, filename, point_size = 1.5, point_alpha = 0.9, shape = 19, base_size = 8, width = 2, height = 2, dpi = 600, compression = "lzw") {
  library(ggplot2)
  pca_object$layers <- pca_object$layers[!sapply(pca_object$layers, function(layer) inherits(layer$geom, "GeomPoint"))]

  pca_object <- pca_object +
    geom_point(size = point_size, alpha = point_alpha, shape = shape, stroke = 0) +
    theme_classic(base_size = base_size) +
    theme(
      plot.title = element_blank(),
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none"
    )
  print(pca_object)
  ggsave(filename = filename, plot = pca_object, width = width, height = height, units = "in", dpi = dpi, compression = compression)
  return(pca_object)
}

# ==========================================
# 3. Photo Index Functions
# ==========================================

calculate_photo_hog_proportion <- function(hog_counts, photo_hogs, metadata) {
  library(DESeq2)
  library(dplyr)
  
  photo_hogs <- intersect(photo_hogs, rownames(hog_counts))

  dds <- DESeqDataSetFromMatrix(countData = hog_counts, colData = metadata, design = ~ 1)
  dds <- estimateSizeFactors(dds)
  normalized_counts <- counts(dds, normalized = TRUE)

  total_counts <- colSums(normalized_counts)
  normalized_photo_hog_counts <- normalized_counts[photo_hogs, , drop = FALSE]
  photo_hog_sums <- colSums(normalized_photo_hog_counts)

  photo_percent_counts <- (photo_hog_sums / total_counts) * 100

  lower_bound <- quantile(photo_percent_counts, 0.05, na.rm = TRUE)
  upper_bound <- quantile(photo_percent_counts, 0.95, na.rm = TRUE)

  photo_percent_df <- data.frame(
    sample_id = colnames(hog_counts),
    photo_percent_counts = photo_percent_counts
  ) %>%
    mutate(ps_index = (photo_percent_counts - lower_bound) / (upper_bound - lower_bound)) %>%
    mutate(ps_index = pmin(pmax(ps_index, 0), 1))

  attr(photo_percent_df, "photo_hogs_used") <- photo_hogs
  attr(photo_percent_df, "n_photo_hogs_used") <- length(photo_hogs)
  attr(photo_percent_df, "n_total_ogs") <- nrow(hog_counts)

  return(photo_percent_df)
}

generate_ps_gdd_plot <- function(df, target_season, week_limits, week_breaks) {
  library(ggplot2)
  library(dplyr)
  
  scale_fact <- 100 
  season_df <- df %>% filter(season == target_season)
  
  ggplot(season_df, aes(x = week)) +
    geom_jitter(aes(y = ps_index), color = "green4", alpha = 0.3, width = 0.25, size = 0.8, stroke = 0.3) +
    stat_summary(aes(y = ps_index, group = 1), fun = mean, geom = "line", color = "green4", linewidth = 1) +
    stat_summary(aes(y = weekly_gdd_pulse / scale_fact, group = 1), fun = mean, geom = "line", color = "#8B4513", linewidth = 1) +
    stat_summary(aes(y = ps_index), fun = mean, geom = "point", shape = 21, color = "black", fill = "green3", size = 2.5, stroke = 0.8) +
    stat_summary(aes(y = weekly_gdd_pulse / scale_fact), fun = mean, geom = "point", shape = 21, color = "black", fill = "#d95f02", size = 2.5, stroke = 0.8) +
    scale_y_continuous(
      name = "photosynthetic index",
      limits = c(0, 1.1),
      breaks = seq(0, 1, 0.2),
      sec.axis = sec_axis(~ . * scale_fact, name = expression(paste("accumulated weekly GDDs (°C)")), breaks = seq(0, 100, by = 25))
    ) +
    scale_x_continuous(limits = week_limits, breaks = week_breaks) +
    theme_classic(base_size = 10) +
    theme(
      axis.title.y.left = element_text(color = "black"),
      axis.title.y.right = element_text(color = "black"),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      plot.title = element_text(face = "bold.italic")
    ) +
    labs(x = "week")
}

# ==========================================
# 4. GO Enrichment Functions
# ==========================================

perform_GO_enrichment_OG <- function(og_vector, gmt_bp_df, gmt_mf_df, gmt_cc_df, gmt_data_df) {
  library(clusterProfiler)
  
  if (!is.character(og_vector)) {
    stop("Input OGs must be a character vector")
  }
  
  perform_enrichment <- function(og_vector, gmt_df, ontology_label, gmt_data_df) {
    enrichment_results <- enricher(og_vector, TERM2GENE = gmt_df, pAdjustMethod = "BH")
    if (is.null(enrichment_results) || nrow(enrichment_results@result) == 0) {
      return(data.frame())
    }
    enrichment_results_with_desc <- merge(enrichment_results@result, gmt_data_df, 
                                          by.x = "ID", by.y = "goid", all.x = TRUE)
    enrichment_results_with_desc$Ontology <- ontology_label
    return(enrichment_results_with_desc)
  }
  
  bp_results <- perform_enrichment(og_vector, gmt_bp_df, "BP", gmt_data_df)
  mf_results <- perform_enrichment(og_vector, gmt_mf_df, "MF", gmt_data_df)
  cc_results <- perform_enrichment(og_vector, gmt_cc_df, "CC", gmt_data_df)
  
  go_results <- rbind(bp_results, mf_results, cc_results)
  if (!is.null(go_results) && nrow(go_results) > 0) {
    go_results$p.adjust <- p.adjust(go_results$pvalue, method = "BH")
  }
  return(go_results)
}

perform_and_plot_GO_OG <- function(og_vector, gmt_bp_df, gmt_mf_df, gmt_cc_df, gmt_data_df, output_txt) {
  library(ggplot2)
  library(dplyr)
  
  go_results <- perform_GO_enrichment_OG(og_vector, gmt_bp_df, gmt_mf_df, gmt_cc_df, gmt_data_df)
  
  if (nrow(go_results) > 0) {
    write.table(go_results, file = output_txt, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("GO enrichment results saved to:", output_txt, "\n")
  }
  
  if (nrow(go_results) == 0) return(list(go_results = go_results, lollipop_plot = NULL, miniplot = NULL))
  
  go_data <- go_results %>%
    mutate(
      p.adjust = as.numeric(p.adjust),
      logP = -log10(p.adjust),
      short_desc = sapply(strsplit(desc, " "), function(x) paste(head(x, 6), collapse = " ")),
      label = paste0(short_desc, " (", ID, ")")
    ) %>%
    filter(!is.na(logP) & is.finite(logP))

  top_by_ontology <- go_data %>%
    group_by(Ontology) %>%
    slice_max(order_by = logP, n = 6) %>%
    ungroup() %>%
    mutate(highlight = TRUE)

  go_data <- go_data %>%
    mutate(
      highlight = ifelse(ID %in% top_by_ontology$ID, "highlight", 
                         ifelse(p.adjust < 0.05, "pale", "nonsignificant"))
    )
  
  ontology_colors <- c("BP" = "orange", "MF" = "brown1", "CC" = "#00bfc4")
  
  p_lollipop <- ggplot(top_by_ontology, aes(x = logP, y = reorder(label, logP), color = Ontology)) +
    geom_segment(aes(xend = 0, yend = label), size = 0.8) +
    geom_point(aes(size = Count), alpha = 0.8) +
    scale_size_continuous(trans = "identity", breaks = pretty(top_by_ontology$Count, n = 5)) + 
    scale_color_manual(values = ontology_colors) + 
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "darkslategray") +
    labs(x = "-log10(Adjusted p-value)", y = "GO Term (Description + GO ID)", size = "OG Count", color = "Ontology") +
    theme_classic() +
    theme(
      axis.text.y = element_text(size = 8),  
      axis.title.x = element_text(size = 8),
      legend.position = "right",
      legend.text = element_text(size = 8),  
      legend.title = element_text(size = 8),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  pale_colors <- c("BP" = "#FFD9B3", "MF" = "#E8B2B2", "CC" = "#AEECEF")
  highlight_colors <- ontology_colors  
  nonsignificant_color <- "gray80"  

  go_data$color <- ifelse(go_data$highlight == "highlight", highlight_colors[go_data$Ontology], 
                          ifelse(go_data$highlight == "pale", pale_colors[go_data$Ontology], 
                                 nonsignificant_color))

  p_bubble <- ggplot(go_data, aes(x = logP, y = reorder(label, logP), color = color)) +
    geom_point(alpha = 0.8, size = 1) +  
    scale_color_identity() +  
    facet_wrap(~Ontology, scales = "free_y", ncol = 1) +  
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "darkslategrey", linewidth = 0.8) +  
    scale_x_continuous(limits = c(0, max(go_data$logP, na.rm = TRUE) * 1.1)) +  
    labs(x = "-log10(Adjusted p-value)", y = NULL) +
    theme_classic() +
    theme(
        strip.text = element_text(size = 8, face = "bold"),  
        axis.text.y = element_blank(),  
        axis.title.y = element_blank(),  
        axis.ticks.y = element_blank(),  
        axis.text.x = element_text(size = 8),  
        legend.position = "none",
        plot.margin = margin(5, 5, 5, 5)
      )

  return(list(go_results = go_results, lollipop_plot = p_lollipop, miniplot = p_bubble))
}

perform_GO_enrichment_OG_specific_universe <- function(og_vector, gmt_bp_df, gmt_mf_df, gmt_cc_df, gmt_data_df, universe_ogs, minGS = 5, maxGS = 2000, pAdjustMethod = "BH") {
  library(clusterProfiler)
  
  stopifnot(is.character(og_vector))
  og_vector <- intersect(og_vector, universe_ogs)
  if (length(og_vector) < minGS) return(data.frame())

  perform_enrichment <- function(og_vec, gmt_df, ontology_label) {
    gmt_df_u <- gmt_df %>%
      dplyr::filter(gene %in% universe_ogs) %>%
      dplyr::distinct(term, gene)

    if (nrow(gmt_df_u) == 0) return(data.frame())

    res <- clusterProfiler::enricher(
      gene          = og_vec,
      TERM2GENE     = gmt_df_u,
      universe      = universe_ogs,
      pAdjustMethod = pAdjustMethod,
      minGSSize     = minGS,
      maxGSSize     = maxGS
    )
    if (is.null(res) || nrow(res@result) == 0) return(data.frame())

    out <- res@result
    out <- out %>%
      dplyr::left_join(gmt_data_df %>% dplyr::select(goid = goid, desc),
                       by = c("ID" = "goid")) %>%
      dplyr::mutate(
        Ontology   = ontology_label,
        logP       = -log10(p.adjust)
      )
    out
  }

  bp_results <- perform_enrichment(og_vector, gmt_bp_df, "BP")
  mf_results <- perform_enrichment(og_vector, gmt_mf_df, "MF")
  cc_results <- perform_enrichment(og_vector, gmt_cc_df, "CC")

  go_results <- dplyr::bind_rows(bp_results, mf_results, cc_results)
  if (nrow(go_results)) go_results$p.adjust <- p.adjust(go_results$pvalue, method = pAdjustMethod)
  go_results
}
