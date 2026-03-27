# 03_plot_UMAP.R
# This script computes and visualizes a Seurat-style UMAP of all libraries
# across both seasons based on the consolidated Objective 03 data objects.

# Command-line arguments for input/output directories
args <- commandArgs(trailingOnly = TRUE)
in_dir <- if (length(args) > 0) args[1] else file.path(getwd(), "03_cross-species_transcriptome_analysis", "data")
out_dir <- if (length(args) > 1) args[2] else file.path(getwd(), "03_cross-species_transcriptome_analysis", "plots")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

library(DESeq2)
library(matrixStats)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(uwot)

# Load consolidated data objects
metadata <- readRDS(file.path(in_dir, "metadata.rds"))
og_counts <- readRDS(file.path(in_dir, "og_counts.rds"))

# Define UMAP function
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
  set.seed(seed)

  # Match metadata and counts
  metadata <- metadata %>% filter(sample_id %in% colnames(count_data))
  count_data <- count_data[, metadata$sample_id, drop = FALSE]
  metadata <- metadata[match(colnames(count_data), metadata$sample_id), ]

  # VST transform
  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = metadata, design = ~ 1)
  dds <- estimateSizeFactors(dds)
  dds <- dds[rowSums(counts(dds, normalized = TRUE)) > 0, ]
  vst_data <- vst(dds, blind = TRUE)
  vst_mat <- assay(vst_data)

  # Gene variance and elbow plot
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

  # Seurat-style scaling
  scaled_top <- t(scale(t(vst_top)))
  scaled_top[scaled_top > scale_max] <- scale_max
  scaled_top[scaled_top < -scale_max] <- -scale_max

  # PCA on scaled data
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

  # UMAP
  umap_res <- umap(pcs_mat, n_neighbors = umap_n_neighbors, min_dist = umap_min_dist,
                   metric = umap_metric, spread = spread)
  umap_df <- as.data.frame(umap_res)
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$sample_id <- rownames(pcs_mat)
  umap_df <- left_join(umap_df, metadata, by = "sample_id")

  # Color scale
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

  # UMAP plot
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

# --- Define Custom Palettes ---
species_palette <- c(
  "Anatherum_virginicum"     = "#1F77B4",
  "Andropogon_gerardii"      = "#FF7F0E",
  "Bothriochloa_laguroides"  = "#2CA02C",
  "Chrysopogon_zizanioides"  = "#D62728",
  "Elionorus_tripsacoides"   = "#9467BD",
  "Miscanthus_giganteus"     = "#8C564B",
  "Panicum_virgatum"         = "#98DF8A",
  "Schizachyrium_scoparium"  = "#E377C2",
  "Sorghastrum_nutans"       = "#7F7F7F",
  "Sorghum_bicolor"          = "#FF9896",
  "Themeda_triandra"         = "#E6AB02",
  "Tripsacum_dactyloides"    = "#17BECF",
  "Zea_diploperennis"        = "#BCBD22",
  "Zea_mays"                 = "#36eb11"
)

accession_palette_named <- c(
  "Anatherum_virginicum"     = "#1F77B4",
  "Andropogon_gerardii"      = "#FF7F0E",
  "Ag_Niagara"               = "#A6CEE3",
  "Bothriochloa_laguroides"  = "#2CA02C",
  "Chrysopogon_zizanioides"  = "#D62728",
  "Elionorus_tripsacoides"   = "#9467BD",
  "Miscanthus_giganteus"     = "#8C564B",
  "Panicum_virgatum_A1"      = "#98DF8A",
  "Panicum_virgatum_Kanlow"  = "#2CA02C",
  "Pv_Liberty"               = "#f69400",
  "Schizachyrium_scoparium"  = "#E377C2",
  "Sorghastrum_nutans"       = "#7F7F7F",
  "Sn_WI"                    = "#993404",
  "Sorghum_bicolor"          = "#FF9896",
  "Sb_BTx623"                = "#67000D",
  "Ss_Ahring"                = "#B3DE69",
  "Td_floridanum"            = "#B2DF8A",
  "Td_IL"                    = "#66C2A5",
  "Td_KS_B6"                 = "#1F78B4",
  "Td_KS_C9"                 = "#6A3D9A",
  "Td_McKain"                = "#A6CEE3",
  "Td_MO"                    = "#B2E2E2",
  "Td_OK"                    = "#63A3CB",
  "Td_Pete"                  = "#0ad2e4",
  "Td_TX"                    = "#084081",
  "Td_Iuka"                  = "#0f93e4",
  "Td_Texas_Sue"             = "#0cc688",
  "Themeda_triandra"         = "#E6AB02",
  "Zea_diploperennis"        = "#BCBD22",
  "Zd_momo"                  = "#006D2C",
  "Zm_lh244"                 = "#36eb11",
  "Zm_lh287Xlh244"           = "#77e204"
)

# Run the UMAP process
result <- run_umap_seurat_style(
  count_data = og_counts,
  metadata = metadata,
  color_by = "tissue",
  species_palette = species_palette,
  accession_palette = accession_palette_named,
  n_top_genes = 3000,
  n_pcs = 6
)

# View diagnostic plots (optional for interactive runs)
# print(result$gene_elbow_plot)
# print(result$pc_elbow_plot)

# --- 1. Base UMAP Plot (No explicitly assigned palette) ---
umap_species <- result$umap_plot +
  geom_point(size = 0.025, alpha = 0.9, stroke = 0) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 14),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )

tiff(file.path(out_dir, "umap_test.tiff"), width = 4, height = 4, units = "in", res = 600, compression = "lzw")
print(umap_species)
dev.off()

# --- 2. UMAP colored by tissue_type palette ---
tissue_type_palette <- c(
  "leaf" = "#30e004",    # Green
  "root" = "#ee6766",    # Purple
  "stem" = "#07cde3"     # Blue
)

umap_tissue_type <- result$umap_plot +
  scale_color_manual(values = tissue_type_palette, name = "Tissue Type") +
  geom_point(size = 0.025, alpha = 0.9, stroke = 0) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 14),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )

tiff(file.path(out_dir, "umap_tissue_type.tiff"), width = 2.5, height = 2.5, units = "in", res = 600, compression = "lzw")
print(umap_tissue_type)
dev.off()

# --- 3. UMAP plot styled by time (default plot color) ---
umap_time <- result$umap_plot +
  geom_point(size = 0.025, alpha = 0.9) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )

tiff(file.path(out_dir, "umap_time.tiff"), width = 2.5, height = 2.5, units = "in", res = 600)
print(umap_time)
dev.off()

# --- 4. UMAP colored by simpler tissue palette ---
tissue_palette <- c(
  "leaf" = "#3BA935",    # green
  "stem" = "#2C82C9",    # blue
  "root" = "#C6A969"     # light brown
)

umap_tissue <- result$umap_plot +
  scale_color_manual(values = tissue_palette) +
  geom_point(size = 0.025, alpha = 0.9) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )

tiff(file.path(out_dir, "umap_tissue.tiff"), width = 2.5, height = 2.5, units = "in", res = 600)
print(umap_tissue)
dev.off()

# --- 5. UMAP colored by season (needs re-running or manually coloring the umap_df data layer) ---
# Wait, the original data was generated with color_by="tissue". If we remap colors to "season", 
# we must use aes(color=season) because the static plot structure uses tissue.
umap_season <- ggplot(result$umap_df, aes(x = UMAP1, y = UMAP2, color = season)) +
  geom_point(size = 0.025, alpha = 0.5) +
  scale_color_manual(values = c("2023" = "#4B0082", "2024" = "#DAA520")) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )

tiff(file.path(out_dir, "umap_season_label.tiff"), width = 2.5, height = 2.5, units = "in", res = 600)
print(umap_season)
dev.off()

message("All Seurat-style UMAP figures generated and saved to ", out_dir)
