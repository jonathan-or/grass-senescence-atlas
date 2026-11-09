# 07_qc_all_libraries.R
# Supplementary Figure 4: library QC on all 3,798 sequenced libraries.
#   A-D  PCA (PC1 vs PC2) colored by tissue, sampling week, species, season
#   E-I  PCA colored by mapping rate, mapped reads, and the QC cutoffs (35% mapping, 100K reads)
#   J    stacked bars of libraries that passed / failed QC by species, tissue, season, week
#
# Inputs (03_cross-species_transcriptome_analysis/data):
#   metadata_all_libraries.rds, og_counts_all_libraries.rds  (from 06_merge_all_libraries.R)
#   metadata.rds                                               (the 2,685 libraries that passed QC)
# The VST and PCA of all libraries are cached in data/ (vst_all_libraries.rds,
# pca_all_libraries.rds); delete them to recompute.
#
# Usage (from the repository root):
#   Rscript 03_cross-species_transcriptome_analysis/07_qc_all_libraries.R [in_dir] [out_dir]

args <- commandArgs(trailingOnly = TRUE)
in_dir  <- if (length(args) > 0) args[1] else file.path(getwd(), "03_cross-species_transcriptome_analysis", "data")
out_dir <- if (length(args) > 1) args[2] else file.path(getwd(), "03_cross-species_transcriptome_analysis", "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(DESeq2)
  library(patchwork)
  library(RColorBrewer)
  library(forcats)
})

meta_whole  <- readRDS(file.path(in_dir, "metadata_all_libraries.rds"))
counts_whole <- readRDS(file.path(in_dir, "og_counts_all_libraries.rds"))
meta_filt   <- readRDS(file.path(in_dir, "metadata.rds"))

# Helper to get VST matrix (from cache or calculate)
get_vst_matrix <- function(count_data, metadata, cache_path) {
  if (file.exists(cache_path)) {
    cat("Loading VST matrix from cache:", cache_path, "\n")
    return(readRDS(cache_path))
  }
  
  cat("Calculating VST matrix (cold run)...\n")
  # Filter out samples with 0 total counts
  valid_samples <- colnames(count_data)[colSums(count_data) > 0]
  count_data <- count_data[, valid_samples, drop = FALSE]
  metadata <- metadata %>% filter(sample_id %in% valid_samples)
  metadata <- metadata[match(colnames(count_data), metadata$sample_id), ]
  
  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = metadata, design = ~ 1)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- dds[rowSums(counts(dds, normalized = TRUE)) > 0, ] 
  vst_data <- vst(dds, blind = TRUE)
  vst_mat <- assay(vst_data)
  
  saveRDS(vst_mat, cache_path)
  return(vst_mat)
}

# Helper to get PCA dataframe (from cache or calculate using VST)
get_pca_df <- function(vst_mat, metadata, cache_path) {
  if (file.exists(cache_path)) {
    cat("Loading PCA coordinates from cache:", cache_path, "\n")
    return(readRDS(cache_path))
  }
  
  cat("Running PCA calculation (cold run)...\n")
  metadata <- metadata %>% filter(sample_id %in% colnames(vst_mat))
  metadata <- metadata[match(colnames(vst_mat), metadata$sample_id), ]
  
  pca_res <- prcomp(t(vst_mat))
  pca_df <- as.data.frame(pca_res$x)
  pca_df$sample_id <- rownames(pca_df)
  pca_df <- left_join(pca_df, metadata, by = "sample_id")
  
  # Attach variance explained as attributes
  variance_explained <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 2)
  attr(pca_df, "variance_explained") <- variance_explained
  
  saveRDS(pca_df, cache_path)
  return(pca_df)
}

# Helper to plot PCA using pre-computed PCA dataframe

# ==========================================
# VST and PCA of all libraries (cached)
# ==========================================
vst_whole  <- get_vst_matrix(counts_whole, meta_whole, file.path(in_dir, "vst_all_libraries.rds"))
pca_df_all <- get_pca_df(vst_whole, meta_whole, file.path(in_dir, "pca_all_libraries.rds"))

# ==========================================
# Panels A-D
# ==========================================
pca_df <- pca_df_all
# Clean labels
pca_df <- pca_df %>%
  mutate(
    species = gsub("_", " ", species),
    season = factor(season),
    week_num = as.numeric(week) # Convert week to numeric for continuous viridis gradient
  )

# --- Define User-Requested Color Palettes ---
# RColorBrewer Set2 colors
set2_colors <- brewer.pal(3, "Set2")
season_palette <- c(
  "2024" = set2_colors[2], # Muted green/teal
  "2023" = set2_colors[1]  # Muted orange/peach
)

tissue_palette <- c(
  "leaf" = "#30e004",    # Green
  "root" = "#ee6766",    # Soft red/pink
  "stem" = "#07cde3"     # Soft blue
)

species_palette <- c(
  "Anatherum virginicum"     = "#1F77B4",
  "Andropogon gerardii"      = "#FF7F0E",
  "Bothriochloa laguroides"  = "#2CA02C",
  "Chrysopogon zizanioides"  = "#D62728",
  "Elionorus tripsacoides"   = "#9467BD",
  "Miscanthus giganteus"     = "#8C564B",
  "Panicum virgatum"         = "#98DF8A",
  "Schizachyrium scoparium"  = "#E377C2",
  "Sorghastrum nutans"       = "#7F7F7F",
  "Sorghum bicolor"          = "#FF9896",
  "Themeda triandra"         = "#E6AB02",
  "Tripsacum dactyloides"    = "#17BECF",
  "Zea diploperennis"        = "#BCBD22",
  "Zea mays"                 = "#36eb11"
)

# Helper plotting function for PCA (both X and Y axes restored, base theme 8)
plot_pca_horizontal <- function(df, color_by, palette = NULL, is_continuous = FALSE, x_lims = NULL, y_lims = NULL) {
  variance_explained <- attr(df, "variance_explained")
  x_pc_var <- variance_explained[1] # PC1
  y_pc_var <- variance_explained[2] # PC2
  
  if (is_continuous) {
    color_scale <- scale_color_viridis_c(option = "viridis", name = "Sampling Week")
  } else {
    name_legend <- case_when(
      color_by == "tissue" ~ "Tissue",
      color_by == "species" ~ "Species",
      color_by == "season" ~ "Season",
      TRUE ~ color_by
    )
    color_scale <- scale_color_manual(values = palette, name = name_legend)
  }
  
  pca_plot <- ggplot(df, aes(x = PC1, y = PC2, color = .data[[color_by]])) +
    geom_point(size = 0.2, alpha = 0.7) +
    color_scale +
    theme_classic(base_size = 8) +  # Base theme size 8
    labs(
      x = paste0("PC1 - ", x_pc_var, "% variance"),
      y = paste0("PC2 - ", y_pc_var, "% variance")
    ) +
    theme(
      # X-Axis (PC1) - Restored with line, ticks, and numbers
      axis.title.x = element_text(size = 6, face = "bold"),
      axis.text.x = element_text(size = 5),
      axis.ticks.x = element_line(linewidth = 0.3),
      axis.line.x = element_line(linewidth = 0.3),
      
      # Y-Axis (PC2) - Restored with line, ticks, and numbers
      axis.title.y = element_text(size = 6, face = "bold"),
      axis.text.y = element_text(size = 5),
      axis.ticks.y = element_line(linewidth = 0.3),
      axis.line.y = element_line(linewidth = 0.3),
      
      plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      
      # Style the legends slightly smaller so they fit nicely
      legend.key.size = unit(0.2, "cm"),
      legend.text = element_text(size = 5),
      legend.title = element_text(size = 6, face = "bold")
    )
  
  if (!is.null(x_lims)) pca_plot <- pca_plot + xlim(x_lims)
  if (!is.null(y_lims)) pca_plot <- pca_plot + ylim(y_lims)
  
  return(pca_plot)
}

# Get exact PCA limits to align panels
x_lims <- range(pca_df$PC1, na.rm = TRUE)
y_lims <- range(pca_df$PC2, na.rm = TRUE)

cat("Generating individual PCA plots...\n")
p_tissue <- plot_pca_horizontal(pca_df, "tissue", tissue_palette, x_lims = x_lims, y_lims = y_lims) +
  labs(title = "Tissue")

p_week <- plot_pca_horizontal(pca_df, "week_num", is_continuous = TRUE, x_lims = x_lims, y_lims = y_lims) +
  labs(title = "Sampling Week")

p_species <- plot_pca_horizontal(pca_df, "species", species_palette, x_lims = x_lims, y_lims = y_lims) +
  labs(title = "Species")

p_season <- plot_pca_horizontal(pca_df, "season", season_palette, x_lims = x_lims, y_lims = y_lims) +
  labs(title = "Season")

# ==========================================
# 1. Save Clean Version (No Legends, 8x2.5)
# ==========================================
cat("Saving clean version (no legends, 8x2.5)...\n")
combined_clean <- (p_tissue | p_week | p_species | p_season) & theme(legend.position = "none")

ggsave(file.path(out_dir, "pca_horizontal_whole.png"), combined_clean, width = 8, height = 2.5, dpi = 300)
tiff(file.path(out_dir, "pca_horizontal_whole.tiff"), width = 8, height = 2.5, units = "in", res = 600)
print(combined_clean)
dev.off()

# ==========================================
# 2. Save Labeled Version (With Legends on Right, 11x2.5)
# ==========================================
cat("Saving labeled version (with legends on right, 11x2.5)...\n")
combined_legend <- (p_tissue | p_week | p_species | p_season) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(out_dir, "pca_horizontal_whole_legend.png"), combined_legend, width = 11, height = 2.5, dpi = 300)
tiff(file.path(out_dir, "pca_horizontal_whole_legend.tiff"), width = 11, height = 2.5, units = "in", res = 600)
print(combined_legend)
dev.off()

cat("Horizontal PCA plot completed successfully!\n")

# ==========================================
# Panels E-I
# ==========================================
pca_df <- pca_df_all
# Calculate mapping rate and QC fields based on CleanReads
pca_df <- pca_df %>%
  mutate(
    mapping_rate = MappedReads / CleanReads,
    mapping_rate = ifelse(is.na(mapping_rate) | is.nan(mapping_rate), 0, mapping_rate),
    
    # Binary variables
    read_cutoff = factor(ifelse(MappedReads >= 100000, "Passed (>= 100K)", "Failed (< 100K)"),
                         levels = c("Passed (>= 100K)", "Failed (< 100K)")),
    mapping_cutoff = factor(ifelse(mapping_rate >= 0.35, "Passed (>= 35%)", "Failed (< 35%)"),
                            levels = c("Passed (>= 35%)", "Failed (< 35%)")),
    overall_qc = factor(ifelse(mapping_rate >= 0.35 & MappedReads >= 100000, "Passed QC", "Failed QC"),
                        levels = c("Passed QC", "Failed QC"))
  )

# --- Define Color Palettes ---
# Soft colorblind-safe QC palette
qc_palette <- c(
  "Passed (>= 100K)" = "#0072B2", "Failed (< 100K)" = "#E69F00",
  "Passed (>= 35%)"  = "#0072B2", "Failed (< 35%)"  = "#E69F00",
  "Passed QC"        = "#0072B2", "Failed QC"        = "#E69F00"
)

# Helper plotting function for PCA (both X and Y axes restored, base theme 8)
plot_pca_qc_horizontal <- function(df, color_by, palette = NULL, is_continuous = FALSE, is_log10 = FALSE, x_lims = NULL, y_lims = NULL) {
  variance_explained <- attr(df, "variance_explained")
  x_pc_var <- variance_explained[1] # PC1
  y_pc_var <- variance_explained[2] # PC2
  
  if (is_continuous) {
    name_leg <- ifelse(color_by == "mapping_rate", "Mapping Rate", "Mapped Reads")
    if (is_log10) {
      color_scale <- scale_color_viridis_c(option = "rocket", direction = -1, trans = "log10", name = name_leg, labels = scales::label_comma())
    } else {
      color_scale <- scale_color_viridis_c(option = "rocket", direction = -1, name = name_leg)
    }
  } else {
    name_legend <- case_when(
      color_by == "read_cutoff" ~ "100K Reads",
      color_by == "mapping_cutoff" ~ "35% Mapping",
      color_by == "overall_qc" ~ "QC Status",
      TRUE ~ color_by
    )
    color_scale <- scale_color_manual(values = palette, name = name_legend)
  }
  
  pca_plot <- ggplot(df, aes(x = PC1, y = PC2, color = .data[[color_by]])) +
    geom_point(size = 0.2, alpha = 0.7) +
    color_scale +
    theme_classic(base_size = 8) +  # Base theme size 8
    labs(
      x = paste0("PC1 - ", x_pc_var, "% variance"),
      y = paste0("PC2 - ", y_pc_var, "% variance")
    ) +
    theme(
      # X-Axis (PC1) - Restored with line, ticks, and size 8 titles
      axis.title.x = element_text(size = 8, face = "bold"),
      axis.text.x = element_text(size = 7),
      axis.ticks.x = element_line(linewidth = 0.3),
      axis.line.x = element_line(linewidth = 0.3),
      
      # Y-Axis (PC2) - Restored with line, ticks, and size 8 titles
      axis.title.y = element_text(size = 8, face = "bold"),
      axis.text.y = element_text(size = 7),
      axis.ticks.y = element_line(linewidth = 0.3),
      axis.line.y = element_line(linewidth = 0.3),
      
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      
      # Style the legends
      legend.key.size = unit(0.25, "cm"),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 7, face = "bold")
    )
  
  if (!is.null(x_lims)) pca_plot <- pca_plot + xlim(x_lims)
  if (!is.null(y_lims)) pca_plot <- pca_plot + ylim(y_lims)
  
  return(pca_plot)
}

# Get exact PCA limits to align panels
x_lims <- range(pca_df$PC1, na.rm = TRUE)
y_lims <- range(pca_df$PC2, na.rm = TRUE)

cat("Generating individual 5-panel PCA plots...\n")
p_map_rate <- plot_pca_qc_horizontal(pca_df, "mapping_rate", is_continuous = TRUE, is_log10 = FALSE, x_lims = x_lims, y_lims = y_lims) +
  labs(title = "Mapping Rate")

p_mapped_reads <- plot_pca_qc_horizontal(pca_df, "MappedReads", is_continuous = TRUE, is_log10 = TRUE, x_lims = x_lims, y_lims = y_lims) +
  labs(title = "Mapped Reads Count")

p_map_cut <- plot_pca_qc_horizontal(pca_df, "mapping_cutoff", qc_palette, x_lims = x_lims, y_lims = y_lims) +
  labs(title = "Mapping Cutoff (35%)")

p_reads <- plot_pca_qc_horizontal(pca_df, "read_cutoff", qc_palette, x_lims = x_lims, y_lims = y_lims) +
  labs(title = "Reads Cutoff (100K)")

p_overall <- plot_pca_qc_horizontal(pca_df, "overall_qc", qc_palette, x_lims = x_lims, y_lims = y_lims) +
  labs(title = "Overall QC Status")

# ==========================================
# 1. Save Clean Version (No Legends, 10x2.5)
# ==========================================
cat("Saving clean 5-panel version (no legends, 10x2.5)...\n")
combined_clean <- (p_map_rate | p_mapped_reads | p_map_cut | p_reads | p_overall) & theme(legend.position = "none")

ggsave(file.path(out_dir, "pca_qc_horizontal_whole.png"), combined_clean, width = 10, height = 2.5, dpi = 300)
tiff(file.path(out_dir, "pca_qc_horizontal_whole.tiff"), width = 10, height = 2.5, units = "in", res = 600)
print(combined_clean)
dev.off()

# ==========================================
# 2. Save Labeled Version (With Legends on Right, 13.5x2.5)
# ==========================================
cat("Saving labeled 5-panel version (with legends on right, 13.5x2.5)...\n")
combined_legend <- (p_map_rate | p_mapped_reads | p_map_cut | p_reads | p_overall) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(out_dir, "pca_qc_horizontal_whole_legend.png"), combined_legend, width = 13.5, height = 2.5, dpi = 300)
tiff(file.path(out_dir, "pca_qc_horizontal_whole_legend.tiff"), width = 13.5, height = 2.5, units = "in", res = 600)
print(combined_legend)
dev.off()

cat("5-Panel Horizontal PCA QC plot completed successfully!\n")

# ==========================================
# Panel J
# ==========================================
# Label kept vs removed
meta_whole <- meta_whole %>%
  mutate(
    status = factor(ifelse(sample_id %in% meta_filt$sample_id, "Kept (Passed QC)", "Removed (Failed QC)"),
                    levels = c("Kept (Passed QC)", "Removed (Failed QC)")),
    # Clean species labels
    species = gsub("_", " ", species),
    # Ensure week is a factor
    week = factor(week)
  )

# Color scheme matching Passed/Kept (blue) and Failed/Removed (orange/yellow)
qc_palette <- c("Kept (Passed QC)" = "#0072B2", "Removed (Failed QC)" = "#E69F00")

# Helper function to generate stacked barplot
make_summary_barplot <- function(df, group_var, fill_palette) {
  # Calculate counts and percentages
  counts_df <- df %>%
    dplyr::count(.data[[group_var]], status, name = "n") %>%
    group_by(.data[[group_var]]) %>%
    mutate(
      percent = n / sum(n),
      label = as.character(n)
    ) %>%
    ungroup()
  
  # For ordering the Y axis by total libraries
  total_counts <- df %>%
    dplyr::count(.data[[group_var]], name = "total")
  
  counts_df <- counts_df %>%
    left_join(total_counts, by = group_var)
  
  # Sort factors by total counts (especially useful for species and tissue)
  if (group_var %in% c("species", "tissue")) {
    counts_df[[group_var]] <- fct_reorder(counts_df[[group_var]], counts_df$total)
  }
  
  plot_out <- ggplot(counts_df, aes(y = .data[[group_var]], x = percent, fill = status)) +
    geom_col(width = 0.8) +
    # Dynamic text coloring: white inside blue (Kept), black inside yellow (Removed)
    geom_text(
      aes(label = label, color = status),
      position = position_stack(vjust = 0.5),
      size = 2.0,
      fontface = "bold"
    ) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
    scale_fill_manual(values = fill_palette) +
    scale_color_manual(values = c("Kept (Passed QC)" = "white", "Removed (Failed QC)" = "black"), guide = "none") +
    theme_classic(base_size = 8) +  # Base theme of 8
    labs(
      x = NULL,
      y = NULL
    ) +
    theme(
      axis.text.y = element_text(size = 7),
      axis.text.x = element_text(size = 6),
      plot.title = element_text(face = "bold", size = 8, hjust = 0.5),
      legend.position = "none"
    )
  
  return(plot_out)
}

cat("Generating summary plots...\n")
p_species <- make_summary_barplot(meta_whole, "species", qc_palette) +
  labs(title = "QC Status by Species")

p_tissue <- make_summary_barplot(meta_whole, "tissue", qc_palette) +
  labs(title = "QC Status by Tissue")

p_season <- make_summary_barplot(meta_whole, "season", qc_palette) +
  labs(title = "QC Status by Season")

p_week <- make_summary_barplot(meta_whole, "week", qc_palette) +
  labs(title = "QC Status by Week")

# Combine side-by-side into a 1x4 horizontal grid (4 columns) without any titles/subtitles/legends
combined_layout <- p_species | p_tissue | p_season | p_week

# Hide all legends across the layout
combined_layout <- combined_layout & theme(legend.position = "none")

cat("Saving plots (8x3 inches)...\n")
ggsave(file.path(out_dir, "library_qc_summary_barplots.png"), combined_layout, width = 8, height = 3, dpi = 300)
tiff(file.path(out_dir, "library_qc_summary_barplots.tiff"), width = 8, height = 3, units = "in", res = 600)
print(combined_layout)
dev.off()

cat("Barplot summary completed successfully!\n")
