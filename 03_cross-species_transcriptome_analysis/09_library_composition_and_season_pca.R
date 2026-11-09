# 09_library_composition_and_season_pca.R
# Supplementary Figure 6: composition of the atlas (2,685 libraries) and seasonal batch effects.
#   A    libraries per species
#   B-C  proportion of libraries per tissue (B) and per season (C), species order as in A
#   D-E  PCA of leaf (D) and root (E) libraries colored by season: PC1 vs PC2, PC3, PC4
# Panels A-B follow end_of_season/02_plot_PCA.R.
#
# Usage (from the repository root):
#   Rscript 03_cross-species_transcriptome_analysis/09_library_composition_and_season_pca.R [in_dir] [out_dir]

args <- commandArgs(trailingOnly = TRUE)
in_dir  <- if (length(args) > 0) args[1] else file.path(getwd(), "03_cross-species_transcriptome_analysis", "data")
out_dir <- if (length(args) > 1) args[2] else file.path(getwd(), "03_cross-species_transcriptome_analysis", "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(forcats)
  library(patchwork)
  library(ggbreak)
  library(DESeq2)
})
source(file.path(getwd(), "03_cross-species_transcriptome_analysis", "Functions.R"))

metadata  <- readRDS(file.path(in_dir, "metadata.rds"))
og_counts <- readRDS(file.path(in_dir, "og_counts.rds"))

species_palette <- c(
  "Anatherum_virginicum"     = "#1F77B4", "Andropogon_gerardii"      = "#FF7F0E",
  "Bothriochloa_laguroides"  = "#2CA02C", "Chrysopogon_zizanioides"  = "#D62728",
  "Elionorus_tripsacoides"   = "#9467BD", "Miscanthus_giganteus"     = "#8C564B",
  "Panicum_virgatum"         = "#98DF8A", "Schizachyrium_scoparium"  = "#E377C2",
  "Sorghastrum_nutans"       = "#7F7F7F", "Sorghum_bicolor"          = "#FF9896",
  "Themeda_triandra"         = "#E6AB02", "Tripsacum_dactyloides"    = "#17BECF",
  "Zea_diploperennis"        = "#BCBD22", "Zea_mays"                 = "#36eb11"
)
tissue_type_palette <- c("leaf" = "#30e004", "root" = "#ee6766", "stem" = "#07cde3")
season_palette      <- c("2023" = "#FC8D62", "2024" = "#66C2A5")  # RColorBrewer Set2

# ==========================================
# A) Libraries per species
# ==========================================
library_counts_species <- metadata %>%
  dplyr::count(species, name = "n") %>%
  mutate(species = fct_reorder(species, n))
species_order <- levels(library_counts_species$species)

species_bar_plot <- ggplot(library_counts_species, aes(x = n, y = species, fill = species)) +
  geom_col(width = 0.9) +
  geom_text(aes(label = n), hjust = -0.2, size = 2.8) +
  scale_fill_manual(values = species_palette) +
  scale_x_continuous(limits = c(0, 1300), expand = c(0, 0)) +
  scale_x_break(c(300, 1050), scales = 0.6) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 10) +
  labs(x = "number of libraries", y = NULL) +
  theme(legend.position = "none", axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = margin(5.5, 40, 5.5, 5.5))

ggsave(file.path(out_dir, "library_counts_species.tiff"), species_bar_plot,
       width = 4.6, height = 3.3, units = "in", dpi = 600, compression = "lzw")

# ==========================================
# B-C) Proportion of libraries per tissue and per season
# ==========================================
# Stacking: the first factor level is drawn at the right-hand end of each bar
proportion_bar_plot <- function(group_var, palette, levels) {
  df <- metadata %>%
    mutate(group = factor(as.character(.data[[group_var]]), levels = levels)) %>%
    dplyr::count(species, group, name = "n") %>%
    group_by(species) %>%
    mutate(percent = n / sum(n)) %>%
    ungroup() %>%
    mutate(species = factor(species, levels = species_order))
  ggplot(df, aes(y = species, x = percent, fill = group)) +
    geom_col(width = 0.9) +
    geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 2.8, color = "black") +
    scale_x_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
    scale_fill_manual(values = palette, breaks = sort(levels), name = NULL) +
    theme_classic(base_size = 10) +
    labs(x = "Proportion of Libraries", y = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          axis.text.x = element_text(size = 9), legend.position = "top",
          plot.margin = margin(5.5, 12, 5.5, 12))
}

ggsave(file.path(out_dir, "library_counts_tissue.tiff"), proportion_bar_plot("tissue", tissue_type_palette, c("leaf", "root", "stem")),
       width = 2.8, height = 3.3, units = "in", dpi = 600, compression = "lzw")
ggsave(file.path(out_dir, "library_counts_season.tiff"), proportion_bar_plot("season", season_palette, c("2024", "2023")),
       width = 2, height = 3.3, units = "in", dpi = 600, compression = "lzw")

# ==========================================
# D-E) Seasonal batch effect: leaf and root PCA colored by season
# ==========================================
season_pca_row <- function(tissue_name) {
  meta_t   <- metadata %>% filter(tissue == tissue_name)
  counts_t <- og_counts[, colnames(og_counts) %in% meta_t$sample_id, drop = FALSE]
  pca <- plot_pca_ggplot2(counts_t, meta_t, "season")
  ve  <- pca$variance_explained
  df  <- pca$pca_df %>% mutate(season = as.character(season))
  cat(tissue_name, ":", nrow(df), "libraries; PC1-4 % variance:", ve[1:4], "\n")
  panels <- lapply(2:4, function(k) {
    ggplot(df, aes(PC1, .data[[paste0("PC", k)]], color = season)) +
      geom_point(size = 0.6, alpha = 0.8) +
      scale_color_manual(values = season_palette) +
      labs(x = paste0("PC1 - ", ve[1], " % variance"), y = paste0("PC", k, " - ", ve[k], " % variance")) +
      theme_classic(base_size = 9) +
      theme(legend.position = "none")
  })
  wrap_plots(panels, nrow = 1) +
    plot_annotation(title = paste0(nrow(df), " ", tissue_name, " RNA-seq libraries | PCA analysis of batch effect across seasons"))
}

ggsave(file.path(out_dir, "pca_season_leaf.tiff"), season_pca_row("leaf"),
       width = 7, height = 2.4, units = "in", dpi = 600, compression = "lzw")
ggsave(file.path(out_dir, "pca_season_root.tiff"), season_pca_row("root"),
       width = 7, height = 2.4, units = "in", dpi = 600, compression = "lzw")

cat("Supplementary Figure 6 panels saved to", out_dir, "\n")
