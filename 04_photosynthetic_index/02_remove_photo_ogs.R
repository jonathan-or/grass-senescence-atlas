# 04_photosynthetic_index/02_remove_photo_ogs.R
# Removes photosynthesis orthogroups to prevent downstream trait circularity during network generation

args <- commandArgs(trailingOnly = TRUE)
in_dir <- if (length(args) > 0) args[1] else file.path(getwd(), "03_cross-species_transcriptome_analysis", "data")
out_dir <- if (length(args) > 1) args[2] else file.path(getwd(), "04_photosynthetic_index", "data")
plot_dir <- if (length(args) > 2) args[3] else file.path(getwd(), "04_photosynthetic_index", "plots")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

source(file.path(getwd(), "03_cross-species_transcriptome_analysis", "Functions.R"))

library(dplyr)
library(tidyr)
library(ggplot2)

cat("Loading prerequisite mapping objects...\n")
metadata <- readRDS(file.path(in_dir, "metadata.rds"))
og_counts <- readRDS(file.path(in_dir, "og_counts.rds"))
leaf_metadata_ps_index <- readRDS(file.path(in_dir, "leaf_metadata_ps_index.rds"))
photo_ogs <- readRDS(file.path(in_dir, "photo_ogs.rds"))

# Step 1: Establish base expressions per tissue
cat("Isolating specific tissue expression domains...\n")
leaf_og_counts <- og_counts[, colnames(og_counts) %in% leaf_metadata_ps_index$sample_id]

root_metadata <- metadata %>% filter(tissue == "root") %>% mutate(week = as.integer(week))
root_og_counts <- og_counts[, colnames(og_counts) %in% root_metadata$sample_id]

stem_metadata <- metadata %>% filter(tissue == "stem") %>% mutate(week = as.integer(week))
stem_og_counts <- og_counts[, colnames(og_counts) %in% stem_metadata$sample_id]

# Step 2: Calculate weekly average ps_index projected chronologically
weekly_ps_index <- leaf_metadata_ps_index %>%
  group_by(week, accession) %>%
  summarize(avg_ps_index = mean(ps_index, na.rm = TRUE), .groups = "drop") %>%
  mutate(week = as.integer(week))

cat("Merging and projecting ps_index chronologically across root and stem structures...\n")
root_metadata_ps_index <- root_metadata %>%
  left_join(weekly_ps_index, by = c("week", "accession")) %>%
  dplyr::rename(ps_index = avg_ps_index) %>%
  drop_na(ps_index)

stem_metadata_ps_index <- stem_metadata %>%
  left_join(weekly_ps_index, by = c("week", "accession")) %>%
  dplyr::rename(ps_index = avg_ps_index) %>%
  drop_na(ps_index)

# Step 3: Eliminate Circular nodes (115 Photo OGs)
cat("Stripping photo_ogs systematically from mapping profiles...\n")

leaf_og_counts_dps <- leaf_og_counts[!rownames(leaf_og_counts) %in% photo_ogs, ]
root_og_counts_dps <- root_og_counts[!rownames(root_og_counts) %in% photo_ogs, ]
stem_og_counts_dps <- stem_og_counts[!rownames(stem_og_counts) %in% photo_ogs, ]

# Validate column mapping safety
root_og_counts_dps <- root_og_counts_dps[, root_metadata_ps_index$sample_id, drop = FALSE]
stem_og_counts_dps <- stem_og_counts_dps[, stem_metadata_ps_index$sample_id, drop = FALSE]

cat("Dimensions of stripped leaf matrices:", dim(leaf_og_counts_dps), "\n")
cat("Dimensions of stripped root matrices:", dim(root_og_counts_dps), "\n")
cat("Dimensions of stripped stem matrices:", dim(stem_og_counts_dps), "\n")

# Step 4: Export visual evaluations (PCA) verifying matrix stability
cat("Constructing post-drop PCA evaluations...\n")
pca_leaf <- plot_pca_ggplot2(leaf_og_counts_dps, leaf_metadata_ps_index, "ps_index", species_palette, accession_palette, 1, 2)
pca_root <- plot_pca_ggplot2(root_og_counts_dps, root_metadata_ps_index, "ps_index", species_palette, accession_palette, 1, 2)
pca_stem <- plot_pca_ggplot2(stem_og_counts_dps, stem_metadata_ps_index, "ps_index", species_palette, accession_palette, 1, 2)

ggsave(file.path(plot_dir, "pca_leaf_dps.tiff"), pca_leaf$pca_plot, width=5, height=4, dpi=600, compression="lzw")
ggsave(file.path(plot_dir, "pca_root_dps.tiff"), pca_root$pca_plot, width=5, height=4, dpi=600, compression="lzw")
ggsave(file.path(plot_dir, "pca_stem_dps.tiff"), pca_stem$pca_plot, width=5, height=4, dpi=600, compression="lzw")

# Step 5: Master Array Export
cat("Exporting 6 definitive structural arrays locally for WGCNA networks...\n")
saveRDS(leaf_og_counts_dps, file.path(out_dir, "leaf_og_counts_dps.rds"))
saveRDS(leaf_metadata_ps_index, file.path(out_dir, "leaf_metadata_ps_index.rds"))

saveRDS(root_og_counts_dps, file.path(out_dir, "root_og_counts_dps.rds"))
saveRDS(root_metadata_ps_index, file.path(out_dir, "root_metadata_ps_index.rds"))

saveRDS(stem_og_counts_dps, file.path(out_dir, "stem_og_counts_dps.rds"))
saveRDS(stem_metadata_ps_index, file.path(out_dir, "stem_metadata_ps_index.rds"))

cat("Success! 06_gene_network_analysis mapped and ready.\n")
