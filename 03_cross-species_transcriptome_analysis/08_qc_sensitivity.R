# 08_qc_sensitivity.R
# Supplementary Figure 5: sensitivity of the global PCA to the library QC thresholds.
#   A  increasing mapping-rate thresholds (unfiltered, >35%, >50%, >75%)
#   B  increasing mapped-read thresholds (unfiltered, >100K, >200K, >500K)
# Rows are colored by tissue and by sampling week. Each PCA uses the 5,000 most variable
# orthogroups of the all-library VST (data/vst_all_libraries.rds, from 07_qc_all_libraries.R).
# Filtered panels start from the 2,685 libraries that passed QC (>35% mapping, >100K reads).
#
# Usage (from the repository root):
#   Rscript 03_cross-species_transcriptome_analysis/08_qc_sensitivity.R [in_dir] [out_dir]

args <- commandArgs(trailingOnly = TRUE)
in_dir  <- if (length(args) > 0) args[1] else file.path(getwd(), "03_cross-species_transcriptome_analysis", "data")
out_dir <- if (length(args) > 1) args[2] else file.path(getwd(), "03_cross-species_transcriptome_analysis", "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(matrixStats)
})

vst_mat    <- readRDS(file.path(in_dir, "vst_all_libraries.rds"))
meta_whole <- readRDS(file.path(in_dir, "metadata_all_libraries.rds"))
meta_filt  <- readRDS(file.path(in_dir, "metadata.rds"))

meta_whole <- meta_whole %>%
  mutate(
    species = gsub("_", " ", species),
    week_num = as.numeric(week)
  )

meta_filt <- meta_filt %>%
  mutate(
    mapping_rate = MappedReads / CleanReads,
    mapping_rate = ifelse(is.na(mapping_rate) | is.nan(mapping_rate), 0, mapping_rate),
    species = gsub("_", " ", species),
    week_num = as.numeric(week)
  )

tissue_palette <- c("leaf" = "#30e004", "root" = "#ee6766", "stem" = "#07cde3")

# Fast PCA function using top 5,000 variable features (standard DESeq2 plotPCA approach)
run_pca_subset_fast <- function(vst, metadata, subset_samples, ntop = 5000) {
  vst_sub <- vst[, subset_samples, drop = FALSE]
  
  # Select top ntop most variable features
  rv <- rowVars(vst_sub)
  select_idx <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  vst_top <- vst_sub[select_idx, , drop = FALSE]
  
  pca_res <- prcomp(t(vst_top))
  pca_df <- as.data.frame(pca_res$x)
  pca_df$sample_id <- rownames(pca_df)
  pca_df <- left_join(pca_df, metadata, by = "sample_id")
  var_exp <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 2)
  attr(pca_df, "variance_explained") <- var_exp
  return(pca_df)
}

samples_whole <- colnames(vst_mat)
samples_mr35  <- meta_filt$sample_id
samples_mr50  <- meta_filt %>% filter(mapping_rate >= 0.50) %>% pull(sample_id)
samples_mr75  <- meta_filt %>% filter(mapping_rate >= 0.75) %>% pull(sample_id)

samples_rc100k <- meta_filt$sample_id
samples_rc200k <- meta_filt %>% filter(MappedReads >= 200000) %>% pull(sample_id)
samples_rc500k <- meta_filt %>% filter(MappedReads >= 500000) %>% pull(sample_id)

cat("Libraries per threshold:", length(samples_whole), length(samples_mr35), length(samples_mr50),
    length(samples_mr75), length(samples_rc200k), length(samples_rc500k), "\n")
cat("Running fast PCAs on subsets...\n")
pca_whole  <- run_pca_subset_fast(vst_mat, meta_whole, samples_whole)
pca_mr35   <- run_pca_subset_fast(vst_mat, meta_filt, samples_mr35)
pca_mr50   <- run_pca_subset_fast(vst_mat, meta_filt, samples_mr50)
pca_mr75   <- run_pca_subset_fast(vst_mat, meta_filt, samples_mr75)

pca_rc100k <- run_pca_subset_fast(vst_mat, meta_filt, samples_rc100k)
pca_rc200k <- run_pca_subset_fast(vst_mat, meta_filt, samples_rc200k)
pca_rc500k <- run_pca_subset_fast(vst_mat, meta_filt, samples_rc500k)

plot_single_pca <- function(df, title_str, color_by = "tissue") {
  var_exp <- attr(df, "variance_explained")
  
  if (color_by == "tissue") {
    p <- ggplot(df, aes(x = PC1, y = PC2, color = tissue)) +
      geom_point(size = 0.08, alpha = 0.55) +
      scale_color_manual(values = tissue_palette, name = "Tissue")
  } else {
    p <- ggplot(df, aes(x = PC1, y = PC2, color = week_num)) +
      geom_point(size = 0.08, alpha = 0.55) +
      scale_color_viridis_c(option = "viridis", name = "Sampling Week")
  }
  
  p <- p +
    labs(
      title = title_str,
      x = paste0("PC1 (", var_exp[1], "%)"),
      y = paste0("PC2 (", var_exp[2], "%)")
    ) +
    theme_classic(base_size = 5.5) +
    theme(
      plot.title = element_text(size = 5.5, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 4.5, face = "bold"),
      axis.text = element_text(size = 3.5),
      axis.ticks = element_line(linewidth = 0.1),
      axis.line = element_line(linewidth = 0.1),
      legend.position = "none",
      plot.margin = margin(1.5, 1.5, 1.5, 1.5, "pt")
    )
  return(p)
}

# --- Row 1: Mapping Rate Sensitivity (Tissue) ---
p_mr1_t <- plot_single_pca(pca_whole, paste0("Unfiltered (", format(length(samples_whole), big.mark = ","), ")"), "tissue")
p_mr2_t <- plot_single_pca(pca_mr35,  paste0("Mapping Rate > 35% (", format(length(samples_mr35), big.mark = ","), ")"), "tissue")
p_mr3_t <- plot_single_pca(pca_mr50,  paste0("Mapping Rate > 50% (", format(length(samples_mr50), big.mark = ","), ")"), "tissue")
p_mr4_t <- plot_single_pca(pca_mr75,  paste0("Mapping Rate > 75% (", format(length(samples_mr75), big.mark = ","), ")"), "tissue")

# --- Row 2: Mapping Rate Sensitivity (Sampling Week) ---
p_mr1_w <- plot_single_pca(pca_whole, paste0("Unfiltered (", format(length(samples_whole), big.mark = ","), ")"), "week_num")
p_mr2_w <- plot_single_pca(pca_mr35,  paste0("Mapping Rate > 35% (", format(length(samples_mr35), big.mark = ","), ")"), "week_num")
p_mr3_w <- plot_single_pca(pca_mr50,  paste0("Mapping Rate > 50% (", format(length(samples_mr50), big.mark = ","), ")"), "week_num")
p_mr4_w <- plot_single_pca(pca_mr75,  paste0("Mapping Rate > 75% (", format(length(samples_mr75), big.mark = ","), ")"), "week_num")

# --- Row 3: Read Count Sensitivity (Tissue) ---
p_rc1_t <- plot_single_pca(pca_whole,  paste0("Unfiltered (", format(length(samples_whole), big.mark = ","), ")"), "tissue")
p_rc2_t <- plot_single_pca(pca_rc100k, paste0("Mapped Reads > 100K (", format(length(samples_rc100k), big.mark = ","), ")"), "tissue")
p_rc3_t <- plot_single_pca(pca_rc200k, paste0("Mapped Reads > 200K (", format(length(samples_rc200k), big.mark = ","), ")"), "tissue")
p_rc4_t <- plot_single_pca(pca_rc500k, paste0("Mapped Reads > 500K (", format(length(samples_rc500k), big.mark = ","), ")"), "tissue")

# --- Row 4: Read Count Sensitivity (Sampling Week) ---
p_rc1_w <- plot_single_pca(pca_whole,  paste0("Unfiltered (", format(length(samples_whole), big.mark = ","), ")"), "week_num")
p_rc2_w <- plot_single_pca(pca_rc100k, paste0("Mapped Reads > 100K (", format(length(samples_rc100k), big.mark = ","), ")"), "week_num")
p_rc3_w <- plot_single_pca(pca_rc200k, paste0("Mapped Reads > 200K (", format(length(samples_rc200k), big.mark = ","), ")"), "week_num")
p_rc4_w <- plot_single_pca(pca_rc500k, paste0("Mapped Reads > 500K (", format(length(samples_rc500k), big.mark = ","), ")"), "week_num")

# Combine into 4-row x 4-col 8x8 in layout
combined_4x4 <- (
  (p_mr1_t | p_mr2_t | p_mr3_t | p_mr4_t) /
  (p_mr1_w | p_mr2_w | p_mr3_w | p_mr4_w) /
  (p_rc1_t | p_rc2_t | p_rc3_t | p_rc4_t) /
  (p_rc1_w | p_rc2_w | p_rc3_w | p_rc4_w)
) + plot_annotation(
  title = "PCA Sensitivity Analysis to QC Thresholds (Colored by Tissue and Sampling Week)",
  theme = theme(plot.title = element_text(size = 7.5, face = "bold", hjust = 0.5))
)

out_png  <- file.path(out_dir, "pca_qc_sensitivity_combined_8x8.png")
out_tiff <- file.path(out_dir, "pca_qc_sensitivity_combined_8x8.tiff")

cat("Saving 4x4 grid figures (8x8 inches)...\n")
ggsave(out_png, combined_4x4, width = 8, height = 8, units = "in", dpi = 300)
ggsave(out_tiff, combined_4x4, width = 8, height = 8, units = "in", dpi = 600)

cat("Successfully generated 4x4 grid PCA sensitivity figure!\n")
