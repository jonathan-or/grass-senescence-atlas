# 04b_perform_GO_enrichment.R
# This script performs Gene Ontology (GO) enrichment analysis on orthogroups.
# It uses pre-processed GMT and annotation mapping objects derived from PANNZER outputs.

args <- commandArgs(trailingOnly = TRUE)
in_dir <- if (length(args) > 0) args[1] else file.path(getwd(), "03_cross-species_transcriptome_analysis", "data")
out_dir <- if (length(args) > 1) args[2] else file.path(getwd(), "03_cross-species_transcriptome_analysis", "plots")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load core functions
source(file.path(getwd(), "03_cross-species_transcriptome_analysis", "Functions.R"))

library(clusterProfiler)
library(ggplot2)
library(dplyr)

# 1. Load GO Annotation Dependencies (.rds)
cat("Loading GO annotation files...\n")
gmt_bp_df <- readRDS(file.path(in_dir, "og_go_df_BP.rds"))
gmt_mf_df <- readRDS(file.path(in_dir, "og_go_df_MF.rds"))
gmt_cc_df <- readRDS(file.path(in_dir, "og_go_df_CC.rds"))
gmt_data_df <- readRDS(file.path(in_dir, "gmt_data.rds"))

# 2. Example Analysis
# In a full run, this would be supplied by expression cluster definitions 
# or differential expression results. Here we generate a mock subset of orthogroups 
# to demonstrate the pipeline functions out-of-the-box.

cat("Demonstrating GO enrichment with a subset of orthogroups...\n")
og_counts <- readRDS(file.path(in_dir, "og_counts.rds"))

# Define background universe (all expressed orthogroups)
universe_ogs <- rownames(og_counts)

# Define a foreground set (mocked: top 100 expressed OGs by mean count across all samples)
og_means <- rowMeans(og_counts)
foreground_ogs <- names(sort(og_means, decreasing = TRUE)[1:100])

# 3. Perform Enrichment using specific universe
cat("Running enrichment against specific universe...\n")
enrichment_results <- perform_GO_enrichment_OG_specific_universe(
  og_vector = foreground_ogs,
  gmt_bp_df = gmt_bp_df,
  gmt_mf_df = gmt_mf_df,
  gmt_cc_df = gmt_cc_df,
  gmt_data_df = gmt_data_df,
  universe_ogs = universe_ogs,
  pAdjustMethod = "BH"
)

# 4. Perform and plot (Generic background)
cat("Generating visual plots...\n")
output_txt <- file.path(out_dir, "test_GO_enrichment_results.txt")
plot_results <- perform_and_plot_GO_OG(
  og_vector = foreground_ogs,
  gmt_bp_df = gmt_bp_df,
  gmt_mf_df = gmt_mf_df,
  gmt_cc_df = gmt_cc_df,
  gmt_data_df = gmt_data_df,
  output_txt = output_txt
)

# Export standard Lollipop Plot
if (!is.null(plot_results$lollipop_plot)) {
  tiff(file.path(out_dir, "GO_lollipop_plot.tiff"), width = 6, height = 4, units = "in", res = 600, compression = "lzw")
  print(plot_results$lollipop_plot)
  dev.off()
}

# Export Mini Bubble Plot
if (!is.null(plot_results$miniplot)) {
  tiff(file.path(out_dir, "GO_mini_bubble_plot.tiff"), width = 4, height = 4, units = "in", res = 600, compression = "lzw")
  print(plot_results$miniplot)
  dev.off()
}

cat("GO Enrichment completed successfully.\n")
