# 05_pseudotime_analysis.R
# This script performs pseudotime trajectory modeling on leaf tissue using PhenoPath.

args <- commandArgs(trailingOnly = TRUE)
in_dir <- if (length(args) > 0) args[1] else file.path(getwd(), "03_cross-species_transcriptome_analysis", "data")
out_dir <- if (length(args) > 1) args[2] else file.path(getwd(), "03_cross-species_transcriptome_analysis", "plots")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load core functions
source(file.path(getwd(), "03_cross-species_transcriptome_analysis", "Functions.R"))

library(dplyr)
library(ggplot2)
library(tidyr)
library(DESeq2)

# Set up specific custom palettes requested
species_palette <- c(
  "Anatherum_virginicum"     = "#1F77B4", "Andropogon_gerardii"      = "#FF7F0E",
  "Bothriochloa_laguroides"  = "#2CA02C", "Chrysopogon_zizanioides"  = "#D62728",
  "Elionorus_tripsacoides"   = "#9467BD", "Miscanthus_giganteus"     = "#8C564B",
  "Panicum_virgatum"         = "#98DF8A", "Schizachyrium_scoparium"  = "#E377C2",
  "Sorghastrum_nutans"       = "#7F7F7F", "Sorghum_bicolor"          = "#FF9896",
  "Themeda_triandra"         = "#E6AB02", "Tripsacum_dactyloides"    = "#17BECF",
  "Zea_diploperennis"        = "#BCBD22", "Zea_mays"                 = "#36eb11"
)

# 1. Load objects & filter for leaf libraries
cat("Loading consolidated objects and dropping non-leaf libraries...\n")
metadata <- readRDS(file.path(in_dir, "metadata.rds"))
og_counts <- readRDS(file.path(in_dir, "og_counts.rds"))

leaf_metadata <- metadata %>% filter(tissue == "leaf")
leaf_og_counts <- og_counts[, colnames(og_counts) %in% leaf_metadata$sample_id, drop = FALSE]

cat("Leaf libraries selected:", ncol(leaf_og_counts), "\n")

# 2. Initial PCAs
cat("Generating Initial Leaf PCAs...\n")
pca_leaf_species <- plot_pca_ggplot2(leaf_og_counts, leaf_metadata, "species", species_palette)
export_pca_plot(pca_leaf_species$pca_plot, file.path(out_dir, "SupFig2A_pca_leaf_species.tiff"), width = 2, height = 2, shape=19, point_size=0.5, point_alpha=0.8)

pca_leaf_julian <- plot_pca_ggplot2(leaf_og_counts, leaf_metadata, "julian_days", species_palette)
export_pca_plot(pca_leaf_julian$pca_plot, file.path(out_dir, "SupFig2B_pca_leaf_julian_days.tiff"), width = 2, height = 2, shape=19, point_size=0.5, point_alpha=0.8)

# 3. Trajectory Inference (PhenoPath)
cat("Checking for cached PhenoPath model...\n")
phenopath_cache <- file.path(in_dir, "leaf_pseudotime_model.rds")

if (file.exists(phenopath_cache)) {
  cat("Loading existing trajectory model...\n")
  leaf_pseudotime_result <- readRDS(phenopath_cache)
} else {
  cat("Running Bayesian trajectory inference (this may take a while)...\n")
  leaf_pseudotime_result <- calculate_pseudotime(leaf_og_counts, leaf_metadata, covariate = "species")
  saveRDS(leaf_pseudotime_result, phenopath_cache)
}

pseudotime_df <- leaf_pseudotime_result$pseudotime_df
fit_leaf <- leaf_pseudotime_result$fit

# Attach Pseudotime to metadata
leaf_metadata <- merge(leaf_metadata, pseudotime_df, by.x = "sample_id", by.y = "SampleID")

# 4. Styled Figure 2B PCA via inferred Pseudotime
pca_leaf_pseudotime <- plot_pca_ggplot2(leaf_og_counts, leaf_metadata, "Pseudotime")
p2b <- export_pca_plot(pca_leaf_pseudotime$pca_plot, file.path(out_dir, "Figure2B_pca_pseudotime.tiff"), width = 2.5, height = 2.5, point_size=1.2)
p2b_styled <- p2b + 
  geom_point(aes(fill = Pseudotime), shape = 21, color = "black", size = 1.2, stroke = 0.1) +
  scale_fill_viridis_c(option = "magma", direction = -1) + 
  theme_classic(base_size = 10) + theme(plot.title=element_blank(), axis.title=element_text(size=10), legend.position="right")

ggsave(file.path(out_dir, "Figure2B_pca_pseudotime_styled.tiff"), plot = p2b_styled, width = 2.5, height = 2.5, units = "in", dpi = 600, compression = "lzw")

# 5. Chronological Time vs Pseudotime Regression (SupFig 2D)
cat("Testing Correlation between Julian Days and Pseudotime...\n")
df <- leaf_metadata %>% filter(!is.na(julian_days), !is.na(Pseudotime))
cor_global <- cor.test(df$julian_days, df$Pseudotime, method = "pearson")
r_val <- round(unname(cor_global$estimate), 2)
p_val <- signif(cor_global$p.value, 2)

pseudotime_corr <- ggplot(df, aes(x = julian_days, y = Pseudotime)) +
  geom_point(shape = 21, fill = "deepskyblue3", color = "black", stroke = 0.3, size = 0.7, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, fill = "orange", linetype = "dashed", linewidth = 0.3, alpha = 0.3, color = "red") +
  labs(x = "time (julian days)", y = "pseudotime", title = "") +
  annotate("text", x = min(df$julian_days, na.rm=T), y = max(df$Pseudotime, na.rm=T), label = paste0("r = ", r_val, "\nP = ", p_val), hjust = 0, vjust = 1, size = 2.5) +
  scale_x_continuous(breaks = seq(200, 300, by = 25), limits = c(200, 300)) + 
  theme_classic(base_size = 9) + theme(legend.position = "none")

ggsave(file.path(out_dir, "SupFig2D_leaf_pseudotime_corr.tiff"), pseudotime_corr, width = 2, height = 2.2, units = "in", dpi = 600, compression = "lzw")

# 6. Marker Trajectories (Figure 2C)
cat("Plotting Marker Gene Trajectories...\n")
target_genes <- c("SAG12" = "OG0008818", "SGR1" = "OG0007133", "RBCS1" = "OG0004822", "PEPC1" = "OG0017925")
gene_colors <- c("SAG12" = "purple", "SGR1" = "#d95f02", "RBCS1" = "forestgreen", "PEPC1" = "chartreuse3")

dds <- DESeqDataSetFromMatrix(countData = leaf_og_counts, colData = leaf_metadata, design = ~1)
dds <- estimateSizeFactors(dds)
vst_mat <- assay(vst(dds, blind = TRUE))

expr_vals <- vst_mat[target_genes, ]
scaled_expr <- t(apply(expr_vals, 1, scale))
colnames(scaled_expr) <- colnames(expr_vals)

expr_df <- as.data.frame(t(scaled_expr))
colnames(expr_df) <- names(target_genes)
expr_df$sample_id <- rownames(expr_df)

plot_df <- expr_df %>%
  pivot_longer(cols = all_of(names(target_genes)), names_to = "Gene", values_to = "Expression") %>%
  left_join(leaf_metadata %>% select(sample_id, Pseudotime), by = "sample_id")

pt_markers <- ggplot(plot_df, aes(x = Pseudotime, y = Expression, color = Gene, fill = Gene)) +
  geom_point(alpha = 0.6, size = 0.8, stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x, linewidth = 1, se = TRUE, alpha=0.7) +
  scale_color_manual(values = gene_colors) + scale_fill_manual(values = gene_colors) +
  theme_classic(base_size = 10) +
  labs(x = "pseudotime", y = "expression (z-score)", color = "orthogroup", fill = "orthogroup") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(file.path(out_dir, "Figure2C_pseudotime_markers.tiff"), pt_markers, width = 2.7, height = 2.5, units = "in", dpi = 600, compression = "lzw")

# 7. GO Enrichment Pipeline for Pseudotime Extremes
cat("Performing Pseudotime Modulation GO enrichment...\n")
pseudotime_effects_df <- data.frame(Gene = fit_leaf$feature_names, PseudotimeEffect = fit_leaf$m_lambda)
top_threshold <- quantile(pseudotime_effects_df$PseudotimeEffect, 0.99)
bottom_threshold <- quantile(pseudotime_effects_df$PseudotimeEffect, 0.01)

# Plot histogram of modulations
pt_mod_leaf <- ggplot(pseudotime_effects_df, aes(x = PseudotimeEffect)) +
  geom_histogram(binwidth = 0.2, fill = "deepskyblue3", alpha = 0.7, color = "black", size=0.2) +
  geom_vline(xintercept = c(top_threshold, bottom_threshold), linetype = "dashed", color = "red") +
  labs(x = "expression modulation (a.u.)", y = "orthogroup frequency") +
  annotate("text", x = top_threshold + 0.5, y = 500, label = "1%", color = "red", angle = 90) +
  annotate("text", x = bottom_threshold - 0.5, y = 500, label = "1%", color = "red", angle = 90) +
  theme_classic(base_size = 10)

ggsave(file.path(out_dir, "SupFig2E_pseudotime_modulation_genes.tiff"), pt_mod_leaf, width = 2.5, height = 2, units = "in", dpi = 600, compression = "lzw")

# Extract the 1% most positive and negative OGs
top_ogs_vector <- pseudotime_effects_df$Gene[pseudotime_effects_df$PseudotimeEffect >= top_threshold]
bottom_ogs_vector <- pseudotime_effects_df$Gene[pseudotime_effects_df$PseudotimeEffect <= bottom_threshold]

cat("Evaluating Negative Modulation (Photosynthesis shift)...\n")
gmt_bp_df <- readRDS(file.path(in_dir, "og_go_df_BP.rds"))
gmt_mf_df <- readRDS(file.path(in_dir, "og_go_df_MF.rds"))
gmt_cc_df <- readRDS(file.path(in_dir, "og_go_df_CC.rds"))
gmt_data_df <- readRDS(file.path(in_dir, "gmt_data.rds"))

go_neg_res <- perform_and_plot_GO_OG(bottom_ogs_vector, gmt_bp_df, gmt_mf_df, gmt_cc_df, gmt_data_df, file.path(out_dir, "leaf_neg_OG_GO_enrichment.txt"))
if (!is.null(go_neg_res$lollipop_plot)) {
  lp_plot <- go_neg_res$lollipop_plot + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme_classic(base_size = 14)
  ggsave(file.path(out_dir, "SupFig2F_GO_neg_modulation_OG.tiff"), lp_plot, width = 7.5, height = 3.6, units = "in", dpi = 600, compression = "lzw")
  ggsave(file.path(out_dir, "SupFig2G_GO_neg_modulation_small_OG.tiff"), go_neg_res$miniplot, width = 3, height = 1, units = "in", dpi = 600, compression = "lzw")
}

cat("Evaluating Positive Modulation (Senescence shift)...\n")
go_pos_res <- perform_and_plot_GO_OG(top_ogs_vector, gmt_bp_df, gmt_mf_df, gmt_cc_df, gmt_data_df, file.path(out_dir, "leaf_pos_OG_GO_enrichment.txt"))
if (!is.null(go_pos_res$lollipop_plot)) {
  lp_pos_plot <- go_pos_res$lollipop_plot + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme_classic(base_size = 14)
  ggsave(file.path(out_dir, "GO_pos_modulation_OG.tiff"), lp_pos_plot, width = 7.5, height = 3.6, units = "in", dpi = 600, compression = "lzw")
}
if (!is.null(go_pos_res$miniplot)) {
  ggsave(file.path(out_dir, "GO_pos_modulation_small_OG.tiff"), go_pos_res$miniplot, width = 1, height = 3, units = "in", dpi = 600, compression = "lzw")
}

cat("Pseudotime and Trajectory Assessment Modeling Complete.\n")

cat("Identifying the core photosynthetic orthogroup set for downstream index calculations...\n")
photo_go_terms <- c(
  "GO:0009535", "GO:0015979", "GO:0009523", "GO:0009522", "GO:0009765",
  "GO:0055035", "GO:0009941", "GO:0010287", "GO:0016168",
  "GO:0009416", "GO:0009654", "GO:0018298", "GO:0009532", "GO:0019904",
  "GO:0019253", "GO:0022900", "GO:0010207", "GO:0009773", "GO:0010598",
  "GO:0015995", "GO:0009570", "GO:0031977", "GO:0006094", "GO:0030388",
  "GO:0010206", "GO:0042170", "GO:0045156", "GO:0010319",
  "GO:0009579", "GO:0009644", "GO:0010196", "GO:0016620", "GO:0009534",
  "GO:0050661", "GO:0031969", "GO:0009658", "GO:0009642", "GO:0009543",
  "GO:0051537", "GO:0051536", "GO:0042651", "GO:0009055","GO:0015977"
)
neg_go_res <- read.table(file.path(out_dir, "leaf_neg_OG_GO_enrichment.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
photo_go_results <- neg_go_res %>%
  filter(ID %in% photo_go_terms) %>%
  pull(geneID)

photo_og_vector <- unique(unlist(strsplit(photo_go_results, "/")))
saveRDS(photo_og_vector, file.path(in_dir, "photo_ogs.rds"))

# Map to Zm ids
gene_og_map <- readRDS(file.path(in_dir, "gene_og_mapping.rds"))
zm_photo_map <- gene_og_map %>% 
  filter(Orthogroup %in% photo_og_vector & grepl("^Zm", GeneID))
write.table(zm_photo_map, file.path(in_dir, "photo_ogs_Zm.txt"), sep="\t", quote=FALSE, row.names=FALSE)
cat("Exported photo_ogs.rds and photo_ogs_Zm.txt to data directory.\n")
