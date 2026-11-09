# 05_maize_gradient_source_sink.R
# Supplementary Figure 11: expression of the leaf source and sink gene sets along the
# base-to-tip maize leaf gradient (15 sections; Wang et al. 2014 samples in Zea_mays_counts.txt).
#   A  % of the transcriptome from maize members of the source and sink orthogroups
#   B  as A, plus the ps_index genes and sink + ps_index genes
#   C  continuous DESeq2 model (~ leaf section): up / down / not significant / unexpressed
#      counts per gene set (table)
#   D  z-scores of the source genes significantly downregulated towards the tip
#
# Gene sets (maize members, from gene_og_mapping.rds) of the leaf co-expression modules
# (06_gene_network_analysis/results/06a_gene_module_assignments_leaf.csv):
#   source = m2_blue + m7_black + m3_brown (1,040 OGs); sink = m10_purple + m4_yellow + m6_red
# ps_index genes: 05_N_quant/data/photo_genes_Zm.txt
#
# Usage (from the repository root):
#   Rscript 05_N_quant/05_maize_gradient_source_sink.R [plot_dir]

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
in_dir   <- file.path(getwd(), "05_N_quant", "data")
plot_dir <- if (length(args) > 0) args[1] else file.path(getwd(), "05_N_quant", "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

# ==========================================
# 1. Inputs
# ==========================================
metadata <- read.csv(file.path(in_dir, "metadata.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
counts   <- read.delim(file.path(in_dir, "Zea_mays_counts.txt"), row.names = 1, stringsAsFactors = FALSE)

# Collapse transcripts to genes
rownames(counts) <- sub("_T0\\d{2}$", "", rownames(counts))
counts <- aggregate(counts, by = list(Gene = rownames(counts)), FUN = sum)
rownames(counts) <- counts$Gene
counts$Gene <- NULL

# Wang et al. leaf sections 1-15
wang_meta   <- metadata %>% filter(batch == "wang", dev_stage != 16)
wang_counts <- counts[, colnames(counts) %in% wang_meta$sample_id]
wang_meta   <- wang_meta[match(colnames(wang_counts), wang_meta$sample_id), ]
wang_meta$dev_stage     <- as.numeric(wang_meta$dev_stage)
wang_meta$dev_stage_num <- wang_meta$dev_stage

# Source / sink gene sets from the leaf network
modules  <- read_csv(file.path(getwd(), "06_gene_network_analysis", "results", "06a_gene_module_assignments_leaf.csv"),
                     show_col_types = FALSE)
zm_map   <- readRDS(file.path(getwd(), "03_cross-species_transcriptome_analysis", "data", "gene_og_mapping.rds")) %>%
  filter(Species == "Zea_mays")
maize_genes_of <- function(colors) unique(zm_map$GeneID[zm_map$Orthogroup %in% modules$OG[modules$module_color %in% colors]])
source_genes <- maize_genes_of(c("blue", "black", "brown"))
sink_genes   <- maize_genes_of(c("purple", "yellow", "red"))
photo_genes  <- read.table(file.path(in_dir, "photo_genes_Zm.txt"), header = FALSE, stringsAsFactors = FALSE)$V1
cat("Gene sets: source", length(source_genes), "| sink", length(sink_genes), "| ps_index", length(photo_genes), "\n")

# ==========================================
# 2. Panels A-B: % of the transcriptome (size-factor normalized counts)
# ==========================================
dds_sf <- DESeqDataSetFromMatrix(countData = wang_counts, colData = wang_meta, design = ~ 1)
dds_sf <- dds_sf[rowSums(counts(dds_sf)) > 10, ]
dds_sf <- estimateSizeFactors(dds_sf)
norm_counts <- counts(dds_sf, normalized = TRUE)

sets <- list(
  "Source Genes"    = intersect(source_genes, rownames(norm_counts)),
  "Sink Genes"      = intersect(sink_genes, rownames(norm_counts)),
  "PS_Index Genes"  = intersect(photo_genes, rownames(norm_counts))
)
sets[["Sink + PS_Index"]] <- unique(c(sets[["Sink Genes"]], sets[["PS_Index Genes"]]))
cat("Expressed genes:", paste(names(sets), lengths(sets), sep = " = ", collapse = "; "), "\n")

df_pct_long <- bind_rows(lapply(names(sets), function(s) {
  tibble(sample_id = wang_meta$sample_id, dev_stage = wang_meta$dev_stage, Gene_Set = s,
         pct_count = colSums(norm_counts[sets[[s]], wang_meta$sample_id]) / colSums(norm_counts[, wang_meta$sample_id]) * 100)
}))
summary_pct <- df_pct_long %>%
  group_by(dev_stage, Gene_Set) %>%
  summarise(n = n(), mean_pct = mean(pct_count), se_pct = sd(pct_count) / sqrt(n),
            ymin = mean_pct - se_pct, ymax = mean_pct + se_pct, .groups = "drop")

set_colors <- c("Sink + PS_Index" = "#FF7F0E", "PS_Index Genes" = "#2CA02C",
                "Sink Genes" = "#9467BD", "Source Genes" = "#1F77B4")

plot_pct <- function(keep, y_max, y_step) {
  pts <- df_pct_long %>% filter(Gene_Set %in% keep) %>% mutate(Gene_Set = factor(Gene_Set, levels = keep))
  smr <- summary_pct %>% filter(Gene_Set %in% keep) %>% mutate(Gene_Set = factor(Gene_Set, levels = keep))
  ggplot() +
    geom_point(data = pts, aes(dev_stage, pct_count, color = Gene_Set),
               alpha = 0.35, size = 0.6, position = position_jitter(width = 0.15, height = 0)) +
    geom_ribbon(data = smr, aes(dev_stage, ymin = ymin, ymax = ymax, fill = Gene_Set), alpha = 0.2, color = NA) +
    geom_line(data = smr, aes(dev_stage, mean_pct, color = Gene_Set), linewidth = 0.65) +
    geom_point(data = smr, aes(dev_stage, mean_pct, color = Gene_Set), size = 1.1) +
    scale_color_manual(values = set_colors, name = "Gene Set") +
    scale_fill_manual(values = set_colors, name = "Gene Set") +
    scale_x_continuous(breaks = 1:15, limits = c(0.7, 15.3)) +
    scale_y_continuous(limits = c(0, y_max), breaks = seq(0, y_max, by = y_step)) +
    labs(x = "leaf section", y = "percentage (%) of total transcriptome") +
    theme_classic(base_size = 8) +
    theme(axis.title = element_text(size = 7.5), axis.text.x = element_text(size = 5, color = "black"),
          axis.text.y = element_text(size = 6.5, color = "black"), plot.margin = margin(4, 4, 4, 4, "pt"),
          legend.title = element_text(size = 7), legend.text = element_text(size = 6),
          legend.key.size = unit(9, "pt"))
}

ggsave(file.path(plot_dir, "wang_source_sink_percent_expressed_legend.png"),
       plot_pct(c("Source Genes", "Sink Genes"), 12, 3), width = 3.2, height = 2.5, units = "in", dpi = 300)
ggsave(file.path(plot_dir, "wang_all4_percent_expressed_legend.png"),
       plot_pct(c("Sink + PS_Index", "PS_Index Genes", "Sink Genes", "Source Genes"), 50, 10),
       width = 3.4, height = 2.5, units = "in", dpi = 300)

# ==========================================
# 3. Panel C: continuous gradient model (~ leaf section)
# ==========================================
dds <- DESeqDataSetFromMatrix(countData = wang_counts, colData = wang_meta, design = ~ dev_stage_num)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)
res_grad <- as.data.frame(results(dds, name = "dev_stage_num")) %>%
  rownames_to_column("geneid") %>%
  filter(!is.na(padj))

classify <- function(genes, label) {
  tested <- intersect(genes, res_grad$geneid)
  r <- res_grad %>% filter(geneid %in% tested)
  n <- c(up   = sum(r$log2FoldChange > 0 & r$padj < 0.05),
         down = sum(r$log2FoldChange < 0 & r$padj < 0.05),
         ns   = sum(r$padj >= 0.05),
         unexpressed = length(setdiff(genes, tested)))
  tibble(Gene_Set = label,
         Category = c("Upregulated towards tip (LFC > 0, padj < 0.05)",
                      "Downregulated towards tip (LFC < 0, padj < 0.05)",
                      "Not significant (padj >= 0.05)",
                      "Not expressed / zero-count in dataset"),
         Gene_Count = unname(n),
         Pct_of_Set = c(round(100 * n[1:3] / length(tested), 1), 0))
}
de_table <- bind_rows(classify(source_genes, "Source Genes"), classify(sink_genes, "Sink Genes"))
print(as.data.frame(de_table))
write_csv(de_table, file.path(plot_dir, "de_source_sink_intersections_continuous.csv"))

# ==========================================
# 4. Panel D: z-scores of source genes downregulated towards the tip
# ==========================================
source_tested <- intersect(source_genes, res_grad$geneid)
source_down   <- res_grad %>% filter(geneid %in% source_tested, log2FoldChange < 0, padj < 0.05) %>% pull(geneid)
title_d <- paste0(length(source_down), " genes | ", round(100 * length(source_down) / length(source_tested)), "% of source set")

vst_z <- t(scale(t(assay(vst(dds, blind = FALSE)))))
df_long <- as.data.frame(vst_z[source_down, ]) %>%
  rownames_to_column("geneid") %>%
  pivot_longer(-geneid, names_to = "sample_id", values_to = "zscore") %>%
  left_join(wang_meta %>% select(sample_id, dev_stage), by = "sample_id")
df_mean <- df_long %>% group_by(sample_id, dev_stage) %>% summarise(zscore = mean(zscore), .groups = "drop")

p_d <- ggplot() +
  geom_jitter(data = df_long, aes(dev_stage, zscore), color = "#1F77B4", alpha = 0.03, size = 0.3, width = 0.07) +
  geom_boxplot(data = df_long, aes(dev_stage, zscore, group = factor(dev_stage)), outlier.shape = NA,
               fill = "white", color = "#1F77B4", alpha = 0.65, width = 0.35, linewidth = 0.4) +
  geom_smooth(data = df_mean, aes(dev_stage, zscore), method = "loess", formula = y ~ x, color = "#0C4C8A",
              fill = "#1F77B4", se = TRUE, linewidth = 0.9, alpha = 0.22) +
  scale_x_continuous(breaks = 1:15, limits = c(0.7, 15.3)) +
  scale_y_continuous(limits = c(-3.5, 3.5)) +
  labs(title = title_d, x = "leaf section", y = "gene expression (z-score)") +
  theme_classic(base_size = 8) +
  theme(plot.title = element_text(size = 8, hjust = 0.5), axis.title = element_text(size = 7.5),
        axis.text.x = element_text(size = 5, color = "black"), axis.text.y = element_text(size = 6.5, color = "black"),
        plot.margin = margin(4, 4, 4, 4, "pt"))
ggsave(file.path(plot_dir, "wang_source_subset_zscore_all_sig_2.5in.png"), p_d, width = 2.5, height = 2.5, units = "in", dpi = 300)

cat("Supplementary Figure 11 panels saved to", plot_dir, "\n")
