############################################################
## Reference underground network (root + stem) — WGCNA
## - Based on end_of_season/08_wgcna_undg_MGC.R
## - Uses pre-computed objects for reproducibility
############################################################

suppressPackageStartupMessages({
  library(WGCNA)
  library(tidyverse)
  library(matrixStats)
  library(ggrepel)
  library(uwot)
  library(patchwork)
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(scales)
  library(stringr)
  library(plotly)
  library(DESeq2)
})

# --- Reproducibility & WGCNA threading ---
set.seed(42)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

source(file.path(getwd(), "03_cross-species_transcriptome_analysis", "Functions.R"))

# =========================
# 0) Paths
# =========================
ps_dir   <- file.path(getwd(), "04_photosynthetic_index", "data")
data_dir <- file.path(getwd(), "06_gene_network_analysis", "data")
annot_dir <- file.path(getwd(), "03_cross-species_transcriptome_analysis", "data")
out_dir  <- file.path(getwd(), "06_gene_network_analysis", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =========================
# 1) Load inputs
# =========================
cat("Loading input objects...\n")
root_metadata_ps_index <- readRDS(file.path(ps_dir, "root_metadata_ps_index.rds"))
stem_metadata_ps_index <- readRDS(file.path(ps_dir, "stem_metadata_ps_index.rds"))
root_og_counts_dps     <- readRDS(file.path(ps_dir, "root_og_counts_dps.rds"))
stem_og_counts_dps     <- readRDS(file.path(ps_dir, "stem_og_counts_dps.rds"))

gene_og_mapping        <- readRDS(file.path(annot_dir, "gene_og_mapping.rds"))
og_go_df_BP            <- readRDS(file.path(annot_dir, "og_go_df_BP.rds"))
og_go_df_MF            <- readRDS(file.path(annot_dir, "og_go_df_MF.rds"))
og_go_df_CC            <- readRDS(file.path(annot_dir, "og_go_df_CC.rds"))
gmt_data               <- readRDS(file.path(annot_dir, "gmt_data.rds"))

# =========================
# 2) Combine underground metadata & counts
# =========================
root_meta <- root_metadata_ps_index %>%
  mutate(sample_id = as.character(sample_id), tissue = "root")
stem_meta <- stem_metadata_ps_index %>%
  mutate(sample_id = as.character(sample_id), tissue = "stem")

underground_metadata_ps_index <- bind_rows(root_meta, stem_meta)

# Outlier removal
outlier_ids <- c("AC24BRTL002869")
underground_metadata_ps_index <- underground_metadata_ps_index %>%
  dplyr::filter(!sample_id %in% outlier_ids)

# Combine counts
root_counts <- as.matrix(root_og_counts_dps)
stem_counts <- as.matrix(stem_og_counts_dps)
common_genes <- intersect(rownames(root_counts), rownames(stem_counts))

og_counts_underground <- cbind(root_counts[common_genes, ], stem_counts[common_genes, ])

# Keep only samples in metadata
keep_ids <- intersect(colnames(og_counts_underground), underground_metadata_ps_index$sample_id)
og_counts_underground <- og_counts_underground[, keep_ids, drop = FALSE]
underground_metadata_ps_index <- underground_metadata_ps_index %>%
  filter(sample_id %in% keep_ids) %>%
  dplyr::slice(match(colnames(og_counts_underground), sample_id))

stopifnot(identical(colnames(og_counts_underground), underground_metadata_ps_index$sample_id))
message("Underground cohort: ", ncol(og_counts_underground), " samples; ", nrow(og_counts_underground), " OGs.")

# =========================
# 3) PCA & VST (Global)
# =========================
pca_underground_ps_index <- plot_pca_ggplot2(
  metadata = underground_metadata_ps_index,
  count_data = og_counts_underground,
  color_by = "tissue",
  species_palette = species_palette,
  accession_palette = accession_palette,
  x_pc = 1, y_pc = 2
)

vst_matrix_underground <- as.matrix(pca_underground_ps_index$vst_matrix)
mode(vst_matrix_underground) <- "numeric"

# =========================
# 4) Define remobilizer cohort
# =========================
remobilizer_species    <- c("Sorghastrum_nutans")
remobilizer_accessions <- c("Td_Pete", "Td_KS_C9", "Ag_Niagara")

meta_underground_remobilizer <- underground_metadata_ps_index %>%
  filter(species %in% remobilizer_species | accession %in% remobilizer_accessions)

common_ids_rem <- intersect(colnames(vst_matrix_underground), meta_underground_remobilizer$sample_id)
expr_underground_rem <- vst_matrix_underground[, common_ids_rem, drop = FALSE]
meta_underground_remobilizer <- meta_underground_remobilizer %>%
  filter(sample_id %in% common_ids_rem) %>%
  dplyr::slice(match(colnames(expr_underground_rem), sample_id))

# =========================
# 5) Build or load network
# =========================
net_und_path <- file.path(data_dir, "net_und.rds")
DatExpr_und_path <- file.path(data_dir, "DatExpr_und.rds")

if (file.exists(net_und_path) && file.exists(DatExpr_und_path)) {
  cat("Loading pre-computed WGCNA network objects...\n")
  net_und <- readRDS(net_und_path)
  DatExpr_und <- readRDS(DatExpr_und_path)
  softPower_undg <- net_und$power
} else {
  cat("Building WGCNA network from scratch...\n")
  
  # Subset to top 3,000 most variable OGs within cohort
  og_var <- matrixStats::rowVars(expr_underground_rem)
  ord    <- order(og_var, decreasing = TRUE)
  K      <- min(3000L, length(ord))
  sel    <- ord[seq_len(K)]
  expr_und_rem_3k <- expr_underground_rem[sel, , drop = FALSE]
  
  DatExpr_und <- t(expr_und_rem_3k)
  
  # Drop zero-MAD genes
  mad_vals <- matrixStats::colMads(DatExpr_und)
  DatExpr_und <- DatExpr_und[, mad_vals > 0, drop = FALSE]
  
  gsg <- WGCNA::goodSamplesGenes(DatExpr_und, verbose = 3)
  if (!gsg$allOK) {
    DatExpr_und <- DatExpr_und[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
    meta_underground_remobilizer <- dplyr::filter(meta_underground_remobilizer, sample_id %in% rownames(DatExpr_und))
  }
  meta_underground_remobilizer <- meta_underground_remobilizer[match(rownames(DatExpr_und), meta_underground_remobilizer$sample_id), ]
  
  # Soft-threshold selection
  powers <- c(1:14, seq(16, 24, 2))
  sft <- pickSoftThreshold(DatExpr_und, dataIsExpr = TRUE, networkType = "signed",
                           corFnc = "bicor", corOptions = list(maxPOutliers = 0.1),
                           powerVector = powers, verbose = 5)
  
  r2 <- -sign(sft$fitIndices[,3]) * sft$fitIndices[,2]
  cand <- which(r2 >= 0.6)
  softPower_undg <- if (length(cand) > 0) sft$fitIndices[min(cand), 1] else sft$fitIndices[which.max(r2), 1]
  message("Selected soft-threshold β (underground): ", softPower_undg)
  
  net_und <- blockwiseModules(
    DatExpr_und,
    power             = softPower_undg,
    networkType       = "signed",
    corType           = "bicor",
    maxPOutliers      = 0.1,
    TOMType           = "signed",
    deepSplit         = 2,
    minModuleSize     = 50,
    mergeCutHeight    = 0.40,
    pamRespectsDendro = TRUE,
    minKMEtoStay      = 0.3,
    reassignThreshold = 0,
    numericLabels     = TRUE,
    saveTOMs          = FALSE,
    verbose           = 3
  )
  
  saveRDS(net_und, net_und_path)
  saveRDS(DatExpr_und, DatExpr_und_path)
}

# Ensure metadata matches DatExpr
meta_underground_remobilizer <- meta_underground_remobilizer %>%
  filter(sample_id %in% rownames(DatExpr_und)) %>%
  dplyr::slice(match(rownames(DatExpr_und), sample_id))

moduleColors_und     <- labels2colors(net_und$colors)
moduleEigengenes_und <- orderMEs(net_und$MEs)

# Rename to uME_ prefix for underground
me_nums_all   <- as.integer(sub("^ME_?", "", colnames(moduleEigengenes_und)))
me_colors_all <- WGCNA::labels2colors(me_nums_all)
colnames(moduleEigengenes_und) <- paste0("uME_", me_colors_all)

# =========================================
# 6) Module summary & Sizes
# =========================================
modules_und <- table(moduleColors_und) %>%
  sort(decreasing = TRUE) %>%
  as.data.frame() %>%
  dplyr::rename(color = moduleColors_und, size = Freq) %>%
  mutate(
    module_num = ifelse(color %in% c("gray", "grey"), 0,
                        as.numeric(factor(color, levels = setdiff(unique(color), c("gray", "grey"))))),
    module_label = paste0("um", module_num, "_", color)
  ) %>%
  arrange(module_num)

pdf(file.path(out_dir, "03_module_sizes_und.pdf"), width = 2, height = 2.5)
ggplot(modules_und, aes(x = size, y = reorder(module_label, size), fill = color)) +
  geom_col(color = "black", linewidth = 0.2) +
  scale_fill_identity() +
  scale_x_continuous(position = "top") +
  labs(x = "size (OGs)", y = "co-expression module") +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_text(size = 9), plot.margin = margin(5, 5, 5, 5))
dev.off()

# ==========================================
# 7) Correlations with ps_index & tissue
# ==========================================
keep_ps    <- !is.na(meta_underground_remobilizer$ps_index)
MEs_ps_und <- moduleEigengenes_und[keep_ps, , drop = FALSE]
ps_und     <- as.numeric(meta_underground_remobilizer$ps_index[keep_ps])
tissue_num <- ifelse(meta_underground_remobilizer$tissue[keep_ps] == "stem", 1, 0)

ME_ps_cor_und <- apply(MEs_ps_und, 2, function(x) suppressWarnings(cor(x, ps_und, method = "spearman")))
ME_ps_p_und   <- sapply(seq_len(ncol(MEs_ps_und)), function(j) suppressWarnings(cor.test(MEs_ps_und[, j], ps_und, method = "spearman")$p.value))
ME_ps_fdr_und <- p.adjust(ME_ps_p_und, method = "BH")

ME_tissue_cor_und <- apply(MEs_ps_und, 2, function(x) suppressWarnings(cor(x, tissue_num, method = "spearman")))
ME_tissue_p_und   <- sapply(seq_len(ncol(MEs_ps_und)), function(j) suppressWarnings(cor.test(MEs_ps_und[, j], tissue_num, method = "spearman")$p.value))
ME_tissue_fdr_und <- p.adjust(ME_tissue_p_und, method = "BH")

eigengene_names <- colnames(MEs_ps_und)
eig_colors      <- sub("^uME_", "", eigengene_names)

ME_ps_tbl_und <- tibble(
  eigengene    = eigengene_names,
  color        = eig_colors,
  rho_ps       = as.numeric(ME_ps_cor_und),
  fdr_ps       = as.numeric(ME_ps_fdr_und),
  rho_tissue   = as.numeric(ME_tissue_cor_und),
  fdr_tissue   = as.numeric(ME_tissue_fdr_und)
) %>%
  left_join(modules_und %>% select(color, module_num, module_label, size), by = "color") %>%
  mutate(tissue_bias = case_when(fdr_tissue <= 0.05 & rho_tissue > 0 ~ "stem",
                                 fdr_tissue <= 0.05 & rho_tissue < 0 ~ "root",
                                 TRUE ~ "neutral")) %>%
  arrange(fdr_ps, desc(abs(rho_ps)))

write_csv(ME_ps_tbl_und, file.path(out_dir, "04_ME_vs_ps_index_und.csv"))

# Reorder eigengenes by ps_index correlation
ps_order <- ME_ps_tbl_und %>% arrange(desc(rho_ps)) %>% pull(eigengene)
moduleEigengenes_und <- moduleEigengenes_und[, ps_order, drop = FALSE]

# ==========================================
# 8) kME, Hubs, and TOM-UMAP
# ==========================================
kME_und <- signedKME(DatExpr_und, moduleEigengenes_und, outputColumnName = "ukME", 
                     corFnc = "bicor", corOptions = "maxPOutliers=0.1")

gene_annot_und <- tibble(OG = colnames(DatExpr_und), 
                         module_num = net_und$colors, 
                         module_color = labels2colors(net_und$colors)) %>%
  left_join(gene_og_mapping %>% filter(Species == "Zea_mays") %>% select(Orthogroup, GeneID), by = c("OG" = "Orthogroup"))

# Load pre-computed UMAP
umap_csv_path <- file.path(data_dir, "09_TOM_UMAP_embedding_underground.csv")
if (file.exists(umap_csv_path)) {
  cat("Loading pre-computed TOM-UMAP embedding...\n")
  df_umap <- read_csv(umap_csv_path, show_col_types = FALSE)
  if (!"module" %in% names(df_umap)) {
    df_umap <- df_umap %>% mutate(module = sub("_.*$", "", module_label))
  }
} else {
  cat("Computing TOM-UMAP embedding from scratch...\n")
}

# --- Hub highlighting update ---
cat("Updating hub highlights from OG_symbols_underground.txt...\n")
og_symbols <- read_tsv(file.path(data_dir, "OG_symbols_underground.txt"), col_types = "cc") %>%
  dplyr::rename(OG = 1, OG_symbol = 2) %>%
  distinct(OG, .keep_all = TRUE)

df_umap <- df_umap %>%
  select(-any_of(c("is_hub", "OG_symbol", "label_text_color"))) %>%
  mutate(is_hub = gene %in% og_symbols$OG) %>%
  left_join(og_symbols, by = c("gene" = "OG")) %>%
  mutate(label_text_color = ifelse(color == "turquoise", "black", "white"))

# Plot network
module_palette <- c(
  turquoise = "#64f9db", blue="blue", brown="#8E6B23", yellow="#F1C40F",
  green="#2ECC71", red="#E74C3C", black="#2C3E50", magenta="magenta",
  pink="pink", purple="#7D3C98", grey="#B3B3B3", gray="#B3B3B3"
)

p_net_focused <- ggplot(df_umap, aes(UMAP1, UMAP2)) +
  # All OGs
  geom_point(aes(fill = color, size = pmax(kME, 0.05)), color = "white", shape = 21, stroke = 0.2, alpha = 0.8) +
  # Highlighted Target Hubs
  geom_point(data = subset(df_umap, is_hub), aes(UMAP1, UMAP2, fill = color, size = pmax(kME, 0.05)),
             color = "black", shape = 21, stroke = 0.8, alpha = 1) +
  # Labels with conditional text color and no overlap limits
  geom_label_repel(
    data = subset(df_umap, is_hub), 
    aes(UMAP1, UMAP2, label = OG_symbol, fill = color, color = label_text_color),
    size = 3.5, fontface = "bold", box.padding = 0.8, point.padding = 0.3,
    label.size = 0.2, segment.color = "black", seed = 14, max.overlaps = Inf
  ) +
  scale_fill_manual(values = module_palette, guide = "none") +
  scale_color_identity() +
  scale_size_continuous(range = c(0.3, 2.5), guide = "none") +
  labs(x = "UMAP1", y = "UMAP2") +
  coord_equal(expand = TRUE) +
  theme_classic(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

ggsave(file.path(out_dir, "11_TOM_UMAP_network_und_m1_m2.png"), p_net_focused, width = 6, height = 6, dpi = 300)

# ============================================================
# 9) GO enrichment & Dynamics
# ============================================================
cat("Running GO enrichment...\n")
all_go <- map_dfr(unique(gene_annot_und$module_num), function(mnum) {
  ogs <- gene_annot_und %>% filter(module_num == mnum) %>% pull(OG) %>% unique()
  if (length(ogs) < 5) return(tibble())
  go_res <- perform_GO_enrichment_OG(ogs, og_go_df_BP, og_go_df_MF, og_go_df_CC, gmt_data)
  if (is.null(go_res) || nrow(go_res) == 0) return(tibble())
  go_res %>% mutate(module_num = mnum, logP = -log10(p.adjust))
})
write_csv(all_go, file.path(out_dir, "13_GO_enrichment_modules_all_und.csv"))

# Dynamics plot
samples_keep <- meta_underground_remobilizer %>% filter(!is.na(week)) %>% pull(sample_id)
Z_all <- scale(DatExpr_und[samples_keep, ], center = TRUE, scale = TRUE)
df_z_week <- as_tibble(Z_all, rownames = "sample_id") %>%
  pivot_longer(-sample_id, names_to = "gene", values_to = "z") %>%
  left_join(meta_underground_remobilizer %>% select(sample_id, week), by = "sample_id") %>%
  group_by(gene, week) %>%
  summarise(z_med = median(z, na.rm = TRUE), .groups = "drop")

df_plot_dyn <- df_umap %>% select(gene, UMAP1, UMAP2) %>% inner_join(df_z_week, by = "gene")

p_dyn <- ggplot(df_plot_dyn, aes(UMAP1, UMAP2, color = z_med)) +
  geom_point(size = 0.5, alpha = 0.95) +
  scale_color_gradientn(colors = rev(c("#ca1875","#ef2c69","#ff1744","#f06848","#fda25e","#fceba5","#f2f8d5","#e4f2f8")), 
                        limits = c(-2, 2), oob = scales::squish) +
  facet_wrap(~ week, nrow = 2) +
  theme_void()

ggsave(file.path(out_dir, "16_TOM_UMAP_dynamics_allweeks_und.png"), p_dyn, width = 6.5, height = 1.8, dpi = 600)

cat("\n==== Underground WGCNA pipeline complete ====\n")
