############################################################
## Reference LEAF Network (3k OGs) — WGCNA (samples × genes)
## - VST matrix from leaf_og_counts_dps (photo OGs removed)
## - Subset: remobilizers (Sorghastrum, A. gerardii, P. virgatum and 'Td_Pete')
## - Signed network (bicor), modules, eigengenes
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
  library(doParallel)
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
leaf_metadata_ps_index <- readRDS(file.path(ps_dir, "leaf_metadata_ps_index.rds"))
leaf_og_counts_dps     <- readRDS(file.path(ps_dir, "leaf_og_counts_dps.rds"))
gene_og_mapping        <- readRDS(file.path(annot_dir, "gene_og_mapping.rds"))
og_go_df_BP            <- readRDS(file.path(annot_dir, "og_go_df_BP.rds"))
og_go_df_MF            <- readRDS(file.path(annot_dir, "og_go_df_MF.rds"))
og_go_df_CC            <- readRDS(file.path(annot_dir, "og_go_df_CC.rds"))
gmt_data               <- readRDS(file.path(annot_dir, "gmt_data.rds"))

# =========================
# 2) PCA & VST of leaf remobilizers
# =========================
pca_leaf_ps_index <- plot_pca_ggplot2(
  metadata = leaf_metadata_ps_index,
  count_data = leaf_og_counts_dps,
  color_by = "ps_index",
  species_palette = species_palette,
  accession_palette = accession_palette,
  x_pc = 1, y_pc = 2
)

vst_matrix_leaf <- as.matrix(pca_leaf_ps_index$vst_matrix)
mode(vst_matrix_leaf) <- "numeric"

# ==========================================
# 3) Define leaf remobilizer reference cohort
# ==========================================
remobilizer_species <- c("Sorghastrum_nutans", "Andropogon_gerardii", "Panicum_virgatum")

df_meta_leaf <- leaf_metadata_ps_index %>%
  dplyr::mutate(sample_id = as.character(sample_id))

leaf_meta_remobilizer <- df_meta_leaf %>%
  dplyr::filter(
    tissue == "leaf",
    (species %in% remobilizer_species) |
      (species == "Tripsacum_dactyloides" & accession == "Td_Pete")
  )

common_ids_leaf <- intersect(colnames(vst_matrix_leaf), leaf_meta_remobilizer$sample_id)
stopifnot(length(common_ids_leaf) >= 30)

expr_leaf_ref <- vst_matrix_leaf[, common_ids_leaf, drop = FALSE]
leaf_meta_remobilizer <- leaf_meta_remobilizer %>%
  dplyr::filter(sample_id %in% common_ids_leaf) %>%
  dplyr::slice(match(colnames(expr_leaf_ref), sample_id))

stopifnot(identical(colnames(expr_leaf_ref), leaf_meta_remobilizer$sample_id))
message("Reference cohort (leaf_remobilizers): ", ncol(expr_leaf_ref), " samples; ", nrow(expr_leaf_ref), " OGs.")

# ==========================================
# 4) UMAP of leaf remobilizer libraries
# ==========================================
result_umap_rem <- run_umap_seurat_style(
  count_data = leaf_og_counts_dps,
  metadata   = leaf_meta_remobilizer,
  color_by   = "ps_index",
  species_palette = species_palette,
  n_top_genes = 3000,
  n_pcs = 6
)

panel_theme <- theme_classic(base_size = 10) +
  theme(plot.title = element_blank(), axis.title = element_text(size = 9),
        axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "none")

umap_df <- result_umap_rem$umap_df
p_psindex <- ggplot(umap_df, aes(UMAP1, UMAP2)) +
  geom_point(aes(fill = ps_index), shape = 21, size = 0.8,
             color = "darkslategray", stroke = 0.1, alpha = 1) +
  scale_fill_gradientn(colors = c("purple","orange","limegreen","darkgreen")) +
  coord_fixed() + labs(x = "UMAP1", y = "UMAP2", fill = "ps_index") +
  panel_theme

tiff(file.path(out_dir, "umap_leaf_rem_psindex.tiff"), width = 2.5, height = 2, units = "in", res = 600)
print(p_psindex)
dev.off()

# =====================================================
# 5) Build or load network
# =====================================================
net_leaf_path <- file.path(data_dir, "net_leaf.rds")
DatExpr_leaf_path <- file.path(data_dir, "DatExpr_leaf.rds")

if (file.exists(net_leaf_path) && file.exists(DatExpr_leaf_path)) {
  cat("Loading pre-computed WGCNA network objects...\n")
  net_leaf <- readRDS(net_leaf_path)
  DatExpr_leaf <- readRDS(DatExpr_leaf_path)

  # Recover softPower from the network object
  softPower_leaf <- net_leaf$power
  if (is.null(softPower_leaf)) {
    # Fallback: re-derive from sft
    powers <- c(1:14, seq(16, 24, 2))
    sft <- pickSoftThreshold(DatExpr_leaf, dataIsExpr = TRUE, networkType = "signed",
                             corFnc = "bicor", corOptions = list(maxPOutliers = 0.1),
                             powerVector = powers, verbose = 3)
    r2 <- -sign(sft$fitIndices[,3]) * sft$fitIndices[,2]
    cand <- which(r2 >= 0.6)
    softPower_leaf <- if (length(cand) > 0) sft$fitIndices[min(cand), 1] else sft$fitIndices[which.max(r2), 1]
  }

} else {
  cat("Building WGCNA network from scratch (this may take a while)...\n")

  # Subset to top 3,000 most variable OGs
  og_var <- matrixStats::rowVars(expr_leaf_ref)
  ord <- order(og_var, decreasing = TRUE)
  K <- min(3000L, nrow(expr_leaf_ref))
  sel <- ord[seq_len(K)]
  expr_leaf_ref_3k <- expr_leaf_ref[sel, , drop = FALSE]

  DatExpr_leaf <- t(expr_leaf_ref_3k)
  stopifnot(identical(rownames(DatExpr_leaf), colnames(expr_leaf_ref_3k)))

  # Drop zero-MAD genes
  mad_vals <- matrixStats::colMads(DatExpr_leaf)
  DatExpr_leaf <- DatExpr_leaf[, mad_vals > 0, drop = FALSE]

  gsg <- goodSamplesGenes(DatExpr_leaf, verbose = 3)
  if (!gsg$allOK) {
    DatExpr_leaf <- DatExpr_leaf[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
    leaf_meta_remobilizer <- dplyr::filter(leaf_meta_remobilizer, sample_id %in% rownames(DatExpr_leaf))
  }

  leaf_meta_remobilizer <- leaf_meta_remobilizer[match(rownames(DatExpr_leaf), leaf_meta_remobilizer$sample_id), , drop = FALSE]
  stopifnot(identical(rownames(DatExpr_leaf), leaf_meta_remobilizer$sample_id))

  # Soft-threshold selection
  powers <- c(1:14, seq(16, 24, 2))
  sft <- pickSoftThreshold(DatExpr_leaf, dataIsExpr = TRUE, networkType = "signed",
                           corFnc = "bicor", corOptions = list(maxPOutliers = 0.1),
                           powerVector = powers, verbose = 5)

  pdf(file.path(out_dir, "01_soft_threshold_leaf.pdf"), width = 6, height = 4)
  par(mfrow = c(1, 2))
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2",
       type="n", main="Scale independence")
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=0.8)
  abline(h=0.6, col="red", lty=2)
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)", ylab="Mean connectivity", type="n", main="Mean connectivity")
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.8)
  dev.off()

  r2 <- -sign(sft$fitIndices[,3]) * sft$fitIndices[,2]
  cand <- which(r2 >= 0.6)
  softPower_leaf <- if (length(cand) > 0) sft$fitIndices[min(cand), 1] else sft$fitIndices[which.max(r2), 1]
  message("Selected soft-threshold β (leaf): ", softPower_leaf)

  # Build network
  net_leaf <- blockwiseModules(
    DatExpr_leaf,
    power             = softPower_leaf,
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

  # Save for future runs
  saveRDS(net_leaf, net_leaf_path)
  saveRDS(DatExpr_leaf, DatExpr_leaf_path)
}

# Ensure metadata matches DatExpr
leaf_meta_remobilizer <- leaf_meta_remobilizer %>%
  dplyr::filter(sample_id %in% rownames(DatExpr_leaf))
leaf_meta_remobilizer <- leaf_meta_remobilizer[match(rownames(DatExpr_leaf), leaf_meta_remobilizer$sample_id), , drop = FALSE]

moduleColors_leaf     <- labels2colors(net_leaf$colors)
moduleEigengenes_leaf <- orderMEs(net_leaf$MEs)

# Dendrogram
pdf(file.path(out_dir, "02_gene_dendrogram_modules_leaf.pdf"), width = 3, height = 3)
plotDendroAndColors(net_leaf$dendrograms[[1]],
                    moduleColors_leaf[net_leaf$blockGenes[[1]]],
                    "modules", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main = "")
dev.off()

# ===================================================
# 6) Module eigengene correlation with ps_index
# ===================================================
# Rename MEs to ME_color and reorder by module size
me_nums_all   <- as.integer(sub("^ME_?", "", colnames(moduleEigengenes_leaf)))
me_colors_all <- WGCNA::labels2colors(me_nums_all)
colnames(moduleEigengenes_leaf) <- paste0("ME_", me_colors_all)

sizes_by_color <- table(WGCNA::labels2colors(net_leaf$colors))
sizes_by_color <- as.integer(sizes_by_color)[match(names(sizes_by_color), names(sizes_by_color))]
names(sizes_by_color) <- names(table(WGCNA::labels2colors(net_leaf$colors)))

me_colors <- sub("^ME_?", "", colnames(moduleEigengenes_leaf))
me_sizes  <- sizes_by_color[me_colors]
ord_idx   <- order(me_sizes, decreasing = TRUE, na.last = TRUE)
moduleEigengenes_leaf <- moduleEigengenes_leaf[, ord_idx, drop = FALSE]

message("==== Leaf network complete ====")
message("Samples: ", nrow(DatExpr_leaf), "; Genes: ", ncol(DatExpr_leaf))
message("Soft power selected: ", softPower_leaf)
message("Modules detected: ", length(unique(net_leaf$colors)))

# Module summary
modules_leaf <- as.data.frame(table(moduleColors_leaf)) %>%
  dplyr::rename(color = moduleColors_leaf, size = Freq) %>%
  arrange(desc(size)) %>%
  mutate(
    module_num = ifelse(color %in% c("gray", "grey"), 0,
                        as.numeric(factor(color, levels = setdiff(unique(color), c("gray", "grey"))))),
    module_label = paste0("m", module_num, "_", color)
  ) %>%
  arrange(module_num)

print(modules_leaf)

# Module sizes bar plot
pdf(file.path(out_dir, "03_module_sizes_leaf.pdf"), width = 2, height = 3)
ggplot(modules_leaf, aes(x = size, y = reorder(module_label, size), fill = color)) +
  geom_col(color = "black", linewidth = 0.3) +
  scale_fill_identity() +
  scale_x_continuous(position = "top") +
  labs(x = "size (OGs)", y = "co-expression module") +
  theme_classic(base_size = 10) +
  theme(axis.text.y = element_text(size = 9), panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(), plot.margin = margin(5, 5, 5, 5))
dev.off()

# --- Correlation analysis ---
keep_ps <- !is.na(leaf_meta_remobilizer$ps_index)
MEs_ps_leaf <- moduleEigengenes_leaf[keep_ps, , drop = FALSE]
ps_leaf     <- as.numeric(leaf_meta_remobilizer$ps_index[keep_ps])

ME_ps_cor_leaf <- apply(MEs_ps_leaf, 2, function(x) suppressWarnings(cor(x, ps_leaf, method = "spearman")))
ME_ps_p_leaf   <- sapply(seq_len(ncol(MEs_ps_leaf)), function(j) suppressWarnings(cor.test(MEs_ps_leaf[, j], ps_leaf, method = "spearman")$p.value))
ME_ps_fdr_leaf <- p.adjust(ME_ps_p_leaf, method = "BH")

eigengene_names <- colnames(MEs_ps_leaf)
eig_colors      <- sub("^ME_?", "", eigengene_names)

col_num_map  <- tibble(module_num = unique(net_leaf$colors)) |>
  mutate(module_color = WGCNA::labels2colors(module_num)) |>
  distinct(module_color, module_num)
num_by_color <- setNames(col_num_map$module_num, col_num_map$module_color)

module_num   <- unname(num_by_color[eig_colors])
module_short <- paste0("m", module_num)
module_label <- paste0(module_short, "_", eig_colors)

size_lookup <- setNames(modules_leaf$size, modules_leaf$color)
sizes_vec   <- unname(size_lookup[eig_colors])

ME_ps_tbl_leaf <- tibble(
  eigengene = eigengene_names, module_num = module_num, module = module_short,
  color = eig_colors, module_label = module_label,
  rho_ps = as.numeric(ME_ps_cor_leaf), p_ps = as.numeric(ME_ps_p_leaf),
  fdr_ps = as.numeric(ME_ps_fdr_leaf), size = sizes_vec
) |> dplyr::arrange(fdr_ps, dplyr::desc(abs(rho_ps)))

readr::write_csv(ME_ps_tbl_leaf, file.path(out_dir, "04_ME_vs_ps_index_leaf.csv"))

# Reorder by ps_index correlation
ps_order <- ME_ps_tbl_leaf %>% dplyr::arrange(dplyr::desc(rho_ps)) %>% dplyr::pull(eigengene)
moduleEigengenes_leaf <- moduleEigengenes_leaf[, ps_order, drop = FALSE]
MEs_ps_leaf <- MEs_ps_leaf[, ps_order, drop = FALSE]
ME_ps_cor_leaf <- ME_ps_cor_leaf[ps_order]

# ME correlation heatmap
MEcorr_leaf <- cor(moduleEigengenes_leaf, use = "pairwise.complete.obs")
MEcorr_leaf <- MEcorr_leaf[ps_order, ps_order]

pdf(file.path(out_dir, "05_ME_correlation_heatmap_leaf.pdf"), width = 5, height = 4.2)
labeledHeatmap(
  Matrix = MEcorr_leaf, xLabels = ps_order, yLabels = ps_order,
  colors = blueWhiteRed(50), zlim = c(-1, 1), textMatrix = NULL,
  main = "eigengene correlation", xLabelsAngle = 90
)
dev.off()

# ==========================================
# 7) kME, annotations, exports
# ==========================================
kME_leaf <- signedKME(DatExpr_leaf, moduleEigengenes_leaf,
                      outputColumnName = "kME",
                      corFnc = "bicor", corOptions = "maxPOutliers=0.1")

gene_annot_leaf <- tibble(
  OG = colnames(DatExpr_leaf),
  module_num = net_leaf$colors,
  module_color = labels2colors(net_leaf$colors)
)

kME_tbl_leaf <- gene_annot_leaf %>%
  dplyr::bind_cols(as_tibble(kME_leaf, .name_repair = "minimal")) %>%
  dplyr::relocate(OG, module_color, module_num)

ME_df_leaf <- as_tibble(moduleEigengenes_leaf, rownames = "sample_id")

gene_annot_leaf_zea <- gene_annot_leaf %>%
  left_join(
    gene_og_mapping %>% filter(Species == "Zea_mays") %>% select(Orthogroup, GeneID),
    by = c("OG" = "Orthogroup")
  )

readr::write_csv(gene_annot_leaf,     file.path(out_dir, "06a_gene_module_assignments_leaf.csv"))
readr::write_csv(gene_annot_leaf_zea, file.path(out_dir, "06b_gene_module_assignments_with_Zea.csv"))
readr::write_csv(kME_tbl_leaf,        file.path(out_dir, "07_gene_kME_table_leaf.csv"))
readr::write_csv(ME_df_leaf,          file.path(out_dir, "08_module_eigengenes_samples_leaf.csv"))

message("==== kME and annotations exported ====")

# ==========================================
# 8) Hub gene export
# ==========================================
hubs_per_module <- 4L
include_grey    <- FALSE

kme_cols_color <- grep("^kME_", names(kME_tbl_leaf), value = TRUE)
kme_long <- kME_tbl_leaf %>%
  tidyr::pivot_longer(cols = all_of(kme_cols_color), names_to = "kME_col", values_to = "kME") %>%
  mutate(module_color_from_col = sub("^kME_", "", kME_col))

kme_self <- kme_long %>%
  filter(module_color == module_color_from_col) %>%
  { if (include_grey) . else filter(., !module_color %in% c("grey","gray")) } %>%
  select(OG, module_color, module_num, kME)

hubs_ranked <- kme_self %>%
  group_by(module_color, module_num) %>%
  arrange(desc(abs(kME)), .by_group = TRUE) %>%
  mutate(hub_rank = dplyr::row_number()) %>%
  slice_head(n = hubs_per_module) %>%
  ungroup()

hubs_with_zea_expanded <- hubs_ranked %>%
  left_join(gene_annot_leaf_zea %>% select(OG, GeneID), by = "OG") %>%
  arrange(module_color, hub_rank, OG, GeneID) %>%
  mutate(module_label = paste0("m", module_num, "_", module_color)) %>%
  relocate(module_label, .after = module_color)

hubs_with_zea_collapsed <- hubs_with_zea_expanded %>%
  group_by(module_color, module_num, OG, kME, hub_rank) %>%
  summarise(Zea_GeneIDs = paste(sort(unique(na.omit(GeneID))), collapse = ";"), .groups = "drop") %>%
  arrange(module_color, hub_rank, OG) %>%
  mutate(module_label = paste0("m", module_num, "_", module_color)) %>%
  relocate(module_label, .after = module_color)

readr::write_csv(hubs_with_zea_expanded,  file.path(out_dir, "10_hubs_with_Zea_expanded.csv"))
readr::write_csv(hubs_with_zea_collapsed, file.path(out_dir, "10_hubs_with_Zea_collapsed.csv"))
message("✓ Hub gene tables exported")

# ===================================================================
# 9) TOM - UMAP network
# ===================================================================
cat("Computing adjacency and TOM (this is the slow step)...\n")

mod_map <- tibble(
  module_num = sort(unique(net_leaf$colors)),
  color      = labels2colors(sort(unique(net_leaf$colors)))
) %>%
  mutate(module = paste0("m", module_num), module_label = paste0(module, "_", color))

genes <- colnames(DatExpr_leaf)

A <- adjacency(DatExpr_leaf, type = "signed", power = softPower_leaf,
               corFnc = "bicor", corOptions = list(maxPOutliers = 0.1))
TOM <- WGCNA::TOMsimilarity(A, TOMType = "signed")
dimnames(TOM) <- list(genes, genes)

module_colors <- labels2colors(net_leaf$colors)
names(module_colors) <- genes
mods_keep <- setdiff(sort(unique(module_colors)), c("grey","gray"))

# kME for hub selection
kME_leaf <- signedKME(DatExpr_leaf, moduleEigengenes_leaf,
                      outputColumnName = "kME",
                      corFnc = "bicor", corOptions = "maxPOutliers=0.1")
kME_mat <- as.matrix(kME_leaf)

# Pick hubs per module
n_hubs <- 4
hubs <- character(0)
for (m in mods_keep) {
  gmod <- names(module_colors)[module_colors == m]
  if (!length(gmod)) next
  kcol <- paste0("kME_", m)
  if (kcol %in% colnames(kME_mat) && any(!is.na(kME_mat[gmod, kcol]))) {
    ord <- order(kME_mat[gmod, kcol], decreasing = TRUE, na.last = NA)
    hubs <- c(hubs, head(gmod[ord], n_hubs))
  } else {
    deg <- rowSums(TOM[gmod, gmod, drop = FALSE], na.rm = TRUE)
    hubs <- c(hubs, head(names(sort(deg, decreasing = TRUE)), n_hubs))
  }
}
hubs <- unique(na.omit(hubs))

# --- UMAP embedding: load pre-computed or compute fresh ---
umap_csv_path <- file.path(data_dir, "09_TOM_UMAP_embedding_leaf.csv")

if (file.exists(umap_csv_path)) {
  cat("Loading pre-computed TOM-UMAP embedding for reproducibility...\n")
  df_umap <- readr::read_csv(umap_csv_path, show_col_types = FALSE)
  # Reconstruct 'module' column if missing (CSV only stores module_label)
  if (!"module" %in% names(df_umap)) {
    df_umap <- df_umap %>% mutate(module = sub("_.*$", "", module_label))
  }
} else {
  cat("Computing TOM-UMAP embedding from scratch...\n")
  keep <- module_colors %in% mods_keep
  genes_keep <- names(module_colors)[keep]
  feature_hubs <- intersect(hubs, genes_keep)

  feature_mat <- TOM[genes_keep, if (length(feature_hubs)) feature_hubs else genes_keep, drop = FALSE]

  set.seed(14)
  emb <- uwot::umap(X = feature_mat, n_neighbors = 50, min_dist = 0.4,
                     n_threads = 12, metric = "cosine", verbose = TRUE)
  ## Note: threading, seed, and dependencies may change UMAP layout across machines
  ## but clustering structure should be retained

  col2num <- setNames(mod_map$module_num, mod_map$color)
  df_umap <- tibble(
    gene = rownames(feature_mat), UMAP1 = emb[,1], UMAP2 = emb[,2],
    color = unname(module_colors[rownames(feature_mat)]),
    module_num = unname(col2num[color])
  ) %>%
    mutate(module = paste0("m", module_num),
           module_label = paste0(module, "_", color),
           is_hub = gene %in% colnames(feature_mat))

  # own-module kME (scaled 0..1 within module)
  df_umap$kME <- NA_real_
  for (m in mods_keep) {
    idx <- which(df_umap$color == m)
    kcol <- paste0("kME_", m)
    if (!(kcol %in% colnames(kME_mat))) next
    df_umap$kME[idx] <- kME_mat[df_umap$gene[idx], kcol]
  }
  df_umap <- df_umap %>%
    group_by(color) %>%
    mutate(kME = (kME - min(kME, na.rm = TRUE)) /
             (max(kME, na.rm = TRUE) - min(kME, na.rm = TRUE))) %>%
    ungroup()
  df_umap$kME[!is.finite(df_umap$kME)] <- 0.05

  write_csv(df_umap %>% select(gene, UMAP1, UMAP2, color, module_num, module_label, is_hub, kME),
            file.path(out_dir, "09_TOM_UMAP_embedding_leaf.csv"))
}

# --- TOM-UMAP point plot ---
module_palette <- c(
  turquoise = "#1ABC9C", blue="blue", brown="#8E6B23", yellow="#F1C40F",
  green="#2ECC71", red="#E74C3C", black="#2C3E50", magenta="magenta",
  pink="pink", purple="#7D3C98", grey="#B3B3B3", gray="#B3B3B3"
)
mods_present <- sort(unique(na.omit(df_umap$color)))
pal_used <- module_palette[mods_present]

p_umap <- ggplot(df_umap, aes(UMAP1, UMAP2)) +
  geom_point(aes(color = color, size = pmax(kME, 0.05)), alpha = 0.9, stroke = 0.1, shape = 16) +
  geom_point(data = subset(df_umap, is_hub),
             aes(UMAP1, UMAP2, fill = color, size = pmax(kME, 0.05)),
             color = "black", shape = 21, stroke = 0.6, alpha = 1) +
  scale_color_manual(values = pal_used, guide = "none") +
  scale_fill_manual(values = pal_used, guide = "none") +
  scale_size_continuous(range = c(0.1, 2.4), guide = "none") +
  coord_equal(expand = FALSE) +
  theme_minimal(base_size = 10)

ggsave(file.path(out_dir, "09_TOM_UMAP_points_leaf.png"), p_umap, width = 4, height = 4, dpi = 300)

# ===================================================================
# 10) TOM-UMAP network with edges and OG symbol labels
# ===================================================================
og_symbols <- readr::read_tsv(file.path(data_dir, "OG_symbols_m2_blue.txt"), col_types = "cc") %>%
  dplyr::rename(OG = 1, OG_symbol = 2) %>%
  mutate(OG = str_trim(OG), OG_symbol = na_if(str_trim(OG_symbol), "")) %>%
  distinct(OG, .keep_all = TRUE)

df_umap <- df_umap %>%
  left_join(og_symbols, by = c("gene" = "OG")) %>%
  mutate(label = coalesce(OG_symbol, gene))

# Sparse TOM edges
top_percent <- 0.05; edge_prop <- 0.03; max_edges <- 150000; within_only <- FALSE
gnames <- rownames(TOM)
ut <- upper.tri(TOM, diag = FALSE)
w  <- TOM[ut]

thr <- stats::quantile(w, probs = 1 - top_percent, na.rm = TRUE)
idx <- which(ut & TOM >= thr, arr.ind = TRUE)
edges <- data.frame(from = gnames[idx[,1]], to = gnames[idx[,2]], weight = TOM[idx], stringsAsFactors = FALSE)

mod_map_vec <- setNames(df_umap$module, df_umap$gene)
col_map <- setNames(df_umap$color, df_umap$gene)
edges$mod_from <- mod_map_vec[edges$from]
edges$mod_to   <- mod_map_vec[edges$to]

set.seed(14)
n_keep <- min(max_edges, ceiling(nrow(edges) * edge_prop))
if (nrow(edges) > n_keep) edges <- edges[sample.int(nrow(edges), n_keep), , drop = FALSE]

edges$edge_col <- ifelse(edges$mod_from == edges$mod_to, col_map[edges$from], "grey80")
w99 <- stats::quantile(edges$weight, 0.99, na.rm = TRUE)
edges$alpha <- pmin(edges$weight / w99, 1)

xy <- df_umap[, c("gene","UMAP1","UMAP2")]
ix_from <- match(edges$from, xy$gene)
ix_to   <- match(edges$to, xy$gene)
edges_plot <- data.frame(
  from = edges$from, to = edges$to,
  x = xy$UMAP1[ix_from], y = xy$UMAP2[ix_from],
  xend = xy$UMAP1[ix_to], yend = xy$UMAP2[ix_to],
  alpha = edges$alpha, edge_col = edges$edge_col, stringsAsFactors = FALSE
)
edges_plot <- edges_plot[complete.cases(edges_plot[, c("x","y","xend","yend")]), , drop = FALSE]

# Update hub flag with curated symbols
df_umap <- df_umap %>% mutate(is_hub = gene %in% og_symbols$OG)

p_net_focused <- ggplot(df_umap, aes(UMAP1, UMAP2)) +
  geom_segment(data = edges_plot,
               aes(x = x, y = y, xend = xend, yend = yend, alpha = alpha),
               color = I(edges_plot$edge_col), linewidth = 0.25,
               lineend = "round", show.legend = FALSE, inherit.aes = FALSE) +
  geom_point(aes(fill = color, size = pmax(kME, 0.05)),
             color = "white", shape = 21, stroke = 0.2, alpha = 0.8) +
  geom_point(data = subset(df_umap, is_hub),
             aes(UMAP1, UMAP2, fill = color, size = pmax(kME, 0.05)),
             color = "black", shape = 21, stroke = 0.8, alpha = 1) +
  geom_label_repel(data = subset(df_umap, is_hub),
                   aes(UMAP1, UMAP2, label = OG_symbol, fill = color),
                   color = "white", size = 3.5, fontface = "bold",
                   box.padding = 0.8, point.padding = 0.3,
                   label.size = 0.2, segment.color = "black", seed = 14,
                   max.overlaps = Inf) +
  scale_fill_manual(values = pal_used, guide = "none") +
  scale_size_continuous(range = c(0.3, 2.5), guide = "none") +
  scale_alpha_continuous(range = c(0.2, 0.8), guide = "none") +
  labs(x = "UMAP1", y = "UMAP2") +
  coord_equal(expand = TRUE) +
  theme_classic(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

ggsave(file.path(out_dir, "11_TOM_UMAP_network_leaf_m2_blue.png"),
       p_net_focused, width = 6, height = 6, dpi = 300)

# =========================================
# 11) Eigengenes vs ps_index faceted plot
# =========================================
order_levels <- tibble(
  eigengene = colnames(MEs_ps_leaf), rho_ps = as.numeric(ME_ps_cor_leaf)
) %>%
  mutate(module_color = sub("^ME_", "", eigengene)) %>%
  left_join(modules_leaf %>% select(color, module_label), by = c("module_color" = "color")) %>%
  arrange(desc(rho_ps)) %>% pull(module_label) %>% unique()

ME_long <- MEs_ps_leaf %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample_id") %>%
  tidyr::pivot_longer(cols = -sample_id, names_to = "module", values_to = "ME_value") %>%
  dplyr::left_join(leaf_meta_remobilizer %>% dplyr::select(sample_id, ps_index), by = "sample_id") %>%
  dplyr::filter(!is.na(ps_index)) %>%
  dplyr::mutate(module_color = sub("^ME_?", "", module)) %>%
  dplyr::left_join(modules_leaf %>% dplyr::select(module_color = color, module_num, module_label), by = "module_color") %>%
  dplyr::mutate(module_label = factor(module_label, levels = intersect(order_levels, unique(module_label))))

p_ME_all_psindex <- ggplot(ME_long, aes(x = ps_index, y = ME_value)) +
  geom_point(color = "gray60", alpha = 0.5, size = 0.3) +
  geom_smooth(aes(color = module_color), method = "loess", se = TRUE, formula = y ~ x,
              fill = "slategray2", linewidth = 0.7) +
  scale_x_reverse(limits = c(1, 0)) +
  facet_wrap(~ module_label, ncol = 11, scales = "fixed") +
  labs(x = "ps_index", y = "module eigengene") +
  scale_color_identity() +
  theme_classic(base_size = 9) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1),
        panel.spacing = unit(3, "pt"), legend.position = "none")

tiff(file.path(out_dir, "12_MEs_psindex.tiff"), width = 7.6, height = 1.4, units = "in", res = 600)
print(p_ME_all_psindex)
dev.off()

# ============================================================
# 12) GO enrichment for ALL modules
# ============================================================
cat("Running GO enrichment for all modules...\n")

mod_map_go <- gene_annot_leaf %>%
  distinct(module_num, module_color, OG) %>%
  dplyr::count(module_num, module_color, name = "module_size") %>%
  arrange(desc(module_size)) %>%
  mutate(module_label = paste0("m", module_num, "_", module_color))

module_order_go <- mod_map_go$module_num

ogs_by_module <- gene_annot_leaf %>%
  distinct(module_num, OG) %>%
  group_by(module_num) %>%
  summarise(ogs = list(unique(OG)), .groups = "drop")

parse_ratio <- function(s) {
  if (is.na(s) || s == "" || !grepl("/", s)) return(NA_real_)
  p <- strsplit(s, "/", fixed = TRUE)[[1]]
  suppressWarnings(as.numeric(p[1]) / as.numeric(p[2]))
}

prep_go_df <- function(go_df, module_num) {
  if (is.null(go_df) || !nrow(go_df)) return(tibble())
  nm <- names(go_df)
  grab <- function(x) if (x %in% nm) go_df[[x]] else NA
  tibble(
    module_num = module_num, ID = grab("ID"),
    term_desc = dplyr::coalesce(grab("desc"), grab("Description")),
    Ontology = dplyr::coalesce(grab("Ontology"), grab("ontology")),
    pvalue = suppressWarnings(as.numeric(grab("pvalue"))),
    p.adjust = suppressWarnings(as.numeric(grab("p.adjust"))),
    Count = suppressWarnings(as.numeric(grab("Count"))),
    GeneRatio = grab("GeneRatio"), BgRatio = grab("BgRatio"),
    RichFactor = suppressWarnings(as.numeric(grab("RichFactor")))
  ) %>%
    mutate(
      GeneRatio_num = vapply(GeneRatio, parse_ratio, numeric(1)),
      BgRatio_num = vapply(BgRatio, parse_ratio, numeric(1)),
      fold_enrichment = ifelse(is.finite(GeneRatio_num) & is.finite(BgRatio_num) & BgRatio_num > 0,
                               GeneRatio_num / BgRatio_num, NA_real_),
      logP = -log10(p.adjust),
      label = paste0(term_desc, " (", ID, ")")
    ) %>%
    filter(!is.na(Ontology), is.finite(logP))
}

all_go <- purrr::map_dfr(module_order_go, function(mnum) {
  ogs <- ogs_by_module %>% filter(module_num == mnum) %>% pull(ogs) %>% .[[1]]
  if (length(ogs) < 5) return(tibble())
  go_res <- perform_GO_enrichment_OG(
    og_vector = ogs, gmt_bp_df = og_go_df_BP, gmt_mf_df = og_go_df_MF,
    gmt_cc_df = og_go_df_CC, gmt_data_df = gmt_data
  )
  prep_go_df(go_res, mnum)
})

all_go <- all_go %>%
  left_join(mod_map_go %>% select(module_num, module_color, module_size, module_label), by = "module_num") %>%
  mutate(
    module_simple = factor(paste0("ME", module_num), levels = paste0("ME", module_order_go)),
    module_label = factor(module_label, levels = mod_map_go$module_label)
  )

readr::write_csv(all_go, file.path(out_dir, "15_GO_enrichment_modules_all.csv"))

# GO bar plot for m2_blue
mods_keep_go <- c("m2_blue")
top_n <- 8
ont_keep <- c("BP")

go_plot_df <- all_go %>% filter(Ontology %in% ont_keep, module_label %in% mods_keep_go)

top_by_module <- go_plot_df %>%
  group_by(module_label, module_num) %>%
  arrange(desc(logP), desc(Count), desc(RichFactor), .by_group = TRUE) %>%
  slice_head(n = top_n) %>%
  ungroup() %>%
  mutate(
    label_x = factor(label, levels = rev(unique(label))),
    col_plot = ifelse(p.adjust <= 0.05, module_color, "gray80"),
    module_label = factor(as.character(module_label), levels = as.character(mods_keep_go))
  )

top_by_module_bar <- top_by_module %>%
  mutate(
    term_short = str_trunc(term_desc, 62),
    txt_col = if_else(col_plot %in% c("grey80","gray80","#B3B3B3","turquoise","blue"), "black", "white")
  ) %>%
  group_by(module_label) %>%
  mutate(max_logP = max(logP, na.rm = TRUE), y_lab = pmax(0, logP - 0.05 * max_logP)) %>%
  ungroup()

p_go_bars <- ggplot(top_by_module_bar, aes(x = label_x, y = logP)) +
  geom_col(aes(fill = "white"), width = 0.7, color = "blue", linewidth = 0.4, alpha = 0.7) +
  geom_text(aes(y = 0.03 * max_logP, label = term_short, fontface="bold"),
            color = "black", hjust = 0, size = 2.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "orange", alpha = 0.6) +
  scale_fill_identity() + coord_flip() +
  facet_wrap(~ module_label, nrow = 1, ncol = 2, scales = "free", drop = TRUE) +
  labs(y = expression(-log[10]("p.adj")), x = NULL) +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 10),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "none", panel.spacing = unit(6, "pt"), plot.margin = margin(5, 5, 5, 5))

ggsave(file.path(out_dir, "15_GO_bar_plot_m2_blue.png"), p_go_bars, width = 3.1, height = 2.5, dpi = 600)

# ============================================================
# 13) Network expression dynamics (weekly z-score on TOM-UMAP)
# ============================================================
samples_keep <- leaf_meta_remobilizer %>% filter(!is.na(week)) %>% pull(sample_id)
DatExpr_use <- DatExpr_leaf[samples_keep, , drop = FALSE]
Z_all <- scale(DatExpr_use, center = TRUE, scale = TRUE)

df_z_week <- as_tibble(Z_all, rownames = "sample_id") %>%
  pivot_longer(-sample_id, names_to = "gene", values_to = "z") %>%
  left_join(leaf_meta_remobilizer %>% select(sample_id, week), by = "sample_id") %>%
  filter(!is.na(week)) %>%
  group_by(gene, week) %>%
  summarise(z_med = median(z, na.rm = TRUE), .groups = "drop")

week_levels <- sort(unique(df_z_week$week))
df_plot <- df_umap %>%
  select(gene, UMAP1, UMAP2) %>%
  inner_join(df_z_week, by = "gene") %>%
  mutate(week = factor(week, levels = week_levels), z_capped = pmax(pmin(z_med, 3), -3))

custom_colors <- rev(c("#ca1875","#ef2c69","#ff1744","#f06848",
                       "#fda25e","#fceba5","#f2f8d5","#e4f2f8"))

color_scale <- scale_color_gradientn(
  colors = custom_colors, limits = c(-2, 2), oob = scales::squish,
  guide = guide_colorbar(barwidth = 0.6, barheight = 5, title.position = "top")
)

p_umap_dyn <- ggplot() +
  { if (exists("edges_plot")) geom_segment(
    data = edges_plot,
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE, linewidth = 0.1, alpha = 0.05, color = "grey60", lineend = "round"
  ) else NULL } +
  geom_point(data = df_plot, aes(UMAP1, UMAP2, color = z_capped),
             size = 0.5, alpha = 0.95, stroke = 0) +
  color_scale +
  coord_equal(expand = FALSE) +
  facet_wrap(~ week, nrow = 2, ncol = 7, drop = TRUE) +
  labs(x = NULL, y = NULL, color = "expression") +
  theme_void(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 10, face = "bold"),
        legend.position = "right", plot.margin = margin(5, 5, 5, 5))

ggsave(file.path(out_dir, "16_TOM_UMAP_dynamics_allweeks.png"),
       p_umap_dyn, width = 6, height = 2, dpi = 600)

# Save key objects for downstream module preservation analysis
saveRDS(net_leaf, file.path(data_dir, "net_leaf.rds"))
saveRDS(DatExpr_leaf, file.path(data_dir, "DatExpr_leaf.rds"))

cat("\n==== WGCNA leaf pipeline complete ====\n")
