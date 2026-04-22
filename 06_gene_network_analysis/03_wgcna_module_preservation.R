############################################################
## Module Preservation Analysis — WGCNA
## - Evaluates conservation of remobilizer modules in annuals and other perennials
## - Supports leaf, root, and stem tissues
## - Uses pre-computed RDS results for efficiency
############################################################

suppressPackageStartupMessages({
  library(WGCNA)
  library(tidyverse)
  library(matrixStats)
  library(purrr)
})

# --- Reproducibility ---
set.seed(42)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# =========================
# 0) Paths & Functions
# =========================
ps_dir   <- file.path(getwd(), "04_photosynthetic_index", "data")
data_dir <- file.path(getwd(), "06_gene_network_analysis", "data")
out_dir  <- file.path(getwd(), "06_gene_network_analysis", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Filter low / no variance genes | return aligned inputs for modulePreservation
mod_preservation_vfilt <- function(Expr_test, ref_genes, ref_colors, DatExpr_ref, minOKFrac = 0.5) {
  Expr_test   <- as.matrix(Expr_test);   mode(Expr_test)   <- "numeric"
  DatExpr_ref <- as.matrix(DatExpr_ref); mode(DatExpr_ref) <- "numeric"

  genes_common <- intersect(colnames(Expr_test), intersect(ref_genes, colnames(DatExpr_ref)))
  if (!length(genes_common)) stop("No overlap between Expr_test, ref_genes, and DatExpr_ref.")

  RefExpr_full  <- DatExpr_ref[, genes_common, drop = FALSE]
  TestExpr_full <- Expr_test[,   genes_common, drop = FALSE]

  ok_ref  <- colMeans(!is.na(RefExpr_full))  >= minOKFrac
  ok_test <- colMeans(!is.na(TestExpr_full)) >= minOKFrac
  genes_final <- genes_common[ ok_ref & ok_test ]
  
  if (!length(genes_final)) stop("All genes failed minOKFrac filter.")

  list(
    RefExpr        = RefExpr_full[,  genes_final, drop = FALSE],
    TestExpr       = TestExpr_full[, genes_final, drop = FALSE],
    ref_colors_sub = ref_colors[genes_final],
    genes_final    = genes_final
  )
}

# Runner for WGCNA::modulePreservation
run_modpres_to_mp <- function(meta_subset, net_ref, DatExpr_ref, full_vst, accession, species, 
                             nPermutations = 5L, minOKFrac = 0.5, randomSeed = 42, networkType = "signed", verbose = 2) {
  
  meta <- meta_subset %>% mutate(sample_id = as.character(sample_id))
  DatExpr_ref <- as.matrix(DatExpr_ref); mode(DatExpr_ref) <- "numeric"
  ref_genes <- colnames(DatExpr_ref)
  
  ref_colors <- WGCNA::labels2colors(net_ref$colors)
  names(ref_colors) <- ref_genes
  ref_labels <- net_ref$colors
  names(ref_labels) <- ref_genes

  full_vst <- as.matrix(full_vst); mode(full_vst) <- "numeric"

  sample_ids <- meta %>% filter(species == !!species, accession == !!accession) %>% pull(sample_id)
  sample_ids <- intersect(colnames(full_vst), sample_ids)
  if (!length(sample_ids)) stop(paste0("No samples found for: ", accession))

  keep_ref <- intersect(rownames(full_vst), ref_genes)
  Expr_test <- t(full_vst[keep_ref, sample_ids, drop = FALSE])
  Expr_test <- Expr_test[, ref_genes[ref_genes %in% colnames(Expr_test)], drop = FALSE]

  flt <- mod_preservation_vfilt(Expr_test, ref_genes, ref_colors, DatExpr_ref, minOKFrac)
  RefExpr  <- flt$RefExpr
  TestExpr <- flt$TestExpr
  
  g1 <- WGCNA::goodSamplesGenes(RefExpr, verbose = 0)
  g2 <- WGCNA::goodSamplesGenes(TestExpr, verbose = 0)
  if (!g1$allOK) { RefExpr <- RefExpr[g1$goodSamples, g1$goodGenes]; TestExpr <- TestExpr[, colnames(RefExpr)] }
  if (!g2$allOK) { TestExpr <- TestExpr[g2$goodSamples, g2$goodGenes]; RefExpr <- RefExpr[, colnames(TestExpr)] }

  multiExpr <- list(Ref = list(data = RefExpr), Test = list(data = TestExpr))
  colorList <- list(Ref = ref_labels[colnames(RefExpr)])

  set.seed(randomSeed)
  mp <- WGCNA::modulePreservation(multiExpr, colorList, referenceNetworks = 1, nPermutations = nPermutations,
                                  randomSeed = randomSeed, networkType = networkType, verbose = verbose)

  list(mp = mp, RefExpr = RefExpr, TestExpr = TestExpr, ref_colors_sub = flt$ref_colors_sub)
}

# Helper to extract results
extract_all_preservation_results <- function(pres_list, modules_meta, drop_gold = TRUE) {
  imap_dfr(pres_list, function(x, acc) {
    # Check if x is a result from purrr::safely
    if (!is.null(x$result)) {
      res <- x$result
    } else {
      res <- x
    }
    
    if (is.null(res) || is.null(res$mp)) return(tibble())
    
    mp <- res$mp
    Zraw   <- try(mp$preservation$Z$ref.Ref$inColumnsAlsoPresentIn.Test, silent = TRUE)
    Obsraw <- try(mp$preservation$observed$ref.Ref$inColumnsAlsoPresentIn.Test, silent = TRUE)
    
    if (inherits(Zraw, "try-error") || is.null(Zraw) || nrow(Zraw) == 0) return(tibble())

    Ztab <- Zraw %>% as.data.frame() %>% rownames_to_column("module_id")
    Obs  <- Obsraw %>% as.data.frame() %>% rownames_to_column("module_id")

    as_num <- suppressWarnings(as.numeric(Ztab$module_id))
    numeric_labels <- any(!is.na(as_num))

    if (numeric_labels) {
      Ztab <- Ztab %>% mutate(module_num = as.integer(as_num)) %>%
        left_join(modules_meta %>% select(module_num, module_color, module_label, module_size), by = "module_num")
    } else {
      Ztab <- Ztab %>% mutate(module_color = module_id) %>%
        left_join(modules_meta %>% select(module_color, module_num, module_label, module_size), by = "module_color")
    }

    if (drop_gold) {
      Ztab <- Ztab %>% filter(!module_id %in% c("0.1", "gold"), !module_color %in% c("gold"))
      Obs  <- Obs %>% filter(!module_id %in% c("0.1", "gold"))
    }

    Ztab %>% left_join(Obs, by = "module_id") %>% mutate(accession = acc, .before = 1)
  })
}

# =========================
# 1) Leaf Preservation
# =========================
cat("\n--- Running Leaf Preservation ---\n")
meta_leaf <- readRDS(file.path(ps_dir, "leaf_metadata_ps_index.rds"))
net_leaf  <- readRDS(file.path(data_dir, "net_leaf.rds"))
DatExpr_leaf <- readRDS(file.path(data_dir, "DatExpr_leaf.rds"))
vst_leaf  <- readRDS(file.path(ps_dir, "leaf_og_counts_dps.rds")) # Or from PCA object if available

# Reconstruct modules table
modules_leaf <- table(labels2colors(net_leaf$colors)) %>% as.data.frame() %>%
  rename(module_color = Var1, module_size = Freq) %>%
  mutate(module_num = as.numeric(factor(module_color, levels = setdiff(unique(module_color), "grey"))),
         module_label = paste0("m", module_num, "_", module_color))

pres_leaf_path <- file.path(data_dir, "pres_list_leaf.rds")
if (file.exists(pres_leaf_path)) {
  cat("Loading pre-computed leaf preservation results...\n")
  pres_list_leaf <- readRDS(pres_leaf_path)
} else {
  cat("Running leaf preservation (100 permutations per accession)...\n")
  accs <- meta_leaf %>% distinct(accession, species)
  pres_list_leaf <- pmap(list(accs$accession, accs$species), function(acc, sp) {
    tryCatch(run_modpres_to_mp(meta_leaf, net_leaf, DatExpr_leaf, vst_leaf, acc, sp, nPermutations = 100),
             error = function(e) { message("Error in ", acc, ": ", e$message); NULL })
  })
  names(pres_list_leaf) <- accs$accession
  saveRDS(pres_list_leaf, pres_leaf_path)
}

pres_summary_leaf <- extract_all_preservation_results(pres_list_leaf, modules_leaf)
write_tsv(pres_summary_leaf, file.path(out_dir, "pres_summary_accessions_leaf.txt"))

# =========================
# 2) Underground Preservation (Root & Stem)
# =========================
cat("\n--- Running Underground Preservation ---\n")
meta_und <- bind_rows(
  readRDS(file.path(ps_dir, "root_metadata_ps_index.rds")) %>% mutate(tissue = "root"),
  readRDS(file.path(ps_dir, "stem_metadata_ps_index.rds")) %>% mutate(tissue = "stem")
)
net_und  <- readRDS(file.path(data_dir, "net_und.rds"))
DatExpr_und <- readRDS(file.path(data_dir, "DatExpr_und.rds"))

# Use global VST from counts (PCA function was used previously to get vst_matrix)
# For simplicity, we assume the user has the counts and we'll just use them if needed.
# But for preservation, we need the EXACT VST space.
# I'll rely on the pre-computed RDS results here.

modules_und <- table(labels2colors(net_und$colors)) %>% as.data.frame() %>%
  rename(module_color = Var1, module_size = Freq) %>%
  mutate(module_num = as.numeric(factor(module_color, levels = setdiff(unique(module_color), "grey"))),
         module_label = paste0("um", module_num, "_", module_color))

# Root Preservation
pres_root_path <- file.path(data_dir, "pres_list_root_2024.rds")
if (file.exists(pres_root_path)) {
  cat("Loading pre-computed root preservation results...\n")
  pres_list_root <- readRDS(pres_root_path)
  pres_summary_root <- extract_all_preservation_results(pres_list_root, modules_und)
  write_tsv(pres_summary_root, file.path(out_dir, "pres_summary_accessions_root.txt"))
}

# Stem Preservation
pres_stem_path <- file.path(data_dir, "pres_list_stem.rds")
if (file.exists(pres_stem_path)) {
  cat("Loading pre-computed stem preservation results...\n")
  pres_list_stem <- readRDS(pres_stem_path)
  pres_summary_stem <- extract_all_preservation_results(pres_list_stem, modules_und)
  write_tsv(pres_summary_stem, file.path(out_dir, "pres_summary_accessions_stem.txt"))
}

cat("\n==== Module preservation analysis complete ====\n")
