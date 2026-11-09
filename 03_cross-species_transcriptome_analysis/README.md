# Objective 03: Cross-species Transcriptome Analysis

This directory contains the scripts and structured data to reproduce the cross-species integration, dimensionality reduction (PCA and UMAP), and downstream Gene Ontology (GO) functional enrichment. 

## Directory Structure

### Data & Annotations

- **`data/`**: To enable fast reproducibility without running bulky parsing steps, the final core variables for the cross-season analysis are cached here as `.rds` objects (`metadata.rds`, `og_counts.rds`, `filteredOG.rds`). 
- **`outlier_libraries.txt`**: A comprehensive record of all sequencing libraries that were dropped during the initial QC steps due to poor mapping, low reads, or being severe PCA outliers.
- **`raw_annotations.tar.xz`**: Contains the raw Gene Ontology (GO) assignment files produced by the PANNZER web server (`*_GO.out`) as well as the `.gmt` intermediate flat files. It is heavily compressed to meet GitHub sizes limit (extract via `tar -xf raw_annotations.tar.xz`).

### Processing Code & Analyses

All analysis stages are parameterized and modularized.

1. **`01_QC_and_PCA.R`**
   Script used to process each individual season/batch, executing variance-stabilizing transformations and library composition outlier filtering.

2. **`02_consolidate_seasons.R`**
   Blends the filtered outputs of 2023 and 2024 seasons into the comprehensive core orthogroup expression matrices (`og_counts.rds`) and aligned `metadata.rds`, merging external environmental weather station data. 

3. **`Functions.R`**
   A master reference script containing the generalized visualization operations (like `plot_pca_ggplot2` and `run_umap_seurat_style`) and testing paradigms (like `perform_GO_enrichment_OG`). Sourced by all subsequent downstream scripts to ensure methods remain strictly uniform.

4. **`03_plot_PCA_and_UMAP.R`**
   Extracts the core normalized counts and produces the publication-style UMAP and PCA clusters. It demonstrates the initial identification of 3,000 highly-variable orthogroups before clustering.

### Functional Annotation Workflow 

For complete transparency, the methodology used to translate PANNZER outputs into testing models is included below. 

5. **`04a_develop_GO_annotations.R`**
   R script that ingests the raw PANNZER outputs from `raw_annotations.tar.xz` and parses them into Gene Matrix Transposed (`.gmt`) lists per species.
   
6. **`generate_gmt_data.py`**
   Python helper script that merges parallel species `.gmt` definitions into a unified terminology definition list.

7. **`04b_perform_GO_enrichment.R`**
   The final analytical pipeline. It reads the consolidated `.rds` annotation dictionaries from `data/` and uses hypergeometric testing to compute and plot significant functional enrichment scores via `clusterProfiler` for features of interest.

### Leaf pseudotime

8. **`05_pseudotime_analysis.R`**
   PhenoPath pseudotime of the 1,146 leaf libraries (Figure 2B-E, Supplementary Figure 7). The fitted model is cached in `data/leaf_pseudotime_model.rds`. Also writes the 115 photosynthetic orthogroups (`data/photo_ogs.rds`, `data/photo_ogs_Zm.txt`) used by `04_photosynthetic_index`.

### Library QC and atlas composition (Supplementary Figures 4-6)

The atlas uses the 2,685 libraries that passed QC (mapping rate > 35% and > 100,000 mapped reads; `data/og_counts.rds`, `data/metadata.rds`). The QC figures also need all 3,798 sequenced libraries before QC.

9. **`06_merge_all_libraries.R`**
   Builds the unfiltered dataset of all 3,798 libraries from the per-season outputs of modules 01-02: `data/metadata_all_libraries.rds` and `data/og_counts_all_libraries.rds`. Usage: `Rscript 03_cross-species_transcriptome_analysis/06_merge_all_libraries.R <dir_2023> <dir_2024>`.
10. **`07_qc_all_libraries.R`**
    Supplementary Figure 4: PCA of all libraries colored by tissue, week, species, season (A-D) and by mapping rate, mapped reads and the QC cutoffs (E-I); passed/failed counts by species, tissue, season and week (J). The VST and PCA of all libraries are cached in `data/` (`vst_all_libraries.rds`, `pca_all_libraries.rds`).
11. **`08_qc_sensitivity.R`**
    Supplementary Figure 5: PCA at increasing mapping-rate and read-count thresholds. One library (`AC24BRTL002104`, no mapped reads) has no counts and is dropped from the VST, so the unfiltered PCA has 3,797 libraries.
12. **`09_library_composition_and_season_pca.R`**
    Supplementary Figure 6: libraries per species, tissue and season (A-C) and leaf and root PCA colored by season, PC1 vs PC2-PC4 (D-E).

Scripts take optional `[in_dir] [out_dir]` arguments (defaults: `data/` and `plots/` of this module) and are run from the repository root.
