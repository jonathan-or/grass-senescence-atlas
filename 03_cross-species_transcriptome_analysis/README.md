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
