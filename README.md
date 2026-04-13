# A multi-organ atlas of gene expression during grass senescence 
This repository contains the analysis pipeline and scripts for the study: Dysregulation of underground sink networks limits nitrogen recycling in annual grasses.

## Repository Structure

The repository is organized into directories following the data processing and analysis pipeline described in the paper:

- `01_mapping_and_counting/`: Shell and R scripts for quality control, adapter trimming, STAR alignment, and gene/orthogroup quantification.
- `02_orthology_inference_and_read_aggregation/`: Scripts for orthogroup identification and cross-species count aggregation.
- `03_cross-species_transcriptome_analysis/`: Integrated analysis of the 14-species atlas, including PCA, Z-score expression modeling, and pseudotime trajectory inference.
- `04_photosynthetic_index/`: Derivation of the physiological aging index (`ps_index`).
- `05_N_quant/`: Integration of the ps_index modeling against end-of-season nitrogen combustion tests, free amino acid quantification mappings (UPLC-MS/MS), and independent validations using base-to-tip maize leaf gradients.
- `06_gene_network_analysis/`: WGCNA and regulatory network inference focused on remobilizing perennial accessions.

## Getting Started

Each directory contains its own `README.md` with specific instructions on how to run the pipeline for that objective.

### Prerequisites

Required software is listed in each directory. 

### Reproducibility

Scripts are refactored to use relative paths. Ensure you set the appropriate environment variables for software binary paths (e.g., `STAR_BIN`, `TRIMMOMATIC_JAR`) if they are not in your system's `PATH`.
