# grass-senescence
Scripts and methods for the paper: Dysregulation of underground sink networks limits nitrogen recycling in annual grasses.

## Repository Structure

The project repository is organized into modular directories reflecting the bioinformatic and analytical pipeline:

- `01_mapping_and_counting/`: Shell and R scripts for quality control, adapter trimming, STAR alignment, and gene/orthogroup quantification.
- `02_orthology_inference/`: Scripts for orthogroup identification and cross-species count normalization.
- `03_cross-species_transcriptome_analysis/`: Integrated analysis of the 14-species atlas, including PCA, Z-score expression modeling, and pseudotime trajectory inference.
- `04_photosynthetic_index/`: Calculation of the ps_index and its correlation with physiological nitrogen and amino acid data.
- `05_gene_network_analysis/`: WGCNA and regulatory network inference focused on remobilizing perennial accessions.

## Getting Started

Each directory contains its own `README.md` with specific instructions on how to run the pipeline for that objective.

### Prerequisites

The following software is generally required across the modules:
- STAR (Alignment)
- Trimmomatic (Trimming)
- minimap2 (rRNA filtering)
- R (with Tidyverse and other specialized packages)
- Samtools

### Reproducibility

Scripts are refactored to use relative paths. Ensure you set the appropriate environment variables for software binary paths (e.g., `STAR_BIN`, `TRIMMOMATIC_JAR`) if they are not in your system's `PATH`.
