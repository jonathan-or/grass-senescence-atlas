# grass-senescence
Scripts and methods for the paper: Dysregulation of underground sink networks limits nitrogen recycling in annual grasses.

## Repository Structure

The project is organized by research objective to ensure reproducibility and clarity:

- `mapping_and_counting/`: Bash and R scripts for quality control, trimming, alignment (STAR), and gene/orthogroup counting.
- `orthogroup_inference_and_count_binning/`: (Coming soon) Scripts for orthology analysis.
- `dimensionality_reduction/`: (Coming soon) PCA and other dimension reduction analyses.
- `photosynthetic_index/`: (Coming soon) Calculation of photosynthetic indices.
- `gene_network_analysis/`: (Coming soon) WGCNA and other network-based analyses.

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
