# Gene Network Analysis (`06_gene_network_analysis`)

This module focuses on identifying co-expression networks and regulatory hubs during grass senescence, with a specific focus on remobilizing perennial accessions.

## Scripts

1. **`01_wgcna_leaf.R`**: Performs WGCNA on leaf samples from the remobilizer cohort.
2. **`02_wgcna_underground.R`**: Performs WGCNA on combined root and stem samples from the remobilizer cohort.
3. **`03_wgcna_module_preservation.R`**: Assesses module preservation (Zsummary) across species and tissues. Test expression is the VST used to build each reference network: the leaf VST for leaf, the combined root + stem VST for root (2024 season, plus *S. bicolor* from 2023) and stem (2024 season).
4. **`04_wgcna_preservation_plotting.R`**: Generates publication-ready Zsummary plots for annual species.
5. **`05_wgcna_kME_analysis.R`**: Builds the per-accession kME master table (Supplementary Data SD13), analyzes intramodular connectivity (kME) preservation and produces integrated hub heatmaps.

## Methodology

- **Data Pre-processing**: Expression data is VST-normalized. Photosynthesis-associated orthogroups (`photo_ogs`) are removed to avoid circularity with the `ps_index` pseudotime modeling.
- **Network Construction**: Signed co-expression networks are built using the `WGCNA` package with the `bicor` correlation function.
- **Hub Identification**: Intramodular connectivity (`kME`) is used to identify hub genes. Curated symbols are applied to key senescence hubs.
- **Per-accession kME**: For each accession, module eigengenes are recomputed from that accession's libraries only, and kME is the Pearson correlation of each OG with its own module's eigengene. Leaf kME uses the leaf VST; underground kME uses a root-only VST (all root libraries). OGs with zero variance within an accession (no reads in any of its libraries) are left out of that accession's eigengene and their kME is reported as `NA`.
- **Visualization**:
    - **TOM-UMAP**: A 2D embedding of the Topological Overlap Matrix (TOM) is used to visualize network architecture and module relationships.
    - **Dynamics Plots**: Median z-score expression profiles are projected onto the TOM-UMAP to visualize network-wide senescence progression.
- **Enrichment**: Gene Ontology (GO) enrichment is performed for all detected modules to characterize biological functions.

## Data Reproducibility

To bypass high-memory Topological Overlap Matrix (TOM) computations and ensure consistent UMAP layouts, pre-computed objects are provided in `data/`:

- `net_leaf.rds`, `net_und.rds`: WGCNA network objects.
- `DatExpr_leaf.rds`, `DatExpr_und.rds`: Filtered expression matrices.
- `09_TOM_UMAP_embedding_leaf.csv`, `09_TOM_UMAP_embedding_underground.csv`: Pre-computed UMAP coordinates.
- `pres_list_leaf.rds`, `pres_list_root_2024.rds`, `pres_list_stem.rds`: Pre-computed module preservation results (100 permutations).
- `preservation_list_und.rds`: Underground preservation across all seasons; only *S. bicolor* (`Sb_BTx623`, sampled in 2023) is used, merged into the root results.
- `Supplemental_Data_kME_Master_Table_v2.csv`: Master table of per-accession kME values (leaf m2_blue; underground um1_turquoise and um2_blue) for the five focal accessions.
- `OG_symbols_m2_blue.txt`, `OG_symbols_underground.txt`, `OG_symbols_hubs.txt`: Curated hub gene symbols.

Scripts load these files when present. To regenerate them, delete the file or set `force_recompute <- TRUE` at the top of `03_wgcna_module_preservation.R` or `05_wgcna_kME_analysis.R`. Recomputed preservation results reproduce the provided ones; recomputed kME values match the provided table except that OGs with zero variance in an accession are `NA` (five maize underground OGs, provided as ~1e-16).

## Requirements

- R packages: `WGCNA`, `DESeq2`, `tidyverse`, `matrixStats`, `ggrepel`, `uwot`, `patchwork`, `ggplot2`, `ggpubr`, `clusterProfiler`.
- Versions used for the published analysis: R 4.4.1, Bioconductor 3.20, `WGCNA` 1.73, `DESeq2` 1.46.0, `dplyr` 1.1.4, `matrixStats` 1.5.0.
