# Gene Network Analysis (`06_gene_network_analysis`)

This module focuses on identifying co-expression networks and regulatory hubs during grass senescence, with a specific focus on remobilizing perennial accessions.

## Scripts

1. **`01_wgcna_leaf.R`**: Performs WGCNA on leaf samples from the remobilizer cohort.
2. **`02_wgcna_underground.R`**: Performs WGCNA on combined root and stem samples from the remobilizer cohort.
3. **`03_wgcna_module_preservation.R`**: Assesses module preservation (Zsummary) across species and tissues.
4. **`04_wgcna_preservation_plotting.R`**: Generates publication-ready Zsummary plots for annual species.
5. **`05_wgcna_kME_analysis.R`**: Analyzes intramodular connectivity (kME) preservation and produces integrated hub heatmaps.

## Methodology

- **Data Pre-processing**: Expression data is VST-normalized. Photosynthesis-associated orthogroups (`photo_ogs`) are removed to avoid circularity with the `ps_index` pseudotime modeling.
- **Network Construction**: Signed co-expression networks are built using the `WGCNA` package with the `bicor` correlation function.
- **Hub Identification**: Intramodular connectivity (`kME`) is used to identify hub genes. Curated symbols are applied to key senescence hubs.
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
- `Supplemental_Data_kME_Master_Table_v2.csv`: Master table of kME values across all species and tissues.
- `OG_symbols_m2_blue.txt`, `OG_symbols_underground.txt`, `OG_symbols_hubs.txt`: Curated hub gene symbols.

## Requirements

- R packages: `WGCNA`, `tidyverse`, `matrixStats`, `ggrepel`, `uwot`, `patchwork`, `ggplot2`, `clusterProfiler`.
