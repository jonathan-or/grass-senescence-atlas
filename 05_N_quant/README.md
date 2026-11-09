# 05_N_quant: Nitrogen Quantification and Amino Acid Dynamics

This module integrates end-of-season nitrogen combustion quantification and free amino acid analysis against phenotypic datasets.

## Objectives
- Track the physiological allocation of raw Nitrogen across structures over consecutive weeks.
- Ingest and bind raw field quantification data into cohesive analytical targets.
- Correlate absolute tissue Nitrogen percentages implicitly against relative transcription-derived phenomena like the Photosynthesis Index.
- Monitor temporal profiles of mass-spectroscopy mapped Free Amino Acids, particularly emphasizing Glutamine transport dynamics compared to senescence-activated indicator genes like `SAG12` and `SAG15`.

## Main Scripts
- **`01_combustion_analysis.R`**: Binds legacy N tables and dynamically structures them against the primary leaf metadata. Outputs global allocation models predicting plant source-sink physiological movements into `plots/`.
- **`02_amino_acid_analysis.R`**: Isolates `Tripsacum dactyloides` timelines against raw amino acid distributions. Exports analytical matrices to `Supplementary_Dataset_Tripsacum_AA.csv` and builds normalized `log2` models showing Glutamine tracking.
- **`03_wang_validation.R`**: Validates the photosynthetic index on the maize base-to-tip leaf gradient (Wang et al. 2014) against protein, amino acids, RuBisCO activity, chlorophyll, glucose and starch (Supplementary Figure 9).
- **`04_n_allocation_all_accessions.R`**: Supplementary Figure 10. Underground / aboveground N allocation at weeks 7 and 14 for all eight accessions, with Tukey HSD per accession (A); *Zea mays* %N concentration ratios including kernels (`data/Zm_kernel_N_quant.txt`), with ANOVA and Tukey HSD (B). Statistics are written next to the plots.
- **`05_maize_gradient_source_sink.R`**: Supplementary Figure 11. Share of the transcriptome from the maize members of the leaf source (m2_blue, m7_black, m3_brown) and sink (m10_purple, m4_yellow, m6_red) orthogroups along the maize leaf gradient (A-B), a continuous DESeq2 model along the gradient (C, written as a table) and z-scores of the source genes that decline towards the tip (D). Gene sets are taken from `06_gene_network_analysis/results/06a_gene_module_assignments_leaf.csv`.

Scripts 04 and 05 take an optional output directory (default `plots/`) and are run from the repository root.

## Pre-requisites
Scripts in this module demand upstream evaluation arrays:
- Global `metadata.rds` mappings.
- Output arrays evaluated by `03_cross-species_transcriptome_analysis`.
- Explicit calculations evaluated via `04_photosynthetic_index`, inherently relying on `leaf_metadata_ps_index.rds`.
