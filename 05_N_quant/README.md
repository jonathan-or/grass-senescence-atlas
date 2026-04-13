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

## Pre-requisites
Scripts in this module demand upstream evaluation arrays:
- Global `metadata.rds` mappings.
- Output arrays evaluated by `03_cross-species_transcriptome_analysis`.
- Explicit calculations evaluated via `04_photosynthetic_index`, inherently relying on `leaf_metadata_ps_index.rds`.
