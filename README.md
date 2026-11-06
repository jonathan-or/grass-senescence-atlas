# A transcriptomic atlas of grass senescence reveals divergent underground sink networks promote nitrogen recycling in perennials

Analysis pipeline and scripts for the study:

> Ojeda-Rivera JO, Oren E, Hsu S-K, Lepak N, La T, Zhai J, Stitzer MC, Yobi A, Angelovici R, Buckler ES, Romay MC. *A transcriptomic atlas of grass senescence reveals divergent underground sink networks promote nitrogen recycling in perennials.* Accepted; in production (full citation and DOI to be added on publication).

The atlas comprises 2,685 RNA-seq (BRB-seq) libraries from field-grown plants of 14 grass species in the Panicoideae (Poaceae), with maize and sorghum as annual references. It spans leaves, roots, stalks, and rhizomes across two seasons (2023 and 2024), from mid-growing season to senescence. A photosynthetic index (`ps_index`) aligns the leaf's transition from nitrogen sink to source across species, and co-expression network analysis compares leaf and underground nitrogen-recycling subnetworks between perennials and annuals.

## Data availability

- **Raw reads**: NCBI BioProject [PRJNA1458764](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1458764).
- **Interactive expression atlas**: <https://maize-genetics.github.io/atlas-of-grass-senescence-dashboard/>
- **Supplementary Datasets SD1–SD16** (library metadata, core orthogroups, network modules, module preservation, comparative kME, etc.) are provided with the article.
- **Files referenced in the Methods and hosted here**:
  - rRNA reference used for filtering: `01_mapping_and_counting/Arabidopsis_Zea_rRNA.fasta`
  - Helixer gene annotations for *P. virgatum* (WBC HAP1 v1.0), *M. sinensis* (v7.0), and *S. bicolor* (BTx623): `02_orthology_inference_and_read_aggregation/*_helixer.gff.gz`
  - OrthoFinder orthogroups: `02_orthology_inference_and_read_aggregation/Orthogroups.tsv.gz`

## Repository Structure

The repository is organized into directories following the data processing and analysis pipeline described in the paper:

- `01_mapping_and_counting/`: Shell and R scripts for adapter trimming and QC, rRNA filtering, STAR alignment, and gene-level quantification with featureCounts.
- `02_orthology_inference_and_read_aggregation/`: Scripts for orthogroup inference and aggregation of gene counts into core orthogroup counts.
- `03_cross-species_transcriptome_analysis/`: Integrated analysis of the 14-species atlas, including PCA, UMAP, pseudotime trajectory inference, and GO enrichment.
- `04_photosynthetic_index/`: Derivation of the photosynthetic index (`ps_index`) and removal of photosynthesis-associated orthogroups before network analysis.
- `05_N_quant/`: Integration of the `ps_index` with end-of-season nitrogen combustion tests and free amino acid quantification (UPLC-MS/MS), and validation on an independent base-to-tip maize leaf gradient.
- `06_gene_network_analysis/`: WGCNA co-expression networks (leaf and underground) of N-remobilizing perennials, module preservation in annuals and other perennials, and comparative kME.

## Getting Started

Each directory contains its own `README.md` with specific instructions on how to run the pipeline for that objective.

### Prerequisites

Required software is listed in each directory. Versions used in the study (as reported in the Methods):

| Step | Software |
|---|---|
| Demultiplexing | BRB-seqTools v1.6 |
| Trimming / QC | Trimmomatic v0.36, FastQC |
| rRNA filtering | Minimap2 v2.27 (`-ax sr`), Samtools |
| Alignment | STAR v2.7.10b |
| Gene counts | Rsubread `featureCounts` (R 4.2) |
| Gene annotation (*P. virgatum*, *M. sinensis*, *S. bicolor*) | Helixer v0.3.2 (`land_plant`) |
| Orthology | OrthoFinder v2.5.4, DIAMOND v2.0.15, gffread v0.10.4 |
| Normalization (VST) | DESeq2 v1.46 |
| Pseudotime | PhenoPath v1.30 |
| GO enrichment | clusterProfiler v4.14.4 |
| Co-expression networks | R 4.4.1, WGCNA v1.73, DESeq2 v1.46.0 |

### Reproducibility

Scripts are refactored to use relative paths and are run from the repository root. Ensure you set the appropriate environment variables for software binary paths (e.g., `STAR_BIN`, `TRIMMOMATIC_JAR`) if they are not in your system's `PATH`.

Computationally expensive intermediate objects (e.g., WGCNA networks, module preservation results) are provided as pre-computed files so that downstream figures can be regenerated directly; module READMEs describe how to rebuild them.

## License

Code in this repository is released under the [MIT License](LICENSE).

## Contact

Jonathan O. Ojeda-Rivera (joo29@cornell.edu), Institute for Genomic Diversity, Cornell University.
