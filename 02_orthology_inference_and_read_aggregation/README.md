# Objective 02: Orthology inference and read aggregation

This directory contains scripts for identifying orthogroups across the studied species and aggregating gene-level transcript counts into orthogroup-level counts.

## Methodology

To establish a comparative framework across the 14 studied species, we performed orthogroup inference using OrthoFinder (v2.5.4). Input protein sequences (.faa) were extracted from GFF3 and FASTA files using `gffread` (v0.10.4). For *Zea mays*, a custom Python script filtered the GFF3 file (B73 v5) for primary transcripts (canonical_transcript=1) before protein extraction. OrthoFinder was executed with parameters `-S diamond -M msa -t 36 -a 36`. 

Following inference, the `Orthogroups.tsv` output was filtered to retain a core set of orthogroups present across all 14 species, restricted to those containing ≤75 total gene members to minimize noise. Raw counts from `featureCounts` were then aggregated at the orthogroup level by summing the counts of all constituent genes within each bin.

## Scripts and Files

### 1. `get_faas.sh`
- **Purpose**: Automates protein sequence extraction from GFF/FASTA pairs using `gffread`.

### 2. `primary_transcript_from_gff.py`
- **Purpose**: Filters *Zea mays* GFF3 for primary/canonical transcripts.

### 3. `01_import_counts_and_metadata.R`
- **Purpose**: Illustrates the processing flow to import raw metadata and gene-level counts from the mapping pipeline. Standardizes species names and cleans transcript IDs.
- **Expected Inputs in Base Directory**:
  - `metadata.txt`: Experimental metadata with a `sample_id` column.
  - Subdirectories for each species containing:
    - `{species}_counts.txt` (or `_merged_counts.txt`): Transcript-level quantifications.
    - `{species}_summary_statistics.txt`: Output from the mapping pipeline.

### 4. `02_aggregate_counts_to_orthogroups.R`
- **Purpose**: Filters the OrthoFinder output (`Orthogroups.tsv.gz`) and aggregates counts from the gene level to the orthogroup level based on the imported count data.
- **Expected Inputs in Base Directory**:
  - `Orthogroups.tsv.gz`: Output from OrthoFinder.
  - `02_orthology_inference_and_read_aggregation/count_data_frames.rds`: Output from script `01`.

### 5. `03_plot_orthogroup_summary.R`
- **Purpose**: Supplementary Figure 3 (panels A-D; the species tree and schematic in A were drawn separately): genes per species in core orthogroups, in any orthogroup and annotated in total (A, bars; also written to `plots/orthogroup_gene_sets_per_species.txt`); number of species per orthogroup (B); genes per orthogroup for all orthogroups and the core set (C); share of each library's reads in core orthogroups (D).
- **Inputs**: `Orthogroups.tsv.gz`, `total_gene_counts_named.tsv` (annotated genes per species) and, for panel D, the per-species gene counts from script `01` plus `03_cross-species_transcriptome_analysis/data/og_counts.rds`. Panel D is skipped if the gene counts are not available.
- **Usage**: `Rscript 02_orthology_inference_and_read_aggregation/03_plot_orthogroup_summary.R [counts_rds] [out_dir]` (defaults: `count_data_frames.rds` in this directory and `plots/`).
- **Note**: the gene counts in the published panel A bars (e.g. 38,727 for *P. virgatum*) are 6-11% higher than the genes in the 16,778 core orthogroups used in the analysis (36,613 for *P. virgatum*), which is what this script plots; panels B-D match the published figure.

## Data Structure for Reproducibility

To successfully run scripts `01` and `02`, a subset or full dataset must be structured as follows:

```
[base_dir]/
├── metadata.txt
├── Orthogroups.tsv.gz
├── Zea_mays/
│   ├── Zea_mays_counts.txt
│   └── Zea_mays_summary_statistics.txt
├── Panicum_virgatum/
│   ├── Panicum_virgatum_counts.txt
│   └── Panicum_virgatum_summary_statistics.txt
└── ... (other species)
```

**Note:** For downstream cross-species analysis (Objective 03), pre-processed and normalized 'final' counts and metadata objects will be provided, circumventing the need to manually run these import and aggregation steps unless actively modifying the fundamental mapping pipeline.

## Data

- `Msinensis_..._helixer.gff.gz`: Helixer annotation for *Miscanthus sinensis*.
- `Pvirgatum..._helixer.gff.gz`: Helixer annotation for *Panicum virgatum*.
- `Sorghum_bicolor..._helixer.gff.gz`: Helixer annotation for *Sorghum bicolor* (BTx623).
