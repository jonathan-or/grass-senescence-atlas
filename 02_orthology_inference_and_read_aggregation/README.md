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
- **Purpose**: Imports raw metadata and gene-level counts from the mapping pipeline. Standardizes species names and cleans transcript IDs. Supports both 2023 (pot-in-pot) and 2024 (remobilization) batches.

### 4. `02_aggregate_counts_to_orthogroups.R`
- **Purpose**: Filters the OrthoFinder output and aggregates counts from the gene level to the orthogroup level based on the imported count data.

## Data

- `Msinensis_..._helixer.gff.gz`: Custom Helixer annotation for *Miscanthus sinensis*.
- `Pvirgatum..._helixer.gff.gz`: Custom Helixer annotation for *Panicum virgatum*.
