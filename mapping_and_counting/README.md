# Mapping and Counting Methods

This directory contains the scripts used for read processing, alignment, and generation of gene counts. The full pipeline ensures high-quality read alignments and reproducible quantification of gene expression.

## 1. Data Preparation
Raw reads are initially organized by renaming them to match metadata identifiers (`01_rename_raw_files.sh`) and then structured according to the experimental design (`02_move_raw_reads.sh`). A metadata cleaning step (`00_clean_metadata.sh`) ensures uniformity across samples before downstream analysis.

## 2. Quality Control and Trimming
Raw sequence data quality is evaluated, and adapter sequences along with low-quality ends are removed using **Trimmomatic** (`03_trimming_and_QC.sh`). 
- **Adapters**: Trimmed utilizing the `TruSeq3-SE.fa` adapter file.
- **Parameters**: `ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`.
Following trimming, read quality is verified using **FastQC**. Summary statistics of the trimming step are subsequently computed (`03_2_trimming_stats.sh`).

## 3. rRNA Filtering
To remove ribosomal RNA contamination, clean reads are mapped against an aggregated index of Arabidopsis and Zea mays rRNA sequences (`Arabidopsis_Zea_rRNA.fasta`). This step (`04_rRNA_filtering.sh`) utilizes [tool name, e.g., Bowtie2] to safely discard rRNA reads while keeping high-quality mRNA sequences for alignment.

## 4. Genome Alignment
Filtered reads are aligned to the respective reference genomes using **STAR** (Spliced Transcripts Alignment to a Reference, `05_STAR_alignment.sh`).
- A STAR genome index is first dynamically generated if not already available for the species.
- Alignments are performed on 32/60 threads using basic two-pass mode, and alignments are output as coordinate-sorted BAM files (`--outSAMtype BAM SortedByCoordinate`).

## 5. Read Counting and Summarization
Transcript quantification is carried out directly after alignment via **featureCounts** (from the Subread package, handled by `06_featureCounts.R` scripts). 
Finally, individual sample counts and alignment statistics are consolidated into comprehensive outcome tables (`07_generate_summary_statistics.sh` and `08_merge_counts_and_summaries.sh`) to facilitate downstream statistical modeling and differential expression analysis.
