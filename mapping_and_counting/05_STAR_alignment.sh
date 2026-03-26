#!/bin/bash

# Check if a species name is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <Species Name>"
    exit 1
fi

species=$1

# Set robust base directory (project root)
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
BASE_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)

# Set STAR path - can be overridden by environment variable
STAR_BIN=${STAR_BIN:-"/programs/STAR-2.7.10b/STAR"}

# Set directories based on the species name
metadata="${BASE_DIR}/metadata.txt"
inDIR="${BASE_DIR}/clean_reads/"
outDIR="${BASE_DIR}/${species}/alignments/"
genomeFA="${BASE_DIR}/${species}/${species}.fasta"
annotation="${BASE_DIR}/${species}/${species}_hxr.gff"
genomeDir="${BASE_DIR}/${species}/genomeIndex/${species}/"
sampleIDs="${BASE_DIR}/${species}_sample_ids.txt"

# Ensure output directories exist
mkdir -p "$outDIR"

# Check if genome files exist before proceeding
if [[ ! -f "$genomeFA" || ! -f "$annotation" ]]; then
    echo "Required genome files missing for $species: $genomeFA or $annotation"
    exit 1
fi

# Generate the genome index if not already present
if [ ! -d "${genomeDir}" ]; then
    mkdir -p "$genomeDir"
    echo "Generating genome index for $species..."
    "$STAR_BIN" \
         --runThreadN 60 \
         --runMode genomeGenerate \
         --genomeDir "${genomeDir}" \
         --genomeFastaFiles "$genomeFA" \
         --sjdbGTFfile "$annotation"
else
    echo "Genome index already exists for $species."
fi
echo "Continuing with further processing..."

# Extract relevant sample IDs based on species from metadata
awk -v sp="$species" -F'\t' '$2 == sp {print $1}' "$metadata" > "$sampleIDs"

# Perform alignments using the extracted sample IDs
while read sample_id; do
    fastqFile="${inDIR}${sample_id}_clean.fastq.gz"
    if [ -f "$fastqFile" ]; then
        if [ ! -f "${outDIR}${sample_id}_Aligned.sortedByCoord.out.bam" ]; then
            echo "Running alignment for $sample_id..."
            "$STAR_BIN" \
                 --outFileNamePrefix "${outDIR}${sample_id}_" \
                 --outSAMtype BAM SortedByCoordinate \
                 --outSAMstrandField intronMotif \
                 --genomeDir "${genomeDir}" \
                 --runThreadN 32 \
                 --readFilesCommand zcat \
                 --readFilesIn "$fastqFile" \
                 --twopassMode Basic
        else
            echo "Alignment already completed for $sample_id, skipping..."
        fi
    else
        echo "FASTQ file not found for $sample_id"
    fi
done < "$sampleIDs"

echo "All alignments for $species completed."

# Call the appropriate R script for counting reads based on species
# Assumes R scripts are in the same directory as this script
if [ "$species" = "Zea_mays" ]; then
    Rscript "${SCRIPT_DIR}/06_featureCounts_Zm.R" "$species"
else
    Rscript "${SCRIPT_DIR}/06_featureCounts.R" "$species"
fi

# After counting, run the summary statistics script
bash "${SCRIPT_DIR}/07_generate_summary_statistics.sh" "$species"
