#!/bin/bash

# Set robust base directory (project root)
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
BASE_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)

# Software binaries
GFFREAD_BIN=${GFFREAD_BIN:-"gffread"}

# Directories for input and output
# Default to current directory if not provided
DIRECTORY=${1:-"${SCRIPT_DIR}"}
OUTPUT_DIR=${2:-"${SCRIPT_DIR}/faas"}

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Loop through each .gff file (handling potential .gz as well)
for gff_file in "$DIRECTORY"/*.gff "$DIRECTORY"/*.gff3 "$DIRECTORY"/*.gff.gz "$DIRECTORY"/*.gff3.gz; do
    [ -e "$gff_file" ] || continue
    echo "Processing $gff_file..."

    # Extract the base name
    baseName=$(basename "$gff_file")
    baseName=${baseName%.gz}
    baseName=${baseName%.gff3}
    baseName=${baseName%.gff}

    # Construct the species name or use the base name
    speciesName=$baseName

    # Define the FASTA file path (assuming genome fasta is in the same dir)
    fasta="$DIRECTORY/${speciesName}.fasta"
    if [ ! -f "$fasta" ]; then
        fasta="$DIRECTORY/${speciesName}.fna"
    fi

    # Define the output path for gffread
    gffread_output="$OUTPUT_DIR/${baseName}.faa"

    # Run gffread command (-y for protein extraction)
    # If file is zipped, we might need a temporary unzipped version or use a pipe if gffread supports it
    # For now, assume unzipped for simplicity or add a check
    if [[ "$gff_file" == *.gz ]]; then
        echo "Decompressing $gff_file temporarily..."
        gunzip -c "$gff_file" > temp.gff
        "$GFFREAD_BIN" -y "$gffread_output" -g "$fasta" temp.gff
        rm temp.gff
    else
        "$GFFREAD_BIN" -y "$gffread_output" -g "$fasta" "$gff_file"
    fi

    echo "Protein sequences saved to $gffread_output"
done

echo "All files processed."
