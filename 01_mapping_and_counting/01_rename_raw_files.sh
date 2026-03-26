#!/bin/bash

# Set robust base directory (project root)
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
BASE_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)

# Path to the metadata file
metadata="${BASE_DIR}/metadata.txt"

# Define the base directory for raw reads
raw_reads_dir="${BASE_DIR}/raw_reads"

# Read through the metadata file line by line
while IFS=$'\t' read -r -a line; do
    # Skip header line
    if [[ "${line[0]}" == "sample_id" ]]; then
        continue
    fi

    # Extract information from columns
    sample_id="${line[0]}"
    plate_pos="${line[6]}"
    batch="${line[7]}"

    # Define path to the batch directory
    batch_dir="$raw_reads_dir/$batch/demultiplexed"

    # Check if the directory exists
    if [[ -d "$batch_dir" ]]; then
        # Rename files based on plate position and sample ID
        if [[ -f "$batch_dir/${plate_pos}.fastq.gz" ]]; then
            mv "$batch_dir/${plate_pos}.fastq.gz" "$batch_dir/${sample_id}.fastq.gz"
        else
            echo "File not found: $batch_dir/${plate_pos}.fastq.gz"
        fi
    else
        echo "Batch directory not found: $batch_dir"
    fi
done < "$metadata"

echo "Renaming completed."
