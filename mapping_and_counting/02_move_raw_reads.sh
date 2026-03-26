#!/bin/bash

# Set robust base directory (project root)
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
BASE_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)

# Define the source directories and destination directory
src_dir_09="${BASE_DIR}/raw_reads/AC24BRTL07/demultiplexed"
src_dir_10="${BASE_DIR}/raw_reads/AC24BRTL08/demultiplexed"
src_dir_11="${BASE_DIR}/raw_reads/MC24BRTL02/demultiplexed"
dest_dir="${BASE_DIR}/raw_reads"

# Move only renamed .fastq.gz files (those that do not match plate_pos pattern)
find "$src_dir_09" -type f -name "*.fastq.gz" ! -name '[A-Z][0-9][0-9].fastq.gz' -exec mv {} "$dest_dir" \;
find "$src_dir_10" -type f -name "*.fastq.gz" ! -name '[A-Z][0-9][0-9].fastq.gz' -exec mv {} "$dest_dir" \;
find "$src_dir_11" -type f -name "*.fastq.gz" ! -name '[A-Z][0-9][0-9].fastq.gz' -exec mv {} "$dest_dir" \;

echo "Move completed, plate_pos files were skipped."
