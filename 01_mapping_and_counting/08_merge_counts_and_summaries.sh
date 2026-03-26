#!/bin/bash

# Set robust base directory (project root)
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
BASE_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)

# Define output files
counts_output="${BASE_DIR}/counts_all_species.txt"
summary_output="${BASE_DIR}/summary_statistics_all_species.txt"

# Check if metadata.txt exists in the root
if [[ ! -f "${BASE_DIR}/metadata.txt" ]]; then
    echo "ERROR: metadata.txt not found at ${BASE_DIR}/metadata.txt! Please provide the correct file."
    exit 1
fi

# Extract sample_id from metadata (excluding header)
awk 'NR>1 {print $1}' "${BASE_DIR}/metadata.txt" > "${BASE_DIR}/metadata_samples.txt"

### **Fix: Ensure Summary Statistics File Keeps the Header**
echo "Merging summary statistics..."
summary_files=( "${BASE_DIR}"/*/summary_statistics_*.txt )

if [[ ${#summary_files[@]} -gt 0 ]]; then
    # Extract the header from the first summary statistics file
    head -n 1 "${summary_files[0]}" > $summary_output

    # Append data while skipping headers from subsequent files
    for file in "${summary_files[@]}"; do
        tail -n +2 "$file" >> temp_summary.txt
    done

    # Filter to include only samples present in metadata.txt
    awk 'NR==1 || FNR==NR {a[$1]; next} $1 in a' "${BASE_DIR}/metadata_samples.txt" temp_summary.txt >> $summary_output

    rm -f temp_summary.txt
fi
echo "Summary statistics merged into $summary_output"

### **Merge Counts Files (Keeping Only Sample IDs in Header)**
echo "Merging counts files..."
count_files=( "${BASE_DIR}"/*/*_counts.txt )

if [[ ${#count_files[@]} -gt 0 ]]; then
    # Extract the first column (GeneID) from the first file
    cut -f1 "${count_files[0]}" > gene_list.txt

    # Initialize counts file with gene IDs
    cp gene_list.txt $counts_output

    for file in "${count_files[@]}"; do
        species=$(basename "$file" | sed 's/_counts.txt//')
        echo "Processing $species..."

        # Extract only the counts (without gene column)
        cut -f2- "$file" > temp_counts.txt

        # Extract only sample IDs from the header
        header=$(head -n 1 "$file" | cut -f2-)

        echo -e "GeneID\t${header}" > temp_header.txt

        # Merge by GeneID
        paste $counts_output temp_counts.txt > temp_merged.txt
        mv temp_merged.txt $counts_output
    done

    # Add proper header to the final file
    cat temp_header.txt $counts_output > temp_final_counts.txt
    mv temp_final_counts.txt $counts_output
fi

# Cleanup
rm -f "${BASE_DIR}/metadata_samples.txt" "${BASE_DIR}/gene_list.txt" temp_counts.txt temp_header.txt

echo "Counts merged into $counts_output"
