#!/bin/bash

# Set robust base directory (project root)
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
BASE_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)

# Standard metadata file for the pipeline
cleaned_metadata="${BASE_DIR}/metadata.txt"

if [ ! -f "$original_metadata" ]; then
    echo "Metadata source $original_metadata not found. Checking for existing cleaned metadata..."
    if [ -f "$cleaned_metadata" ]; then
        echo "Found existing $cleaned_metadata. No cleaning needed."
    else
        echo "ERROR: Neither $original_metadata nor $cleaned_metadata found!"
        exit 1
    fi
else
    # Remove carriage return characters (\r) and create metadata.txt
    tr -d '\r' < "$original_metadata" > "$cleaned_metadata"
    echo "Carriage return characters removed. Created: $cleaned_metadata"
fi
