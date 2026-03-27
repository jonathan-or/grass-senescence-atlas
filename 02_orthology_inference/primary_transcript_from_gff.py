import argparse
import sys

def filter_gff(input_file, output_file):
    with open(input_file, 'r') as file:
        lines = file.readlines()

    relevant_ids = set()
    for line in lines:
        parts = line.strip().split('\t')
        if len(parts) > 2 and parts[2] == 'mRNA':
            attributes = parts[8].split(';')
            for attribute in attributes:
                if 'longest=1' in attribute or 'canonical_transcript=1' in attribute:
                    id = attributes[0].split('=')[1]
                    relevant_ids.add(id)
                    break

    filtered_lines = []
    for line in lines:
        parts = line.strip().split('\t')
        if len(parts) > 8:
            attributes = parts[8].split(';')
            parent_id = None
            for attribute in attributes:
                if attribute.startswith('Parent='):
                    parent_id = attribute.split('=')[1]
                    break

            if parent_id and parent_id in relevant_ids:
                filtered_lines.append(line)
            elif len(parts) > 2 and parts[2] == 'mRNA' and parts[8].split(';')[0].split('=')[1] in relevant_ids:
                filtered_lines.append(line)

    with open(output_file, 'w') as file:
        for line in filtered_lines:
            file.write(line)

# Setup command line argument parsing
parser = argparse.ArgumentParser(description="Filter a GFF file to keep longest transcripts and their children. Usage: python3 primary_transcript_from_gff.py <input.gff> <output.gff>.")
parser.add_argument('input_file', type=str, help='Path to the input GFF file')
parser.add_argument('output_file', type=str, help='Path to the output GFF file')

# Parse arguments
args = parser.parse_args()

# Run the function with provided arguments
filter_gff(args.input_file, args.output_file)
