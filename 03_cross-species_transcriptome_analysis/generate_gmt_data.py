#!/usr/bin/env python3
# generate_gmt_data.py
# Synthesize GMT files to capture term identities and descriptions across multiple species.

import os
import glob
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description="Generate consolidated GO descriptions from GMT files.")
    parser.add_argument("--go_dir", default="raw_annotations", help="Directory containing species output GMT files.")
    args = parser.parse_args()

    go_dir = args.go_dir
    gmt_files = glob.glob(os.path.join(go_dir, "output_*.gmt"))
    
    if not gmt_files:
        print(f"No GMT files found in {go_dir}. Execute 04a_develop_GO_annotations.R first.")
        return

    species_list = set(f.split("_")[-1].replace(".gmt", "") for f in gmt_files)
    all_data = []

    for species in species_list:
        print(f"Processing {species}...")
        species_gmt_files = glob.glob(os.path.join(go_dir, f"output_*_{species}.gmt"))

        for file in species_gmt_files:
            with open(file, "r") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 3:
                        goid, desc, genes = parts[0], parts[1], parts[2:]
                        all_data.append((goid, desc, "\t".join(genes)))

    if not all_data:
        print("No valid data parsed.")
        return

    # Convert to DataFrame
    df = pd.DataFrame(all_data, columns=["goid", "desc", "genes"])

    # Split gene strings back to lists, uniquely unify them, and join
    df["genes"] = df["genes"].apply(lambda x: x.split("\t"))
    df = df.groupby(["goid", "desc"], as_index=False).agg({"genes": lambda x: "\t".join(sorted(set(g for sublist in x for g in sublist)))})

    output_file = os.path.join(go_dir, "unified_gmt_data.txt")
    df.to_csv(output_file, sep="\t", index=False)
    print(f"Generated unified GMT descriptor file: {output_file}")
    
    # Optional mapping script generates gmt_data.rds in R,
    # Here we are just exporting the final text artifact.

if __name__ == "__main__":
    main()
