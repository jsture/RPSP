#!/usr/bin/env python3
"""
Script to generate YAML files for all combinations of GLUR sequences from glurs_sequences.csv
and SMILES strings from all_mols.txt, using the template.yaml format.
"""

import csv
import os
from pathlib import Path


def read_glur_sequences(csv_file):
    """
    Read GLUR sequences from CSV file.

    Returns:
        list: List of tuples (gene_name, sequence)
    """
    sequences = []
    with open(csv_file, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            sequences.append((row["gene_name"], row["sequence"]))
    return sequences


def read_smiles_data(txt_file):
    """
    Read SMILES data from text file.

    Returns:
        list: List of tuples (smiles_id, smiles_string)
    """
    smiles_data = []
    with open(txt_file, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                # Parse the line: 'ID','SMILES'
                parts = line.split(",")
                if len(parts) == 2:
                    smiles_id = parts[0].strip("'")
                    smiles_string = parts[1].strip("'")
                    smiles_data.append((smiles_id, smiles_string))
    return smiles_data


def generate_yaml_files():
    """
    Generate YAML files for all combinations of GLUR sequences and SMILES.
    """
    # Define input file paths
    sequences_file = "../mols/glurs_sequences.csv"
    smiles_file = "../mols/all_mols.txt"
    output_dir = "../yamls"

    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Read input data
    print("Reading GLUR sequences...")
    glur_sequences = read_glur_sequences(sequences_file)
    print(f"Found {len(glur_sequences)} GLUR sequences")

    print("Reading SMILES data...")
    smiles_data = read_smiles_data(smiles_file)
    print(f"Found {len(smiles_data)} SMILES entries")

    # Generate all combinations
    total_combinations = len(glur_sequences) * len(smiles_data)
    print(f"Generating {total_combinations} YAML files...")

    file_count = 0
    for gene_name, sequence in glur_sequences:
        for smiles_id, smiles_string in smiles_data:
            # Create filename encoding both identifiers
            # Format: {gene_name}_{smiles_id}.yaml
            filename = f"{gene_name}_{smiles_id}.yaml"
            filepath = os.path.join(output_dir, filename)

            # Write YAML file with custom formatting
            with open(filepath, "w") as f:
                # Write version
                f.write("version: 1\n")
                f.write("sequences:\n")

                # Write protein section
                f.write("  - protein:\n")
                f.write("      id: A\n")
                f.write(f"      sequence: {sequence}\n")

                # Write ligand section
                f.write("  - ligand:\n")
                f.write("      id: B\n")
                f.write(f"      smiles: '{smiles_string}'\n")

                f.write("properties:\n")
                f.write("  - affinity:\n")
                f.write("      binder: B\n")

            file_count += 1

            # Print progress every 50 files
            if file_count % 50 == 0:
                print(f"Generated {file_count}/{total_combinations} files...")

    print(f"\nSuccessfully generated {file_count} YAML files in {output_dir}")

    # Print summary
    print("\nSummary:")
    print(f"- GLUR sequences: {len(glur_sequences)}")
    print(f"- SMILES entries: {len(smiles_data)}")
    print(f"- Total combinations: {total_combinations}")
    print(f"- Files generated: {file_count}")

    # Show some example filenames
    print("\nExample filenames:")
    example_files = []
    for gene_name, _ in glur_sequences[:2]:
        for smiles_id, _ in smiles_data[:2]:
            example_files.append(f"{gene_name}_{smiles_id}.yaml")

    for example in example_files:
        print(f"  {example}")


if __name__ == "__main__":
    generate_yaml_files()
