#!/usr/bin/env python3
"""
Script to process FASTA files in stupid_plan/glurs directory and create a CSV file
with gene names and sequences.
"""

import re
import csv
from pathlib import Path


def extract_gene_name_and_sequence(fasta_file_path):
    """
    Extract gene name and sequence from a FASTA file.

    Args:
        fasta_file_path (str): Path to the FASTA file

    Returns:
        tuple: (gene_name, sequence) or (None, None) if parsing fails
    """
    try:
        with open(fasta_file_path, "r") as f:
            content = f.read().strip()

        lines = content.split("\n")
        if not lines:
            return None, None

        # First line should be the header
        header = lines[0]
        if not header.startswith(">"):
            return None, None

        # Extract gene name using regex to find GN=XXXX pattern
        gene_match = re.search(r"GN=(\w+)", header)
        if not gene_match:
            return None, None

        gene_name = gene_match.group(1)

        # Join all remaining lines to get the sequence
        sequence = "".join(lines[1:])

        return gene_name, sequence

    except Exception as e:
        print(f"Error processing {fasta_file_path}: {e}")
        return None, None


def process_glurs_directory():
    """
    Process all FASTA files in the stupid_plan/glurs directory and create a CSV file.
    """
    # Define paths (relative to project root)
    glurs_dir = Path("../glurs")
    output_file = "../mols/glurs_sequences.csv"

    # Check if directory exists
    if not glurs_dir.exists():
        print(f"Directory {glurs_dir} does not exist!")
        return

    # Get all .fasta files in the directory
    fasta_files = list(glurs_dir.glob("*.fasta"))

    if not fasta_files:
        print(f"No .fasta files found in {glurs_dir}")
        return

    # Process files and collect data
    gene_data = []

    print(f"Processing {len(fasta_files)} FASTA files...")

    for fasta_file in sorted(fasta_files):
        print(f"Processing: {fasta_file.name}")
        gene_name, sequence = extract_gene_name_and_sequence(fasta_file)

        if gene_name and sequence:
            gene_data.append((gene_name, sequence))
            print(f"  Found gene: {gene_name} (sequence length: {len(sequence)})")
        else:
            print(f"  Failed to extract data from {fasta_file.name}")

    # Create output directory if it doesn't exist
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Write to CSV file
    if gene_data:
        with open(output_file, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            # Write header
            writer.writerow(["gene_name", "sequence"])
            # Write data
            for gene_name, sequence in gene_data:
                writer.writerow([gene_name, sequence])

        print(f"\nSuccessfully created {output_file} with {len(gene_data)} entries")

        # Print summary
        print("\nSummary:")
        for gene_name, sequence in gene_data:
            print(f"  {gene_name}: {len(sequence)} amino acids")
    else:
        print("No valid gene data found to write to CSV file")


if __name__ == "__main__":
    process_glurs_directory()
