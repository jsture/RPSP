import os

# useful filepaths
INPUT_FASTA_PATH = f"{os.getcwd()}/secreted.fasta"
OUTPUT_FASTA_PATH = f"{os.getcwd()}/output/"

# Custom constraints for a prohormone:
# - Must have a signal peptide length of at least SIGNAL_PEPTIDE_LENGTH (cleavage sites in the first SIGNAL_PEPTIDE_LENGTH residues do not count)
# - Must have at least MIN_NUM_CLEAVAGE_SITES amino acids between each cleavage site in the sequence
# - Must have at least MIN_CLEAVAGE_SPACING KK/KR/RR/RK cleavage sites
SIGNAL_PEPTIDE_LENGTH = 20
MIN_NUM_CLEAVAGE_SITES = 4
MIN_CLEAVAGE_SPACING = 3
