import constants
from FastaParser import FastaParser


# parse the FASTA file and set the prohormone matches
parser = FastaParser(constants.INPUT_FASTA_PATH)

# determine the tissue expression profiles of the prohormone matches
counts = parser.tissue_specificity(show_graph=True)

# for each match, graph its number of cleavage sites against its length
parser.graph_cleavage_density()
