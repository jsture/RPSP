import constants
from FastaParser import FastaParsergit 

if __name__ == '__main__':
  parser = FastaParser(constants.INPUT_FASTA_PATH)
  
  # counts = parser.tissue_specificity(show_graph=True)
  # parser.graph_cleavage_density()