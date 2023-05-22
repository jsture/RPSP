import fastaparser
import regex
import constants
import os
from tqdm import tqdm
import grequests
from utils import uniprot_from_fasta
from collections import Counter

# ignore matplotlib deprecation warnings
import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)
import matplotlib.pyplot as plt

class FastaParser:
  def __init__(self, filepath):
    '''
    Initializes a FastaParser object.
    The FastaParser class is used to parse a FASTA file containing multiple sequences,
    and determine which of those sequences meet the constraints for a prohormone.
    
    Parameters:
      filepath (str): The path to the FASTA file to be parsed.
    '''

    # initialize member variables
    self.sequences = self.matches = []
    self.cleavage_spacing = constants.MIN_CLEAVAGE_SPACING
    self.num_cleavage_sites = constants.MIN_NUM_CLEAVAGE_SITES

    # custom regex for matching cleavage sites with the correct spacing
    self.prohormone_regex = rf'(?:(?<=.{{{constants.SIGNAL_PEPTIDE_LENGTH}}})(?:(?<!K|R)(?:KK|KR|RR|RK)(?=[^RKILPVH]|(?<=KR)H|$))(?=(?:(?!(?R)).){{{self.cleavage_spacing},}}|$))'

    # open the FASTA file, read the sequences, and select matches based on the current criteria
    try:
      with open(filepath, 'r') as fasta_file:
        print('Reading FASTA file...')
        self.sequences = [(sequence, self.get_cleavage_sites(sequence)) for sequence in tqdm(fastaparser.Reader(fasta_file))]
        self.set_matches()
    except FileNotFoundError:
      print('The specified file cannot be read.')
    
    # writes peptide sequences from proteins to "output.fasta" file
    self.write_fasta("output.fasta")

  @property
  def cleavage_spacing(self):
    '''
    Returns the spacing between cleavage sites.
    '''
    return self._cleavage_spacing

  @cleavage_spacing.setter
  def cleavage_spacing(self, spacing):
    '''
    Sets the minimum cleavage spacing, under the condition that it is greater than the global MIN_CLEAVAGE_SPACING.
    
    Parameters:
      spacing (int): The new minimum cleavage spacing.
    '''
    self._cleavage_spacing = constants.MIN_CLEAVAGE_SPACING if spacing < constants.MIN_CLEAVAGE_SPACING else spacing

    # recompile regex  and update matches so that it corresponds with the new spacing
    self.prohormone_regex = rf'(?:(?<=.{{{constants.SIGNAL_PEPTIDE_LENGTH}}})(?:(?<!K|R)(?:KK|KR|RR|RK)(?=[^RKILPVH]|(?<=KR)H|$))(?=(?:(?!(?R)).){{{self.cleavage_spacing},}}|$))'
    if self.sequences:
      print('Updating matches...')
      self.sequences = [(sequence, self.get_cleavage_sites(sequence)) for sequence in tqdm([seq for seq, _ in self.sequences])]
    self.set_matches()

  @property
  def num_cleavage_sites(self):
    '''
    Returns the number of cleavage sites.
    '''
    return self._num_cleavage_sites

  @num_cleavage_sites.setter
  def num_cleavage_sites(self, num_cleavage_sites):
    '''
    Sets the number of cleavage sites, under the condition that it is greater than the global MIN_NUM_CLEAVAGE_SITES.

    Parameters:
      num_cleavage_sites (int): The new minimum number of cleavage sites.
    '''
    self._num_cleavage_sites = constants.MIN_NUM_CLEAVAGE_SITES if num_cleavage_sites < constants.MIN_NUM_CLEAVAGE_SITES else num_cleavage_sites
    
    # reset the matches to correspond with the new number of cleavage sites
    self.set_matches()

  def get_cleavage_sites(self, sequence):
    '''
    Returns a list of locations of cleavage sites in the sequence.

    Parameters:
      sequence (str): The sequence to be searched for cleavage sites.
    '''
    return list(regex.finditer(self.prohormone_regex, sequence.sequence_as_string()))

  def set_matches(self):
    '''
    Sets the matches attribute to the list of sequences that have the specified number of cleavage sites and the right spacing between them.
    '''
    self.matches = []
    if self.sequences:
      print(f'Finding matches with at least {self.num_cleavage_sites} cleavage sites separated by at least {self.cleavage_spacing} amino acids...')
      for sequence, match_locations in tqdm(self.sequences):
        if len(match_locations) >= self.num_cleavage_sites:
          self.matches.append((sequence, match_locations))
      
      print('Number of matches: ', len(self.matches))

  @property
  def num_matches(self):
    '''
    Returns the number of sequences in the FASTA file that are prohormones.
    '''
    return len(self.matches)

  def graph_cleavage_density(self):
    '''
    Creates a scatterplot of the number of cleavage sites in each sequence vs. the length of the sequence.
    Prohormones with a large number of cleavage sites and a shorter sequence length are typically desired.
    '''
    x, y = zip(*[(len(sequence), len(match)) for sequence, match in tqdm(self.matches)])
    plt.scatter(x, y)
    plt.title('Cleavage Density of Identified Prohormones')
    plt.xlabel('Gene Length (aa)')
    plt.ylabel('Number of Cleavage Sequences')
    plt.show()

  def tissue_specificity(self, show_graph=False):
    '''
    Returns a dictionary of tissue specificity values for each sequence.

    Parameters:
      show_graph (bool): Whether or not to show a graph of the tissue specificity values.
    '''

    # retrieve tissue expression data The Human Protein Atlas
    # API specification can be found at https://www.proteinatlas.org/about/help/dataaccess
    print('Fetching tissue expression annotations...')
    urls = [f'http://www.proteinatlas.org/api/search_download.php?search={uniprot_from_fasta(sequence)}&format=json&columns=g,up,rnats,rnatsm&compress=no' for sequence, _ in self.sequences]
    responses = grequests.map((grequests.get(url) for url in urls)) # fetch the data in parallel

    # for each tissue, calculate the number of sequences in the FASTA file that are expressed in that tissue
    print('Analyzing responses...')
    tissue_expression = Counter()
    for response in responses:
      if response and response.json():
        expression = response.json()[0]['RNA tissue specific nTPM']
        if expression:
          tissue_expression.update(list(expression.keys()))
        else:
          # "wide": the sequence affects a large number of tissues
          tissue_expression.update(['wide'])
      else:
        tissue_expression.update(['no HPA annotation'])
    
    # display a bar graph showing the number of sequences expressed in each tissue
    if show_graph:
      # sort from highest to lowest
      tissues, frequencies = zip(*sorted(tissue_expression.items(), key = lambda x: -x[1]))
      indices = [i for i in range(len(tissues))]
      width = 0.9
      plt.bar(indices, frequencies, width)
      plt.xticks(indices, tissues, rotation='vertical')
      plt.title('Secreted gene tissue expression')
      plt.ylabel('Number of Genes')
      plt.show()
    return tissue_expression
    
  def write_fasta(self, filename):
    '''
    Writes the matches to a FASTA file.
    
    Parameters:
      filename (str): The name of the file to be written.
    '''

    # create directory if it doesn't exist
    if not os.path.isdir(constants.OUTPUT_FASTA_PATH):
      os.mkdir(constants.OUTPUT_FASTA_PATH)
    
    # write the matches to the FASTA file in the output directory
    with open(f'{constants.OUTPUT_FASTA_PATH}{filename}', 'w') as fasta_file:
      writer = fastaparser.Writer(fasta_file)
      print('Writing matching FASTA sequences to file...')
      for match in tqdm(self.matches):
        peptides = self.split_fasta_sequence(match)
        for p in peptides:
          writer.writefasta(p)
  
  def split_fasta_sequence(fasta_sequence_obj, regex_matches):
    
    '''
    Splits fasta sequences up based on regex matches.
    
    Parameters:
      fasta_sequence: str, regex_matches: List[regex.Match]
    '''
    fasta_sequence = fasta_sequence_obj.sequences
    splits = []
    start_index = 0
    i = 0
    if regex_matches:
      for match in regex_matches[1]:
        if i==0:
          split_sequence = regex_matches[0][0:match.span()[0]]
          splits.append(split_sequence)
          start_index=match.span()[1]
          i=1
        else:
          split_sequence = regex_matches[0][start_index:match.span()[0]]
          splits.append(split_sequence)
          start_index=match.span()[1]
          if start_index < len(regex_matches[0]):
            split_sequence = regex_matches[0][start_index:]
            splits.append(split_sequence)  
      else:
        splits.append(regex_matches[0])
    return splits
