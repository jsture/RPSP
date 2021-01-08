import fastaparser
import regex
import constants
import os
import matplotlib.pyplot as plt
from tqdm import tqdm
import grequests
from utils import uniprot_from_fasta
from collections import Counter

class FastaParser:
  def __init__(self, filepath):
    self.sequences = self.matches = []
    self.cleavage_spacing = constants.MIN_CLEAVAGE_SPACING
    self.num_cleavage_sites = constants.MIN_NUM_CLEAVAGE_SITES
    self.prohormone_regex = rf'(?:(?<=.{{{constants.SIGNAL_PEPTIDE_LENGTH}}})(?:(?<!K|R)(?:KK|KR|RR|RK)(?=[^RKILPVH]|(?<=KR)H|$))(?=(?:(?!(?R)).){{{self.cleavage_spacing},}}|$))'

    try:
      with open(filepath, 'r') as fasta_file:
        print('Reading FASTA file...')
        self.sequences = [(sequence, self.get_cleavage_sites(sequence)) for sequence in tqdm(fastaparser.Reader(fasta_file))]
        self.set_matches()
    except FileNotFoundError:
      print('The specified file cannot be read.')

  @property
  def cleavage_spacing(self):
    return self._cleavage_spacing

  @cleavage_spacing.setter
  def cleavage_spacing(self, spacing):
    self._cleavage_spacing = constants.MIN_CLEAVAGE_SPACING if spacing < constants.MIN_CLEAVAGE_SPACING else spacing
    # recompile regex
    self.prohormone_regex = rf'(?:(?<=.{{{constants.SIGNAL_PEPTIDE_LENGTH}}})(?:(?<!K|R)(?:KK|KR|RR|RK)(?=[^RKILPVH]|(?<=KR)H|$))(?=(?:(?!(?R)).){{{self.cleavage_spacing},}}|$))'
    if self.sequences:
      print('Updating matches...')
      self.sequences = [(sequence, self.get_cleavage_sites(sequence)) for sequence in tqdm([seq for seq, _ in self.sequences])]
    self.set_matches()

  @property
  def num_cleavage_sites(self):
    return self._num_cleavage_sites

  @num_cleavage_sites.setter
  def num_cleavage_sites(self, num_cleavage_sites):
    self._num_cleavage_sites = constants.MIN_NUM_CLEAVAGE_SITES if num_cleavage_sites < constants.MIN_NUM_CLEAVAGE_SITES else num_cleavage_sites
    self.set_matches()

  def get_cleavage_sites(self, sequence):
    return list(regex.finditer(self.prohormone_regex, sequence.sequence_as_string()))

  def set_matches(self):
    self.matches = []
    if self.sequences:
      print(f'Finding matches with at least {self.num_cleavage_sites} cleavage sites separated by at least {self.cleavage_spacing} amino acids...')
      for sequence, match_locations in tqdm(self.sequences):
        if len(match_locations) >= self.num_cleavage_sites:
          self.matches.append((sequence, match_locations))
      
      print('Number of matches: ', len(self.matches))

  @property
  def num_matches(self):
    return len(self.matches)

  def graph_cleavage_density(self):
    x, y = zip(*[(len(sequence), len(match)) for sequence, match in tqdm(self.matches)])
    plt.scatter(x, y)
    plt.title('Cleavage Density of Identified Prohormones')
    plt.xlabel('Gene Length (aa)')
    plt.ylabel('Number of Cleavage Sequences')
    plt.show()

  def tissue_specificity(self, show_graph=False):
    print('Fetching tissue expression annotations...')
    urls = [f'http://www.proteinatlas.org/api/search_download.php?search={uniprot_from_fasta(sequence)}&format=json&columns=g,up,rnats,rnatsm&compress=no' for sequence, match in self.sequences]
    responses = grequests.map((grequests.get(url) for url in urls))

    print('Analyzing responses...')
    tissue_expression = Counter()
    for response in responses:
      if response and response.json():
        expression = response.json()[0]['RNA tissue specific NX']
        if expression:
          tissue_expression.update(list(expression.keys()))
        else:
          tissue_expression.update(['wide'])
      else:
        tissue_expression.update(['no HPA annotation'])

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
    if not os.path.isdir(constants.OUTPUT_FASTA_PATH):
      os.mkdir(constants.OUTPUT_FASTA_PATH)
    
    with open(f'{constants.OUTPUT_FASTA_PATH}{filename}', 'w') as fasta_file:
      writer = fastaparser.Writer(fasta_file)
      print('Writing matching FASTA sequences to file...')
      for match in tqdm(self.matches):
        writer.writefasta(match)
