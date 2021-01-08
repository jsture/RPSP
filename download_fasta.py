'''
Python script to fetch and download the reviewed, secreted human genes
from the UniProt API.
'''

import requests
import constants

if __name__ == '__main__':
  url = 'https://www.uniprot.org/uniprot/?query=locations:(location:%22Secreted%20[SL-0243]%22)+AND+reviewed:yes+AND+organism:%22Homo%20sapiens%20(Human)%20[9606]%22&format=fasta'
  response = requests.get(url)
  with open(constants.INPUT_FASTA_PATH, 'w') as f:
    f.write(response.text)