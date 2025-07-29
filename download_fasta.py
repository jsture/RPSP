# Python script to fetch and download the reviewed, secreted human genes from the UniProt API.

import requests
import constants

# make the request to the UniProt API
url = 'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28cc_scl_term%3ASL-0243%29%29+AND+%28reviewed%3Atrue%29+AND+%28model_organism%3A9606%29'
response = requests.get(url)

# save the results to a FASTA file
with open(constants.INPUT_FASTA_PATH, 'w') as f:
  f.write(response.text)
