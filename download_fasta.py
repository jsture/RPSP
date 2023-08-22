# Python script to fetch and download the reviewed, secreted human genes from the UniProt API.

import requests
import constants

# make the request to the UniProt API
url = 'https://rest.uniprot.org/uniprotkb/search?query=(cc_scl_term:SL-0243)%20AND%20(reviewed:true)%20AND%20Human&format=fasta&size=500'
response = requests.get(url)

# save the results to a FASTA file
with open(constants.INPUT_FASTA_PATH, 'w') as f:
  f.write(response.text)
