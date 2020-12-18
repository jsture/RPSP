# pro-hormone-predictor

## Finding Secreted Genes

To get the FASTA files for all reviewed, secreted, human genes (as in <code>secreted.fasta</code>), use the following request to the [UniProtKB API](https://www.uniprot.org/uniprot/):

> https://www.uniprot.org/uniprot/?query=locations:(location:%22Secreted%20[SL-0243]%22)+AND+reviewed:yes+AND+organism:%22Homo%20sapiens%20(Human)%20[9606]%22&format=fasta
