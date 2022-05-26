def uniprot_from_fasta(fasta):
  '''
  Given a sequence in FASTA format, returns the corresponding Uniprot ID.
  
  Parameters:
    fasta (str): The sequence in FASTA format.
  '''
  return fasta.id[3:9]
