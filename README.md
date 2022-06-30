# pro-hormone-predictor

This program predicts whether a secreted gene has prohormone activity based on the number of cleavage sites it contains.

## Finding Secreted Genes

To get the FASTA files for all reviewed, secreted, human genes (as in <code>secreted.fasta</code>), use the following request to the [UniProtKB API](https://www.uniprot.org/uniprot/):

> https://www.uniprot.org/uniprot/?query=locations:(location:%22Secreted%20[SL-0243]%22)+AND+reviewed:yes+AND+organism:%22Homo%20sapiens%20(Human)%20[9606]%22&format=fasta

## Running the Program (Unix/macOS)
First ensure that you have Python 3.x and pip installed. Then, create and activate a virtual environment using the following commands:
```
python3 -m venv env
source env/bin/activate
```
Install all of the requirements using the following command:
```
python3 -m pip install -r requirements.txt
```
You are now set to run any of the code in this repository.
To fetch and download the FASTA file for reviewed (i.e., recreate <code>secreted.fasta</code>), run
```
python3 download_fasta.py
```
It would be good practice to do this periodically to get the most updated set of secreted genes.
To analyze the secreted genes for potential prohormone activity, run the main program, which you can edit for your desired functionality:
```
python3 main.py
```
