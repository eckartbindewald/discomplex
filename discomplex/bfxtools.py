from Bio import SeqIO
import os
import gzip

# Define PROJECT_HOME as two levels up from the current script
PROJECT_HOME = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = f'{PROJECT_HOME}/data/raw'
INTERIM_DIR = f'{PROJECT_HOME}/data/interim'
PROCESSED_DIR = f'{PROJECT_HOME}/data/processed'


def load_fasta_sequences(file:str, id_pos:int=1, sep="|"):
    """
    This function reads a FASTA formatted file and returns a dictionary
    It assumes that the second field in terms of separator '|' is a unique
    accession code that is used as key in the returned dictionary
    
    Args:
        file(str): Filename
        sep(str): The separator character to split sequence names. Defaults to "|"
        id_pos(int): position of index after split(sep) command

    
    """
    sequences = {}
    
    # Check if the file is gzip compressed
    if file.endswith('.gz'):
        open_func = gzip.open
    else:
        open_func = open

    # Open the file accordingly
    print("Reading sequences in FASTA format from ", file)
    with open_func(file, "rt") as file:
        for record in SeqIO.parse(file, "fasta"):
            uniprot_id = record.description.split(sep)[id_pos]
            sequences[uniprot_id] = record
    return sequences