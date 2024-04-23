from Bio import SeqIO
import os
import gzip
import pandas as pd
import tempfile
import io
import subprocess
import requests

# Define PROJECT_HOME as two levels up from the current script
PROJECT_HOME = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = f'{PROJECT_HOME}/data'
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


def read_plain_sequences(
    file:str, format:str="fasta") -> list[str]:
    """
    Reads biomolecular sequences in FASTA and other formats

    Args:
        file(str): Input filename
        format(str): Format of seuqence file. Defaults to 'fasta'.
            Other formats are supported according to biopython.
    """
    result = []
    # Check if the file is gzip compressed
    if file.endswith('.gz'):
        open_func = gzip.open
    else:
        open_func = open
    with open_func(file, "rt") as file:
        for record in SeqIO.parse(file, format):
            result.append(record.seq)
    return result


def parse_hmmer_scan(file_path):
    columns = [
        "target_name", "accession", "query_name", "query_accession",
        "E-value", "score", "bias", "E-value_domain", "score_domain", "bias_domain",
        "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", "description"
    ]
    # Read the file, skipping comments, and applying the column names
    df = pd.read_csv(
        file_path, 
        delim_whitespace=True, 
        comment='#', 
        header=None, 
        names=columns,
        usecols=range(len(columns))  # Adjust number of columns read based on actual output
    )
    return df


def fetch_alphafold_db_http(accession:str,
    pdb_scratch_dir:str,
    base_url:str='https://alphafold.ebi.ac.uk/files',
    fileformat='pdb',
    name_template="AF-{accession}-F1-model_v{version}.{fileformat}"
    ):
    """
    Fetches AlphaFoldDB structure based on ID in order
    of specified methods.
    Returns dictionary, diction key 'cntent' being string of PDB of MMCIF file if found and 'file' being filename
    """
    # file = None
    content = None
    dest_name = "AF-{accession}-F1-model.{fileformat}".format(accession=accession, fileformat=fileformat)
    destfile = os.path.join(pdb_scratch_dir, dest_name)
    print("Destination filename", dest_name)
    print("Destination path", destfile)
    print(f"Working on HTTP-retrieval of AlphaFold structure with accession code {accession}")
    if not os.path.exists(destfile): 
        if not os.path.exists(pdb_scratch_dir):
            print("Creating directory for caching PDB files.")
            os.makedirs(pdb_scratch_dir)
        for version in [4,3,2,1]:
            fname = name_template.format(accession=accession, version=version, fileformat=fileformat)
            url = f"{base_url}{fname}"
            response = requests.get(url)
            if response.status_code == 200:
                with open(destfile, 'wb') as f:
                    f.write(response.content)
                print(f"Successfully downloaded {accession} from {url} and saved to {destfile}!")
                if not os.path.exists(destfile):
                    raise FileNotFoundError('Strange, file should have been written to ' + destfile)
                break # success, no need to try further
        if not os.path.exists(destfile):
            return { 'error':f"Failed to download {accession}" }
    else:
        print("Using cached structure at", destfile)

    return {'content':content, 'file':destfile}


def find_path(filename, root_dir) -> str:
    """
    Returns string if existing under root. Potentially searches for matching
    path.
    """
    assert isinstance(filename, str)
    assert isinstance(root_dir, str)
    filename = os.path.basename(filename)
    if os.path.exists(os.path.join(root_dir, filename)):
        return os.path.join(root_dir, filename)
    print("TODO : implement search")
    return None


def fetch_alphafold_db_file(accession:str,
    file_root:str,
    name_template="AF-{accession}-F1-model_v{version}.{fileformat}",
    version=4,
    fileformat='pdb'):
    """
    Fetches AlphaFoldDB structure based on ID in order
    assuming that the file is under the specified file root.
    Returns string of PDB of MMCIF file if found
    """
    file = None
    content = None
    fname = name_template.format(accession=accession, version=version,  
        fileformat=fileformat)
    file_path = find_path(fname, file_root)
    print("Trying to find path:", fname, file_root, 'found:', file_path)
    if file_path is not None:
        with open(file, 'r') as f:
            content = f.read()
    assert ',' not in file
    return { 'content':content, 'file':file_path }


def fetch_alphafold_db(accession:str,
    methods:list=['file','http'],
    file_root:str=None,
    pdb_scratch_dir:str='pdb_scratch',
    base_url:str="https://alphafold.ebi.ac.uk/files",
    fileformat='pdb'):
    """
    Fetches AlphaFoldDB structure based on ID in order
    of specified methods.
    """
    pdb = None
    for method in methods:
        if method == 'file':
            if file_root is None:
                raise ValueError("Method file-search specified by no root directory for AlphaFold database provided")
            pdb = fetch_alphafold_db_file(accession=accession,
                file_root=file_root,
                fileformat=fileformat)
        elif method == 'http':
            pdb = fetch_alphafold_db_http(accession=accession,
                pdb_scratch_dir=pdb_scratch_dir,
                base_url=base_url,
                fileformat=fileformat)
        if pdb is not None:
            break
    assert isinstance(pdb, dict)
    return pdb


def hmmer_scan(db_location, 
    queryfile=None, queryseq=None) -> pd.DataFrame:
    """
    Runs program `hmmerscan` from HMMER suite of programs
    to find a query sequence in a database of sequences.
    Similar to BLAST but more sensitive and only slightly
    slower.
    """
    if queryseq is None and queryfile is None:
        raise FileNotFoundError("Either query sequence or query file (in FASTA format) must be specified.")

    if queryseq is not None and queryfile is None:
        with tempfile.NamedTemporaryFile('w', delete=False) as query_file:
            query_file.write(queryseq)
            queryfile = query_file.name

    with tempfile.NamedTemporaryFile('w', delete=False) as tblout_file:
        temp_tblout_path = tblout_file.name

    command = ["hmmscan", "--domtblout", temp_tblout_path, db_location, queryfile]
    print("hmmscan command:")
    print(' '.join(command))
    result = subprocess.run(command, text=True, capture_output=True)

    if result.returncode != 0:
        print("Error in running hmmscan:", result.stderr)
        return pd.DataFrame()

    return parse_hmmer_scan(temp_tblout_path)



def test_hammer_scan(queryfile=f"{DATA_DIR}/fixtures/P03265.fa",
    db_location=f"{RAW_DIR}/pfam/Pfam-A.hmm"):
    assert os.path.exists(queryfile)
    assert os.path.exists(f'{db_location}.h3m')
    result = hmmer_scan(db_location, queryfile=queryfile )
    print(result.to_string())
    assert isinstance(result, pd.DataFrame)
    assert len(result) > 0


if __name__ == "__main__":
    test_hammer_scan()

