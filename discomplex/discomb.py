#!/usr/bin/env python

import discomplex as dc
import os
import glob
import sys
from collections import defaultdict
import argparse
import bfxtools
from bfxtools import read_plain_sequences
import pandas as pd
from data_disprotstructure import structurefile_to_features
from test_disorder_from_confidence import run_disorder_prediction

# Define PROJECT_HOME as two levels up from the current script
PROJECT_HOME = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = f'{PROJECT_HOME}/data'
RAW_DIR = f'{PROJECT_HOME}/data/raw'
INTERIM_DIR = f'{PROJECT_HOME}/data/interim'
# INTERIM_DIR = f'{PROJECT_HOME}/data/interim'
PROCESSED_DIR = f'{PROJECT_HOME}/data/processed'
TABLE_DIR = f'{PROJECT_HOME}/tables/disorder'
MODEL_DIR = f'{PROJECT_HOME}/models/disorder'


def find_similar_uniprots(sequence_file, sequence_db):
    pass


def main_legacy(
    sequence_file,
    input_dir=os.path.join(PROCESSED_DIR,'disorder'),
    model_name='dense',
    model_filenames = None,
    model_dir = MODEL_DIR, 
    window_size=5,
    shuffle_columns=False,
    regions=['flexible linker/spacer',
    'molecular adaptor activity', 'molecular function regulator',
    'phosphorylation display site',
    'lipid binding',
    'ion binding',
    'pre-molten globule',
    'nucleic acid binding','disorder', 'protein binding', 'disorder to order']):
    print("Reading test data from directory", input_dir)
    print("Reading models from directory", model_dir)
    for region in regions:
        region_name = region.replace("/","_").replace(" ", "_")
        test_filenames = glob.glob(f"{input_dir}/disorder_{region_name}_*_test.tsv")
        print(f"Test file names for region {region}:", test_filenames)
        if model_filenames is None:
            model_filenames = glob.glob(f"{model_dir}/disorder_{region_name}_*.keras")
        for test_file in test_filenames:
            for model_file in model_filenames:
                print(f"Starting testing model {basename(model_file)} using data from {basename(test_file)} window-size: w{window_size}")
                if not os.path.exists(model_file):
                    print("Warning: could not find model file at", model_file)
                    continue
                try:
                    main(df_test=test_file,model_name=model_name, model_file=model_file,
                    window_size=window_size,
                    shuffle_columns=shuffle_columns)
                except Exception as e:
                    print(f"Encountered exception during evaluation: {e}")
                    raise e





def verify_file_path(path: str, description: str):
    """
    Verifies that a given file path exists.
    Raises a FileNotFoundError with a custom error message if the file does not exist.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"{description} file at {path} could not be found.")

def parse_cli_args(args):
    """
    Parse command line arguments using argparse.
    """
    parser = argparse.ArgumentParser(description="Process the protein sequence files and databases.")

    # Define expected command-line arguments
    parser.add_argument("-i", "--input", required=True, help="Specifies the input protein sequence in FASTA format.")
    parser.add_argument("-d", "--db", required=True, help="Specifies the file location of the sequence database searchable by `hmmscan`.")
    parser.add_argument("-s", "--structures", required=True, help="The directory containing protein 3D structures - typically AlphaFoldDB.")
    parser.add_argument("-m", "--model_dir", default="models", help="The directory containing TensorFlow models used for this application. Defaults to 'models'.")

    # Parse arguments
    return vars(parser.parse_args(args))


def find_homologs_hmmscan(sequence,
    db_location, method='hmmscan') -> pd.DataFrame:
    result = bfxtools.hmmer_scan(db_location, queryseq=sequence)
    if not isinstance(result, pd.DataFrame) or not 'Accession' in result.columns:
        raise("Internal error: did not obtain parseable result table from hmmscan")
    accessions = result['Accession'].tolist()
    uniprots = [accession.split('.')[0] for accession in accessions]
    result['Uniprot'] = uniprots
    return result


def find_closest_uniprot_ids(sequence, method='hmmscan') -> pd.DataFrame:
    """
    Find closest known uniprot sequence id
    """
    if method == 'hmmscan':
        return find_closest_uniprot_ids_hmmscan(sequence)
    else:
        raise ValueError("Unkonwn sequence homolog identifican tion method. Supported is 'hmmscan'")


def predict_disorder_end_to_end(sequence_file:str,
    sequence_db:str,
    model_dir:str, structure_db:str=None,
    seq_method:str='hmmscan',
    alphafold_methods:list=['file']): # 'file', http
    """
    Overall prediction of disordered regions in proteins
    """
    result = {}
    # stage 1: find closest sequence id
    homologs = bfxtools.hmmer_scan(sequence_db, queryfile=sequence_file)
    print(homologs)
    result['homologs'] = homologs
    print("Identified domains of homologous proteins domains:")
    print(result['homologs'])
    pdb:str = None # string containing single-chain PDB structure
    for idx, row in homologs.iterrows():
        accession = row['accession'].split('.')[0]
        pdb = bfxtools.fetch_alphafold_db(accession=accession,
                    methods=alphafold_methods,
                    file_root=structure_db)
        if pdb is not None and 'content' in pdb and pdb['content']:
            break
    # if PDB is not None: parse to BioPython
    print("pdb variable:", pdb)
    pdbfile = None
    if pdb:
        pdbfile = pdb['file']
    if not pdb or not pdbfile:
        raise Exception("Could not find any structure of homologous protein sequence. Either file system root or base-URL of AlphaFold Database has to be specified")
    print("Using PDB file", pdbfile)
    assert ',' not in pdbfile
    features = structurefile_to_features(pdbfile)

    result = run_disorder_prediction(df_test, model_file, scaler_file=None,
        window_size=17,
        model_name='dense',
        shuffle_columns=False, # can be a list of column names
        drop_columns=['Unnamed: 0'], # , 'ResNo'],
        )

    return result


def main():
    # Parse command line arguments
    params = parse_cli_args(sys.argv[1:])  # exclude the script name itself

    # Validate file paths
    verify_file_path(params['input'], "Input protein sequence")
    verify_file_path(params['db'], "Sequence database")
    verify_file_path(params['structures'], "Structure database")
    
    sequence_file = params['input'] # read_plain_sequences('input')
    sequence_db = params['db']
    model_dir = params['model_dir']
    structure_db = params['structures']
    results = []
    for sequence in sequences:
        results.append(predict_disorder_end_to_end(sequence_file,    
            sequence_db=sequence_db, model_dir=model_dir, structure_db=structure_db))
    # Accessing model_dir with default value handled by argparse
    
    print("### RESULTS: ###")
    for i in range(len(results)):
        print(f"### Sequence {i+1}:")
        print(results[i])
        print("-----")
        print()


def test_predict_disorder_end_to_end():
    sequence_file = f"{DATA_DIR}/fixtures/P03265.fa"
    sequence_db = f"{RAW_DIR}/pfam/Pfam-A.hmm"
    # seq = read_plain_sequences(sequence_file)[0]
    # print("Read sequence:", seq)
    model_dir = MODEL_DIR
    structure_db = f'{INTERIM_DIR}/afdb'
    result = predict_disorder_end_to_end(sequence_file,
        sequence_db=sequence_db,
        model_dir=model_dir,
        structure_db=structure_db, seq_method='hmmscan')
    print(result)
    assert(isinstance(result, dict))
    assert(isinstance(result['homologs'], pd.DataFrame))
    assert(len(result['homologs']) > 0)


if __name__ == "__main__":
    test_predict_disorder_end_to_end()
    assert False
    main()