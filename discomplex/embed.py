
import pandas as pd
from Bio import SeqIO
import biolmai # from bioml.ai for folding with ESMfold via API. Needs API key specified in file .env!
import gzip
from Bio.PDB import PDBParser, Selection, NeighborSearch
import sys
import io
import os
from dotenv import load_dotenv
from pairfolding import load_fasta_sequences
import numpy as np
import requests
import time

# Load environment variables from .env file
load_dotenv()

# Now you can use the environment variable
BIOLMAI_TOKEN = os.getenv('BIOLMAI_TOKEN')

# Define PROJECT_HOME as two levels up from the current script
PROJECT_HOME = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = f'{PROJECT_HOME}/data/raw'
INTERIM_DIR = f'{PROJECT_HOME}/data/interim'
PROCESSED_DIR = f'{PROJECT_HOME}/data/processed'


def protein_embed_esm(sequence, model_slug = 'esm2-650m',
    base_url='https://biolm.ai/api/v2'):
    """
    Compute per-residue embeddings using a variant of the ESM-2
    model.
    """    
    url = f"{base_url}/{model_slug}/encode/"

    headers = {
        "Content-Type": "application/json",
        "Authorization": f"Token {BIOLMAI_TOKEN.strip()}",
    }

    data = {
        "params": {
            "include": [
                "per_token", # was 'mean'
                "contacts",
                "logits",
                "attentions"
            ]
        },
        "items": [
            {
                "sequence": sequence, # "MSILVTRPSPAGEELVSRLRTLGQVAWHFPLIEFSPGQQLPQLADQLAALGESDLLFALSQHAVAFAQSQLHQQDRKWPRLPDYFAIGRTTALALHTVSGQKILYPQDREISEVLLQLPELQNIAGKRALILRGNGGRELIGDTLTARGAEVTFCECYQRCAIHYDGAEEAMRWQAREVTMVVVTSGEMLQQLWSLIPQWYREHWLLHCRLLVVSERLAKLARELGWQDIKVADNADNDALLRALQ",

            },
        ]
    }
    print(f"Submitting following request to {base_url}")
    print(data)
    # Make the request!
    response = requests.post(
        url=url,
        headers=headers,
        json=data,
    ).json()
    # print(response)
    print(type(response))
    if isinstance(response,dict) and 'error' in response:
        raise ValueError(f"An error raise in API call. Check format of submitted data but also status of API Token maintained at the biolm.ai user account: {response['error']}.")
    print(len(response))
    print(response.keys())
    if 'results' not in response:
        raise ValueError("Missing expected keyword 'results' in response dictionary: {")
    results = response['results']
    if not isinstance(results, list):
        raise(f"Strange, expected a 'list' datatype for response[results]: {type(results)}")
    if len(results) ==0:
        raise(f"Strange, expected at list one list item response[results]: {len(results)}")
    if len(results) > 1:
        print(f"Minor warning: only the first result will be analyzed, received {len(results)}")
    results = results[0]
    #  user can go further:   embeddings = results['representations']; logits = results['logits'] etc
    return results
    

    # ˚

    #     cls = biolmai.ESMFoldSingleChain() if len(seqs)==1 else biolmai.ESMFoldMultiChain()
    #     if isinstance(seqs, str):
    #         seqs = [seqs]
    #     if len(seqs) > 1:
    #         seqs = [':'.join(seqs)]
    #     resp = cls.predict(seqs)[0]

    #     print("Response object after biolmai/ESMfolding API call:", type(resp))
    #     # print(resp)
    #     print(resp.keys())
    #     if 'error' in resp:
    #         print("Obtained error message:", resp['error'])
    #         return {'response':resp}
    #     results = resp['results']

    #     pdb = results[0]['pdb']
    #     pdb = pdb.replace(r'\n', '\n') # fix issue: initially newline is provided as 2 characters of \ and n, but it should be one character \n
    #     return {'pdb_text':pdb, 'response':resp}


def protein_embed(sequence, model = {
    'name':'esm2-650m', 'family':'esm'}):
    """
    Wrapper function for potentially different methods beyond ESM
    (like AlphaFold, RosettaFold) etc.
    Currently only method 'esm' is implemented.
    """
    if 'family' not in model or model['family'] != 'esm':
        raise ValueError("Currently only method 'esm' is supported model family, but found " + model)
    return protein_embed_esm(sequence, model_slug=model['name'])
    


def run_embed_loop(
        accessions,
        table_gt,
        sequence_dictionary,
        n=1,
        model = {
            'name':'esm2-650m',
            'family':'esm',
            'dim':1280

        },
        retry_max=3,
        wait_per_iteration=10.0):
    if n <= 0:
        n = len(accessions)
    else:
        n = min(n, len(accessions))
    groundtruth_all = [] # list over all positions of all proteins visited in order
    embedding_vecs = [] # list over all embedding vectors for all positions of all visited proteins
    for i in range(n):
        acc = accessions[i]
        print("Working on row", i, )
        assert acc in sequence_dictionary
        seq = str(sequence_dictionary[acc].seq)
        seq_len = len(seq)
        # get all disordered region for this protein chain:
        ref_regions = table_gt[table_gt['acc']== acc].reset_index(drop=True)
        groundtruth_vec = np.zeros([seq_len])
        # we are now defining an outcome vector for binary classification:
        print(ref_regions.head().to_string())
        for j in range(len(ref_regions)):
            start = ref_regions.at[j, 'start'] - 1 # back to 0-based counting for start
            assert start >= 0
            stop = ref_regions.at[j, 'end']
            assert stop <= seq_len
            groundtruth_vec[start:stop]=1.0
        success = False
        for retry in range(retry_max):
            try:
                embedding_results = protein_embed(seq, model=model)
                success = True
                break
            except Exception as e:
                print(f"API call resulted in exception: {e}")
                print("Retry attempt", retry + 1)
                time.sleep(wait_per_iteration)
        if not success:
            print(f"Warning: could not generate results for row {i}: {acc}")
            continue
        print(embedding_results.keys())
        if not isinstance(embedding_results, dict) or 'representations' not in embedding_results:
            raise("Strange, no embeddings defined")
        embeddings = embedding_results['representations']
        if isinstance(embeddings, dict) and len(embeddings) == 1:
            _embed_keys = list(embeddings.keys()) # name of neural network layer like '33' as string
            embeddings = embeddings[_embed_keys[0]]
        if len(embeddings) != seq_len or not isinstance(embeddings, list):
            raise("Inconsistent number of embeddings vectors, expected a lists os lists so that there is one embedding vector per residue")
        for pos in range(len(embeddings)):
            embedding_vecs.append(embeddings[pos])
        groundtruth_all.extend(groundtruth_vec)
        assert len(embedding_vecs) == len(groundtruth_all), f'The length of the ground-truth vector and the number of embedding vectors did not match for {acc}'
        time.sleep(wait_per_iteration)
    return {'embeddings':embedding_vecs, 'groundtruth':groundtruth_all}


def main(groundtruth_file = f'{RAW_DIR}/disprot/disprot_2023-12.tsv',
        sequence_file = 'data/raw/uniprot/uniprot_sprot.fasta.gz',
        data_output_file=f'{PROCESSED_DIR}/disprot/disprot_esm_embed__NUMROWS__.npz',
        n=1,
        model = {
            'name':'esm2-650m',
            'family':'esm',
            'dim':1280
        },
        random_state=42 ):
    if not os.path.exists(groundtruth_file):
        raise FileNotFoundError("Could not find reference file with ground-truth from disprot: " + groundtruth_file)
    tbl_gt = pd.read_csv(groundtruth_file, sep='\t')
    if random_state is not None and random_state >= 0: # shuffle
        tbl_gt = tbl_gt.sample(frac=1.0, replace=False, random_state=random_state)
    sequence_dictionary = load_fasta_sequences(sequence_file)
    acc_uniq = list(tbl_gt['acc'].unique()) # .tolist()
    print(acc_uniq)
    print(type(acc_uniq))
    assert isinstance(acc_uniq, list)
    acc_seq_uniq= list(set(sequence_dictionary.keys()))
    acc_common = list(set(acc_uniq) & set(acc_seq_uniq))
    tbl_gt = tbl_gt[tbl_gt['acc'].isin(acc_common)].reset_index(drop=True)
    if n <= 0:
        n = len(acc_common)
    else:
        n = min(n, len(acc_common))
    results = run_embed_loop(
        accessions = acc_common,
        table_gt=tbl_gt,
        sequence_dictionary=sequence_dictionary,
        n=n,
        model=model)
    assert isinstance(results, dict)
    assert 'embeddings' in results
    assert 'groundtruth' in results
    embedding_vecs = results['embeddings']
    vector_dim = len(embedding_vecs[0]) # length of individual embedding vectors
    groundtruth_all = results['groundtruth']
    num_vectors = len(embedding_vecs)
    embedding_vecs = np.vstack(embedding_vecs)
    assert embedding_vecs.shape == (num_vectors, vector_dim)
    groundtruth_np = np.zeros( [ len(groundtruth_all), 1] )
    print("shape of groundtruth_np:", groundtruth_np.shape)
    print("'shape' of groundtruth_all list:", len(groundtruth_all), type(groundtruth_all[0]))
    assert groundtruth_np.shape == (num_vectors, 1)
    assert groundtruth_np.shape == (len(groundtruth_all), 1)
    groundtruth_np[:, 0] = groundtruth_all
    # store as npz:
    data_output_file = data_output_file.replace("_NUMROWS__", str(len(groundtruth_all)))
    outdir = os.path.dirname(data_output_file)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    print("Storing Embedding vectors and matching ground truth in file:", data_output_file)
    print("One can load the data with:")
    print(f"data = np.load({data_output_file} )")
    print("Feature data:")
    print("X=data['X']")
    print("Ground-truth binary outcome (disprot binding site yes or now)")
    print("y=data['y']")
    np.savez(data_output_file, X=embedding_vecs, y=groundtruth_np)
    assert embedding_vecs.shape[0] == groundtruth_np.shape[0], 'Number of rows does not match between X and y!'


if __name__ == '__main__':
    main(n=100)