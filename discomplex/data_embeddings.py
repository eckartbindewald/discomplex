
import embed
import os
import sys
import numpy as np
import pandas as pd
from bfxtools import load_fasta_sequences
import pickle

DEBUG = True

# Define PROJECT_HOME as two levels up from the current script
PROJECT_HOME = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = f'{PROJECT_HOME}/data/raw'
INTERIM_DIR = f'{PROJECT_HOME}/data/interim'
PROCESSED_DIR = f'{PROJECT_HOME}/data/processed'

def main(
        sequence_file = f'{RAW_DIR}/interpro/protein-matching-IPR000719.fasta',
        outbase=f'{PROCESSED_DIR}/interprot/esm_embed_IPR000719',
        n=10,
        token_mode="mean", # use "per_token" or "mean",
        model = {
            'name':'esm2-650m',
            'family':'esm',
            'dim':1280
        },
        random_state=42 ):
    sequence_dictionary = load_fasta_sequences(sequence_file, id_pos=0)
    print(sequence_dictionary)
    accessions = list(sequence_dictionary.keys())
    for i in range(len(accessions)):
        accessions[i] = accessions[i].split("|")[0]
    print('accessions:', accessions)
    # acc_uniq = list(tbl_gt['acc'].unique()) # .tolist()
    if n <= 0:
        n = len(accessions)
    else:
        n = min(n, len(accessions))

    results = embed.run_embed_loop(
        accessions = accessions,
        sequence_dictionary=sequence_dictionary,
        n=n,
        token_mode=token_mode,
        model=model)
    outdir = os.path.dirname(outbase)
    if not os.path.exists(outdir):
        print("creating directory", outdir)
        os.makedirs(outdir)
    used_accessions = results['accessions']
    pickle_file = f'{outbase}_{len(used_accessions)}.pickle'
    with open(pickle_file, 'wb') as f:
        print("Writing to output file in 'pickle' format:",pickle_file )
        pickle.dump(results, f)
    if DEBUG:
        print("Results returned from embed.run_embed_loop:")
        print(results)
    assert isinstance(results, dict)
    if not 'embeddings' in results:
        raise "Expecting key 'embeddings' in server result"
    
    embedding_vecs = results['embeddings']
    
    print("Type embeddings", type(embedding_vecs), len(embedding_vecs))
    if isinstance(embedding_vecs, list) and len(embedding_vecs) == 0:
        print("no embedding vectors defined!")
        sys.exit(0)
    print("embedding_vecs:", results['embeddings'])
    vector_dim = len(embedding_vecs[0]) # length of individual embedding vectors
    # groundtruth_all = results['groundtruth']
    num_vectors = len(embedding_vecs)

    embedding_vecs = np.vstack(embedding_vecs)
    assert embedding_vecs.shape == (num_vectors, vector_dim)
    # groundtruth_np = np.zeros( [ len(groundtruth_all), 1] )
    # print("shape of groundtruth_np:", groundtruth_np.shape)
    # print("'shape' of groundtruth_all list:", len(groundtruth_all), type(groundtruth_all[0]))
    # assert groundtruth_np.shape == (num_vectors, 1)
    # assert groundtruth_np.shape == (len(groundtruth_all), 1)
    # groundtruth_np[:, 0] = groundtruth_all
    # store as npz:
    # data_output_file = data_output_file.replace("_NUMROWS__", str(len(groundtruth_all)))
     
    data_output_file = f'{outbase}.npz'
    print("Storing Embedding vectors and matching ground truth in file:", data_output_file)
    print("One can load the data with:")
    print(f"data = np.load({data_output_file} )")
    print("Feature data:")
    print("X=data['X']")
    print("Used accesion codes:")
    print("y=data['accessions']")
    print("Used sequences:")
    print("y=data['sequences']")
    np.savez(data_output_file, X=embedding_vecs, accessions=used_accessions)
    # assert embedding_vecs.shape[0] == groundtruth_np.shape[0], 'Number of rows does not match between X and y!'

if __name__ == "__main__":
    main(n=-1)