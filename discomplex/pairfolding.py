
import pandas as pd
from Bio import SeqIO
import biolmai # from bioml.ai for folding with ESMfold via API. Needs API key specified in file .env!
import gzip
from Bio.PDB import PDBParser, Selection, NeighborSearch
import sys
import io
import os
from os.path import dirname

# Define PROJECT_HOME as two levels up from the current script
PROJECT_HOME = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = f'{PROJECT_HOME}/data/raw'
INTERIM_DIR = f'{PROJECT_HOME}/data/interim'
PROCESSED_DIR = f'{PROJECT_HOME}/data/processed'
TABLE_DIR = f'{PROJECT_HOME}/tables/singlechain_folding'


def load_fasta_sequences(fasta_file):
    """
    This function reads a FASTA formatted file and returns a dictionary
    It assumes that the second field in terms of separator '|' is a unique
    accession code that is used as key in the returned dictionary

    Example:
    sq = load_fasta_sequences(filename)
    sq['P12452'] # gets sequence of tht Uniprot id
    """
    sequences = {}
    
    # Check if the file is gzip compressed
    if fasta_file.endswith('.gz'):
        open_func = gzip.open
    else:
        open_func = open

    # Open the file accordingly
    print("Reading sequences in FASTA format from ", fasta_file)
    with open_func(fasta_file, "rt") as file:
        for record in SeqIO.parse(file, "fasta"):
            uniprot_id = record.description.split('|')[1]
            sequences[uniprot_id] = record
    return sequences

def read_tabular(file):
    if not os.path.exists(file):
        raise FileNotFoundError(file)
    df = None
    if file.endswith("tsv"):
        df = pd.read_csv(file, sep='\t')
    elif file.endswith("csv"):
        df = pd.read_csv(file)
    else:
        raise ValueError("Only formats CSV and TSV are supported.")
    return df


def protein_folding_predict_esm(seqs):
    """
    Perform single-chain or multi-chain protein folding using ESMfold
    For this to work, an API token has to be specified in file .env
    which looks like
    BIOLMAI_TOKEN=YOUR_API_KEY
    You can get the API key for free from bioml.ai after registering
    and generating a new key under "settings".
    """
    cls = biolmai.ESMFoldSingleChain() if len(seqs)==1 else biolmai.ESMFoldMultiChain()
    if isinstance(seqs, str):
        seqs = [seqs]
    if len(seqs) > 1:
        seqs = [':'.join(seqs)]
    resp = cls.predict(seqs)[0]

    print("Response object after biolmai/ESMfolding API call:", type(resp))
    # print(resp)
    print(resp.keys())
    if 'error' in resp:
        print("Obtained error message:", resp['error'])
        return {'response':resp}
    results = resp['results']

    pdb = results[0]['pdb']
    pdb = pdb.replace(r'\n', '\n') # fix issue: initially newline is provided as 2 characters of \ and n, but it should be one character \n
    return {'pdb_text':pdb, 'response':resp}


def protein_folding_predict(seqs, method='esm'):
    """
    Generic wrapper function to perform single-chain or multi-chain protein folding.
    Currently supported: method == 'esm'. In the future it could
    be AlphaFold etc.
    Please see for more details under protein_folding_predict_esm for requirements
    in terms of API tokens.
    """
    if method == 'esm':
        return protein_folding_predict_esm(seqs)
    raise ValueError("Suppported folding methods: esm, found instead " + method)


def find_interacting_residues(pdb_filename, cutoff=8.0,
    count_base=0, # 1 for 1-based counting of output regions, 0 for 0-based output
    ):
    """
    This function returns for a multi-strand complex all regions with respect to the
    first chain where its C-beta atom is closer than a cutoff to a C-beta 
    atom from any other chain (C-alpha in case of GLY)
    """
    # Parse the PDB file
    parser = PDBParser()
    print("Reading PDB structure", pdb_filename)
    structure = parser.get_structure("PDB_structure", pdb_filename)

    # Get the first chain
    first_chain = next(structure.get_chains())

    # Extract C-beta atoms (C-alpha for Glycine) from the first chain
    target_atoms = []
    for residue in first_chain:
        if residue.get_id()[0] == " ":  # Ignore hetero/water residues
            if "CB" in residue or residue.get_resname() == "GLY":
                atom_name = "CB" if "CB" in residue else "CA"  # Use CA for Glycine
                target_atoms.append(residue[atom_name])

    # Get all other chains
    other_chains = [chain for chain in structure.get_chains() if chain != first_chain]

    # Extract all C-beta atoms from other chains
    other_atoms = []
    for chain in other_chains:
        for residue in chain:
            if residue.get_id()[0] == " ":  # Ignore hetero/water residues
                if "CB" in residue:
                    other_atoms.append(residue["CB"])
                elif 'CA' in residue:
                    other_atoms.append(residue["CA"])

    # Search for interactions
    ns = NeighborSearch(other_atoms)
    interacting_regions = set()
    for atom in target_atoms:
        close_atoms = ns.search(atom.get_coord(), cutoff, "R")  # Search within the cutoff radius
        if close_atoms:
            interacting_regions.add(atom.get_parent().get_id()[1])

    # Convert to regions (start, end) tuples
    interacting_regions = sorted(list(interacting_regions))
    regions = []
    if interacting_regions:
        start = interacting_regions[0]
        end = start
        for pos in interacting_regions[1:]:
            if pos == end + 1:
                end = pos
            else:
                regions.append((start + count_base, end + count_base))
                start = end = pos
        regions.append((start + count_base, end + count_base))

    return regions


def fold_pairs(result_file = f'{INTERIM_DIR}/ppi_verified.tsv',
        sequence_file = f'{RAW_DIR}/uniprot/uniprot_sprot.fasta.gz',
        pdb_out_dir = f"{PROCESSED_DIR}/pdb_pairs_out",
        binding_output_file=f'{PROCESSED_DIR}/binding_predictions_tmp.tsv',
        n=10,
        folding_method='esm',
        random_state=42):
    if not os.path.exists(pdb_out_dir):
        print("Creating directory", pdb_out_dir)
        os.makedirs(pdb_out_dir)
    tbl = pd.read_csv(result_file, sep='\t')
    if random_state is not None and random_state > 0:
        tbl = tbl.sample(frac=1, random_state=random_state).reset_index(drop=True)
    sequence_dictionary = load_fasta_sequences(sequence_file)
    if n <= 0:
        n = len(tbl)
    
    binding_starts = []
    binding_ends = []
    ids1 = []
    ids2 = []
    for i in range(n):
        id1 = tbl.at[i,'uniprot1']
        id2 = tbl.at[i,'uniprot2']
        # disprot = tbl.at[i,'disprot']
        if id1 in sequence_dictionary and id2 in sequence_dictionary:
            seq1 = str(sequence_dictionary[id1].seq)
            seq2 = str(sequence_dictionary[id2].seq)
            structure_shortname = f'{id1}_{id2}'
            # print(seq1)
            # print(seq2)
            # print()
            seqs = [seq1, seq2]
            folding_results = protein_folding_predict(seqs, method=folding_method)
            pdb_output_file = pdb_out_dir + "/" + structure_shortname + f"_prediction_{folding_method}.pdb"    
            if not os.path.exists(pdb_output_file):
                if folding_results is not None and 'pdb_text' in folding_results and folding_results['pdb_text']: 
                    folded_pdb = folding_results['pdb_text']
                    with open(pdb_output_file, 'w') as f:
                        print("Writing PDB structure to file", pdb_output_file)
                        f.write(folded_pdb)
            if os.path.exists(pdb_output_file):
                regions = find_interacting_residues(pdb_output_file, cutoff=8.0,
                count_base=0) # 1 for 1-based counting of output regions, 0 for 0-based output
                for region in regions:
                    binding_starts.append(region[0])
                    binding_ends.append(region[1])
                    ids1.append(id1)
                    ids2.append(id2)
                    region_results = pd.DataFrame(data={'id1':ids1,'id2':ids2, 'binding_start':binding_starts, 'binding_end':binding_ends})
                    region_results.to_csv(binding_output_file + "_loop_tmp.tsv", sep='\t')
            else:
                print("Id pair", id1, id2, 'No PPi structure prediction obtained! not both in verified part of Uniprot/Swissprot!')

    region_results = pd.DataFrame(data={'id1':ids1,'id2':ids2, 'binding_start':binding_starts, 'binding_end':binding_ends})
    if not os.path.exists(os.path.dirname(binding_output_file)):
        os.makedirs(os.path.dirname(binding_output_file)) # create needed directory
    print("Writing predicted regions to file", binding_output_file)
    region_results.to_csv(binding_output_file, sep='\t')


def fold_single_chain(ids,
        ref_sequences = f'{RAW_DIR}/uniprot/uniprot_sprot.fasta.gz',
        pdb_out_dir = f"{PROCESSED_DIR}/pdb_singles_out",
        prediction_table_file=f"{TABLE_DIR}/disorder_predictions_folding.tsv",
        n=10,
        folding_method='esm', # esm or dmp
        random_state=42):
    """
    Single-chain folding of a set of uniprot_ids
    """
    if not os.path.exists(dirname(prediction_table_file)):
        print("Creating directory" , dirname(prediction_table_file))
        os.makedirs(dirname(prediction_table_file))
    if isinstance(ref_sequences, str):
        ref_sequences = load_fasta_sequences(ref_sequences)
    n = len(ids)
    ids = list(ids)
    # binding_starts = []
    # binding_ends = []
    # ids1 = []
    folded_paths=[]
    for i in range(n):
        id1 = ids[i]
        # disprot = tbl.at[i,'disprot']
        if id1 in ref_sequences:
            seq1 = str(ref_sequences[id1].seq)
            structure_shortname = f'{id1}'
            # print(seq1)
            # print(seq2)
            # print()
            seqs = [seq1]
            folding_results = protein_folding_predict(seqs, method=folding_method)
            pdb_output_file = pdb_out_dir + "/" + structure_shortname + f"_prediction_{folding_method}.pdb"    
            if not os.path.exists(pdb_output_file):
                if folding_results is not None and 'pdb_text' in folding_results and folding_results['pdb_text']: 
                    folded_pdb = folding_results['pdb_text']
                    with open(pdb_output_file, 'w') as f:
                        print("Writing PDB structure to file", pdb_output_file)
                        f.write(folded_pdb)
                if os.path.exists(pdb_output_file):
                    folded_paths.append(pdb_output_file)
    #         if os.path.exists(pdb_output_file):
    #             pdbdf = pdb_to_dataframe(pdb_output_file, count_base=0) # 1 for 1-based counting of output regions, 0 for 0-based output
    #             for region in regions:
    #                 binding_starts.append(region[0])
    #                 binding_ends.append(region[1])
    #                 ids1.append(id1)
    #                 ids2.append(id2)
    #                 region_results = pd.DataFrame(data={'id1':ids1,'id2':ids2, 'binding_start':binding_starts, 'binding_end':binding_ends})
    #                 region_results.to_csv(binding_output_file + "_loop_tmp.tsv", sep='\t')
    #         else:
    #             print("Id pair", id1, id2, 'No PPi structure prediction obtained! not both in verified part of Uniprot/Swissprot!')

    # dirname(prediction_table_file)
    # region_results = pd.DataFrame(data={'id1':ids1,'id2':ids2, 'binding_start':binding_starts, 'binding_end':binding_ends})
    # if not os.path.exists(os.path.dirname(binding_output_file)):
    #     os.makedirs(os.path.dirname(binding_output_file)) # create needed directory
    # print("Writing predicted regions to file", binding_output_file)
    # region_results.to_csv(binding_output_file, sep='\t')
    return {'folded_pdb_paths':folded_paths}

def main(
        instructions=f"{RAW_DIR}/disprot/disprot_2023-12.tsv",
        id_col='acc',
        mode = 'single',
        ppi_file = f'{INTERIM_DIR}/ppi_verified.tsv',
        ref_sequences = f'{RAW_DIR}/uniprot/uniprot_sprot.fasta.gz',
        output_dir = f'{PROCESSED_DIR}/folding',
        # output_table_file=None,
        n=10,
        folding_method='esm',
        random_state=42):
        
    if isinstance(instructions, str):
        instructions = read_tabular(instructions)
    if not id_col in instructions.columns:
        raise ValueError(f"Column {id_col} not found in instruction table")
    if not os.path.exists(output_dir):
        print("Creating output directory:", output_dir)
        os.makedirs(output_dir)
    if mode == 'single':
        pdb_out_dir = f"{output_dir}/pdb_singlechain_out"
        output_table_file =f'{output_dir}/disorder_predictions_tmp.tsv'
    elif mode == 'pair':
        pdb_out_dir = f"{output_dir}/pdb_pairs_out"
        output_table_file = f'{output_dir}/binding_predictions_tmp.tsv'
    ids = instructions[id_col].unique()
    if not os.path.exists(pdb_out_dir):
        print("Creating PDB output directory:", pdb_out_dir)
        os.makedirs(pdb_out_dir)
    if not isinstance(instructions, pd.DataFrame):
        raise ValueError("information about what to fold has to be in tabular format.")
    if mode == 'single':
        result = fold_single_chain(ids,
            ref_sequences = ref_sequences,
            pdb_out_dir = pdb_out_dir,
            prediction_table_file=output_table_file,
            n=10,
            folding_method='esm', # esm or dmp
            random_state=42)
    print("Paths folded PDB structures:", result)

if __name__ == '__main__':
    main(n=10)
