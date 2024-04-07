import os
import pandas as pd
import urllib.request

# Define PROJECT_HOME as two levels up from the current script
PROJECT_HOME = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = f'{PROJECT_HOME}/data/raw'

FILENAMES_RAW = {
    'complexportal':'idmapping_ComplexPortal.tsv',
    'pdb':'idmapping_PDB.tsv',
    'disprot':'idmapping_disprot.tsv'
}


def process_data(raw_names:dict=FILENAMES_RAW ):
    """Processes data from specified raw data files.

    This function reads data from a set of predefined files (or those provided in the
    `raw_names` dictionary). These files are expected to contain mappings between different
    identifiers like UniProt and DisProt IDs. It then identifies associations between
    these identifiers within protein complexes.

    Args:
        raw_names: A dictionary where keys are directory names and values are file names.
                   These files should be tab-separated and contain identifier mappings.
                   Default is FILENAMES_RAW, a predefined dictionary.

    Returns:
        A pandas DataFrame containing columns for complex IDs, paired UniProt IDs,
        corresponding DisProt IDs, and associated UniProt IDs for each complex.
        The DataFrame is structured with columns: 'complex', 'uniprot1', 'disprot1',
        and 'uniprot2'.

    Raises:
        AssertionError: If the specified directories or files do not exist.

    Note:
        The function assumes the existence of specific directories and files as defined in RAW_DIR
        and FILENAMES_RAW. It also assumes a specific file format and structure for the input data.
    """
    tables = {}
    for d in raw_names.keys():
        indir = os.path.join(RAW_DIR, d)
        infile = os.path.join(indir, raw_names[d])
        print("Data input file:", infile)
        assert os.path.exists(indir)
        assert os.path.exists(infile)
        
        tables[d] = pd.read_csv(infile, sep='\t', header=None)
        if tables[d].shape[1] == 3:
            tables[d] = tables[d].iloc[:,[0,2]]
        tables[d].columns = ['from_id', 'to_id']
        print("Read file", infile, 'with shape', tables[d].shape)
        print(tables[d].head(3))

    # strategy:
    # loop over all complexes from ComplexPortal
    # find all uniprot ids
    # from there find all disprot ids
    cp = tables['complexportal']
    dp = tables['disprot']
    dlist = []
    for cpi in range(len(cp)):  # loop over all complexes
        cpid = cp.at[cpi, 'to_id']
        cp2 = cp[cp['to_id'] == cpid].reset_index(drop=True)
        dp2 = dp[dp['from_id'].isin(cp2['from_id'])].reset_index(drop=True)

        for j in range(len(dp2['from_id'])):
            unip = dp2.at[j, 'from_id']
            dpid = dp2.at[j, 'to_id']
            for unip2 in cp2['from_id'].tolist():
                if unip2 == unip:
                    continue
                print("Working on complex", cpid,unip,f'({dpid})', unip2)
                dlist.append({'complex':cpid,'uniprot1':unip,'disprot1':dpid, 'uniprot2': unip2})
    return pd.DataFrame(dlist)


def download_data(
    uniprot_fasta_url=\
    'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz',
    target_dir = os.path.join(RAW_DIR, 'uniprot')):
    '''
    Download standard data files (in this case UniProt/SwissProt protein sequences)
    and store them under data/raw
    '''
    # Create target directory if it doesn't exist
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)
    # Extract the filename from the URL
    file_name = uniprot_fasta_url.split('/')[-1]
    # Full path for the downloaded file
    file_path = os.path.join(target_dir, file_name)
    # Downloading the file from `uniprot_fasta_url`
    print(f"Downloading {file_name} to {file_path}")
    urllib.request.urlretrieve(uniprot_fasta_url, file_path)
    print(f"Download complete: {file_path}")


if __name__ == '__main__':
    df = process_data()
    outfile = 'data/intermediate/ppi_verified.tsv'
    print("Writing", len(df), "rows to output file", outfile)
    df.to_csv(outfile, sep='\t')
    ids = df['uniprot1'].unique().tolist()
    ids_str= '\n'.join(ids) + '\n'
    outfile = 'data/intermediate/disprot_unitprot_ids.txt'
    print("writing",len(ids), 'uniprot ids supported by ComplexPortal and Disprot to', outfile)
    with open(outfile, 'w') as f:
        f.write(ids_str)
    download_data()
    