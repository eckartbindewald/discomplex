"""
This file generates training data for prediting 
disordered regions in protein such that:
columns: feature data for an amino acid, ground truth for same position, plus reference information as text

"""
import pandas as pd
import os
import requests
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
import matplotlib.pyplot as plt
import seaborn as sns
import freesasa
import numpy as np
import math
from stride import run_stride
from params import aa_encoding

# import MDAnalysis as mda
# import psntools.core as core
# from MDAnalysis.analysis import secondary_structure

# Define PROJECT_HOME as two levels up from the current script
PROJECT_HOME = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = f'{PROJECT_HOME}/data/raw'
INTERIM_DIR = f'{PROJECT_HOME}/data/interim'
PROCESSED_DIR = f'{PROJECT_HOME}/data/processed/disorder'
# TABLE_DIR = f'{PROJECT_HOME}/tables/disorder'
PROTEIN_DOWNLOAD_LIMIT = 100000
STRIDE = "bin/stride"
DIGITS=5

# def extract_sasa(residue_areas):
#     sasa_data = {}
#     for residue_id, residue_area in residue_areas.items():
#         try:
#             # Assume the object has a 'total' attribute or method
#             sasa_data[residue_id] = residue_area.total
#         except AttributeError:
#             print(f"Failed to extract SASA for residue {residue_id}")
#     return sasa_data


def heatmap(correlation_matrix):
    # Heatmap of the correlation matrix
    plt.figure(figsize=(10, 8))
    sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt=".2f", linewidths=.5)
    plt.title('Correlation Heatmap')
    plt.show()


def regions_to_residues(starts, ends):
    """ Convert regions to a list of residues """
    result = []
    for start, end in zip(starts, ends):
        result.extend(range(start, end + 1))
    return result


def structure_to_features(structure, secstruct:pd.DataFrame,
    disprot:pd.DataFrame=None, verbose=True):
    """
    Computes features for disorder prediction for a given
    PDB structure
    structure: PDB structure of single-chain protein as parsed by BioPython
    secstructure: DataFrame in inhouse format as geerated by run_stride(pdbfile, 'stride')
    disprot(pd.DataFrame): DataFrame corresponding to reference ground truth
        data, with columns as provided by the Disprot database.
    """
    # mda_u = mda.Universe(destfile)
    # sse = secondary_structure.SegmentSecondaryStructure.Universe(mda_u)
    # psn = core.PSN("matrix.dat", mda_u)
    # print(psn)
    
    print(secstruct.head())
    if 'chain' not in secstruct.columns:
        raise Exception("Strange input format:", type(secstruct), str(secstruct))
    pdb_chain_names = secstruct['Chain'].unique()
    if verbose:
        print("PDB Chain names:", pdb_chain_names)
    assert isinstance(secstruct, pd.DataFrame), 'Expected data frame but got ' + str(secstruct) 
    assert len(secstruct) > 0
    # print(secstruct.head())
    # print(secstruct.shape)

    # for res in mda_u.residues:
    #     print(f"Residue {res.resname} {res.resid} has secondary structure: {sse.secondary_structure[res.index]}")
    # assert False
    # Calculate SASA using FreeSASA
    # structure_sasa = freesasa.Structure(destfile)
    # result = freesasa.calc(structure_sasa)

    # sasa_result, dummy = freesasa.calcBioPDB(structure)
            
    # print(
    #     sasa_result.residueAreas()['A']['1'].polar*sasa_scale,
    #     sasa_result.residueAreas()['A']['1'].relativePolar,
    #     sasa_result.residueAreas()['A']['1'].apolar*sasa_scale,
    #     sasa_result.residueAreas()['A']['1'].relativeApolar,
    #     sasa_result.residueAreas()['A']['1'].mainChain*sasa_scale,
    #     sasa_result.residueAreas()['A']['1'].relativeMainChain,
    #     sasa_result.residueAreas()['A']['1'].sideChain*sasa_scale,
    #     sasa_result.residueAreas()['A']['1'].relativeSideChain,
    #     sasa_result.residueAreas()['A']['1'].total*sasa_scale,
    #     sasa_result.residueAreas()['A']['1'].relativeTotal,
    #     sasa_result.residueAreas()['A']['1'].residueType,
    #     sasa_result.residueAreas()['A']['1'].hasRelativeAreas)
    # print(extract_sasa(sasa_result.residueAreas()))
    # print(obj.structures)
    # print(obj.relativePolar(), obj.relativeApolar(),
    #    obj.relativeMainChain(), obj.elativeSideChain())
    # print(sasa_result.residueAreas()['A']['1'])
    data_list = []
    # Example of extracting CA atoms, requires Bio.PDB and processing
    for model in structure:
        for chain in model:
            chain_id = chain.id
            assert chain_id in pdb_chain_names, 'Unexpected chain name:' + chain_id + " : " + ','.join(pdb_chain_names)
            for residue in chain:
                chain_length = len(chain) # strand name
                res_name = residue.resname.strip() 
                if 'CA' in residue:
                    ca_atom = residue['CA']
                    residue_id = residue.get_id()[1]  # Getting residue number
                    rel_residue_id = residue_id/chain_length
                    log_residue_id = math.log10(residue_id)
                    log_rev_residue_id = math.log10(chain_length-residue_id+1)
                    aa_code= seq1(res_name) # one-letter code
                    b_factor = ca_atom.get_bfactor()
                    sec = secstruct[(secstruct['ResPdbId']==residue_id) \
                        & (secstruct['Chain']==chain_id)].reset_index(drop=True)
                    # print("Accessing chain id, residue id:", chain_id, residue_id, type(residue_id))
                    # sasa_residue = sasa_result.residueAreas()[chain_id][str(residue_id)]
                    # polar = sasa_residue.polar*sasa_scale
                    # rel_polar = sasa_residue.relativePolar
                    # apolar = sasa_residue.apolar*sasa_scale
                    # rel_apolar = sasa_residue.relativeApolar
                    # main_chain = sasa_residue.mainChain*sasa_scale
                    # rel_main_chain = sasa_residue.relativeMainChain
                    # side_chain = sasa_residue.sideChain*sasa_scale
                    # rel_side_chain = sasa_residue.relativeSideChain
                    # total = sasa_residue.total*sasa_scale
                    # rel_total = sasa_residue.relativeTotal
                    # # sasa_residue.residueType
                    # has_rel_areas = int(sasa_residue.hasRelativeAreas)

                    structured = 1  # Default to structured
                    # Check if this residue is within any disordered region
                    if disprot is not None:
                        disregions = disprot[(disprot['acc'] == acc) & (disprot['start'] <= residue_id) & (disprot['end'] >= residue_id)]
                        if not disregions.empty:
                            structured = 0  # Mark as disordered
                    # sasa = sasa_per_residue[f'residue {chain.id} {residue_id}']
                    # Append dictionary of data to the list
                    if len(sec) == 1:
                        for _, r in sec.iterrows():
                            aa_size, polarity,charge,hydrophobicity,aromaticity=aa_encoding.get(aa_code,(0,0,0,0,0))
                            data_list.append({'Accession': acc,
                                'AACode':aa_code,
                                'aa_size':aa_size,
                                'polarity':polarity,
                                'charge':charge,
                                'hydrophobicity':hydrophobicity,
                                'aromaticity':aromaticity,
                                'Confidence': b_factor/100,
                                'ResNo': residue_id/100,
                                'RelResNo':round(rel_residue_id, DIGITS),
                                'lgResNo':round(log_residue_id,DIGITS),
                                'lgRevResNo':round(log_rev_residue_id, DIGITS),
                                'Phi':round(r['Phi']*angle_scale, DIGITS),
                                'Psi':round(r['Psi']*angle_scale, DIGITS),
                                'SecStruct':r['StructCode'],
                                # 'Polar':polar, 'RelPolar':rel_polar,
                                # 'Apolar':apolar, 'RelApolar':rel_apolar,
                                # 'MainChain':main_chain, 'RelMainChain':rel_main_chain,
                                # 'SideChain':side_chain, 'RelSideChain':rel_side_chain,
                                'TotalArea':round(r['Area']*area_scale, DIGITS),
                                'structured': structured
                                })
                    else:
                        print("Warning: could not identify secondary structure - is binary of program 'stride' installed? ",
                            len(sec), res_name, residue_id, chain_id)
                        if verbose:
                            print(secstruct.head())
                            print(secstruct.tail())
                            print((secstruct['ResPdbId']==residue_id))
                            print(secstruct['Chain']==chain_id)
                            print((secstruct['ResPdbId']==residue_id) & (secstruct['Chain']==chain_id))
                            print(secstruct[(secstruct['ResPdbId']==residue_id) \
                                & secstruct['Chain']==chain_id].reset_index(drop=True))
                            print(len(secstruct[(secstruct['ResPdbId']==residue_id) \
                                & (secstruct['Chain']==chain_id)].reset_index(drop=True)))
                            
                        secstruct.to_csv('tmp_secstruct.tsv', sep='\t')
                        assert False
    return { 'disorder': pd.DataFrame(data_list)}


def structurefile_to_features(pdbfile, 
    disprot:pd.DataFrame=None, verbose=True,
    accession="unknown"):
    """
    Converts PDB/MMCIF file to a tabular representation of features
    Wrapper for function structure_to_features
    """
    assert pdbfile is not None
    assert isinstance(pdbfile, str)
    assert ',' not in pdbfile
    # Read PDB file
    parser = PDBParser()
    structure = parser.get_structure(accession, pdbfile)
    secstruct:pd.DataFrame = run_stride(pdbfile, 'bin/stride')
    return structure_to_features(structure, secstruct, disprot=disprot, verbose=verbose)


def run_disprot_structure(disprot23,
    n=-1, # max number of accession codes
    no_download=False,
    area_scale = 0.01,
    angle_scale = 0.01, # angles between -1.8 and + 1.8 instead of the pi stuff :-)
    base_url = "https://alphafold.ebi.ac.uk/files/",
    file_template="AF-{acc}-F1-model_v{version}.pdb",
    file_template2="AF-{acc}-F1-model.pdb",
    scratchdir=os.path.join(INTERIM_DIR, 'afdb'),
    verbose=True):

    # scratchdir = os.path.join(PROJECT_HOME, 'data-raw/afdb')
    if not os.path.exists(scratchdir):
        os.makedirs(scratchdir, exist_ok=True)
        print(f"Created scratch directory at: {scratchdir}")


    # Load dataset (example using a CSV file, adjust as necessary)
    # print(disprot23.head(20))  # Similar to R's kable function for displaying data frames
    # Initialize a list to store data dictionaries
    data_list = []
    # Accession codes (unique IDs)
    accessions = sorted(disprot23['acc'].unique())
    
    if len(accessions) > PROTEIN_DOWNLOAD_LIMIT:
        print(f"Clipping maximum number of accession codes from {len(accessions)} to ", PROTEIN_DOWNLOAD_LIMIT)
        accessions = accessions[:PROTEIN_DOWNLOAD_LIMIT]
    used_accessions = []
    # Main loop for downloading and analyzing protein structures
    protein_count = 0
    if n <= 0 or n > len(accessions):
        n = len(accessions)
    for acc_iter in range(n):
        acc = accessions[acc_iter]

        dest_name = file_template2.format(acc=acc)
        assert ',' not in dest_name
         # f"AF-{acc}-F1-model.pdb"
        destfile = os.path.join(scratchdir, dest_name)
        print(f"Working on accession code {acc} : {acc_iter} out of {n} proteins")
        if not os.path.exists(destfile):
            if no_download:
                print(f"deactivated download, have to skip {acc}")
                continue
            for version in [4,3,2,1]:
                assert isinstance(acc, str)
                assert isinstance(version, str) or isinstance(version, int)
                fname = file_template.format(acc=acc, version=version) # f"AF-{acc}-F1-model_v{version}.pdb"
                assert ',' not in fname
                url = f"{base_url}{fname}"
                response = requests.get(url)
                if response.status_code == 200:
                    with open(destfile, 'wb') as f:
                        f.write(response.content)
                    print(f"Successfully downloaded {acc} from {url} and saved to {destfile}!")
                    if not os.path.exists(destfile):
                        raise FileNotFoundError('Strange, file should have been written to ' + destfile)
                    break # success, no need to try further
            if not os.path.exists(destfile):
                print(f"Failed to download {acc}")
                continue
        else:
            print("Using cached structure at", destfile)
        protein_count += 1
        used_accessions.append(acc)
        accession_result = structurefile_to_features(destfile, disprot=disprot23)['disprot']
        for idx, row in accession_result.iterrows():
            data_list.append(row)

    return {'disorder':pd.DataFrame(data_list),
        'accessions':used_accessions}

def plot_violin(ca_atoms):
    # Plotting (Example using Matplotlib)
    ca_atoms['Structured'] = ca_atoms['structured'].map({0: 'Disordered', 1: 'Structured'})
    plt.figure(figsize=(10, 6))
    ax = sns.violinplot(x='Structured', y='Confidence', data=ca_atoms)
    ax.set_xlabel('Local Protein Structure')
    ax.set_ylabel('AlphaFold Confidence Value')
    plt.show()

def plot_structure_features(df,x,y, hue='structured'):
    # Create a scatter plot
    plt.figure(figsize=(10, 6))
    scatter_plot = sns.scatterplot(data=df, x=x, y=y, hue=hue, palette='viridis', s=100)

    # Enhancing the plot
    scatter_plot.set_title(f'Scatter Plot of {x} vs {y} Colored by Structured')
    plt.legend(title='Structured')

    # Show plot
    plt.show()


def plot_heatmap(correlation_matrix):
    # Heatmap of the correlation matrix
    plt.figure(figsize=(10, 8))
    sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt=".2f", linewidths=.5)
    plt.title('Correlation Heatmap')
    plt.show()


def main(disprot_file=f'{RAW_DIR}/disprot/disprot_2023-12.tsv',
    n=-1,
    show=False, allowed_regions=['flexible linker/spacer',
    'molecular adaptor activity', 'molecular function regulator',
    'phosphorylation display site',
    'lipid binding','ion binding',
    'pre-molten globule','nucleic acid binding',
    'disorder', 'protein binding', 'disorder to order'],
    accession_col='acc',
    random_state=42,
    train_frac=0.8,
    count_test_region_min=5,
    no_download=False,
    verbose=False): # 
    """
    Overall counts from original disprot data:
    disorder 6292
    protein binding 1413
    disorder to order 738
    flexible linker/spacer 377
    molecular adaptor activity 271
    molecular function regulator 205
    phosphorylation display site 170
    lipid binding 99
    ion binding 90
    pre-molten globule 82
    nucleic acid binding 80
    molecular condensate scaffold activity 77
    DNA binding 75
    molecular function inhibitor activity 71
    molecular function activator activity 70
    order to disorder 69
    amyloid fibril formation 65
    RNA binding 60
    """
    disprot_orig = pd.read_csv(disprot_file, sep='\t')
    print(disprot_orig.columns)
    disprot_orig = disprot_orig.sample(frac=1, replace=False, random_state=random_state)
    assert isinstance(disprot_orig, pd.DataFrame)
    if n > 0 and n < len(disprot_orig):
        print("Truncating reference file from", len(disprot_orig), 'to', n, 'data rows.')
        disprot_orig = disprot_orig.head(n=n)
    accessions_uniq = disprot_orig[accession_col].unique()
    
    print(accessions_uniq)
    n_acc = len(accessions_uniq)
    n_acc_train = round(train_frac * n_acc)
    accessions_uniq_train = accessions_uniq[:n_acc_train]
    accessions_uniq_test = accessions_uniq[n_acc_train:]
    print("number of accession codes (total, train, test):", n_acc, n_acc_train)
    disprot_orig_train = disprot_orig[disprot_orig[accession_col].isin(accessions_uniq_train)]
    disprot_orig_test = disprot_orig[disprot_orig[accession_col].isin(accessions_uniq_test)]
    print(disprot_orig['term_name'].value_counts().to_string())
    for allowed_region in allowed_regions:
        allowed_region_name = allowed_region.replace("/", "_").replace(" ", "_")
        disprot23_train = disprot_orig_train[disprot_orig_train['term_name'].isin([allowed_region])].sort_values(by=[accession_col]).reset_index(drop=True)
        disprot23_test = disprot_orig_test[disprot_orig_test['term_name'].isin([allowed_region])].sort_values(by=[accession_col]).reset_index(drop=True)
        if len(disprot23_test) < count_test_region_min:
            print("Insufficient number of test cases for region", allowed_region, ":", len(disprot23_test))
            continue
        outfile_train = os.path.join(PROCESSED_DIR, f'disorder_{allowed_region_name}_{len(disprot23_train)}_{len(disprot23_test)}_train.tsv')
        outfile_test = os.path.join(PROCESSED_DIR, f'disorder_{allowed_region_name}_{len(disprot23_train)}_{len(disprot23_test)}_test.tsv')
        if os.path.exists(outfile_train) and os.path.exists(outfile_test):
            print("Both train and test data files already exists, skipping:")
            print(outfile_train)
            print(outfile_test)
            continue
        print("Working on region", allowed_region, outfile_train, outfile_test, '...')
        results_train = run_disprot_structure(disprot23_train, n=n, no_download=no_download, verbose=verbose)
        results_test = run_disprot_structure(disprot23_test, n=n, no_download=no_download, verbose=verbose)
        disorder_train = results_train['disorder']
        disorder_test = results_test['disorder']
        print("keys from results_train:", list(results_train.keys()))
        # protein_count = len(results_train["accessions"])
        if not os.path.exists(PROCESSED_DIR):
            print("Creating output directory for processed train&test data:", PROCESSED_DIR)
            os.makedirs(PROCESSED_DIR)

        print("Writing to output file", outfile_train, outfile_test)
        disorder_train.to_csv(outfile_train, sep='\t')
        disorder_test.to_csv(outfile_test, sep='\t')
        
        numeric_df = disorder_train.select_dtypes(include=[np.number])
        # corr_mtx = numeric_df.corr()
    # if show:
    #     plot_heatmap(corr_mtx)
    #     plot_violin(disorder)
        # plot_structure_features(disorder,x='Polar', y='Confidence')
        # plot_structure_features(disorder,x='Apolar', y='Confidence')
        # plot_structure_features(disorder,x='Polar', y='Apolar')
        # plot_structure_features(disorder,x='Total', y='Confidence')
        # plot_structure_features(disorder,x='RelPolar', y='Confidence')
        # plot_structure_features(disorder,x='RelApolar', y='Confidence')
        # plot_structure_features(disorder,x='RelPolar', y='RelApolar')
        # plot_structure_features(disorder,x='RelTotal', y='Confidence')

if __name__ == '__main__':
    main(n=-1, no_download=False)