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

# Define PROJECT_HOME as two levels up from the current script
PROJECT_HOME = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW_DIR = f'{PROJECT_HOME}/data/raw'
INTERIM_DIR = f'{PROJECT_HOME}/data/interim'
PROCESSED_DIR = f'{PROJECT_HOME}/data/processed'
TABLE_DIR = f'{PROJECT_HOME}/tables/disorder'
PROTEIN_DOWNLOAD_LIMIT = 10000
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

def run_disprot_structure(disprot23,
    no_download=True,
    sasa_scale = 0.01,
    URL_BASE = "https://alphafold.ebi.ac.uk/files/",
    scratchdir=os.path.join(PROJECT_HOME, 'data-raw/afdb')):

    # scratchdir = os.path.join(PROJECT_HOME, 'data-raw/afdb')
    if not os.path.exists(scratchdir):
        os.makedirs(scratchdir, exist_ok=True)
        print(f"Created scratch directory at: {scratchdir}")

    # Parameters


    # Load dataset (example using a CSV file, adjust as necessary)
    # print(disprot23.head(20))  # Similar to R's kable function for displaying data frames
    # Initialize a list to store data dictionaries
    data_list = []
    # Accession codes (unique IDs)
    accessions = sorted(disprot23['acc'].unique())[:PROTEIN_DOWNLOAD_LIMIT]
    used_accessions = []
    # Main loop for downloading and analyzing protein structures
    protein_count = 0
    
    for acc in accessions:

        dest_name = f"AF-{acc}-F1-model.pdb"
        destfile = os.path.join(scratchdir, dest_name)
        if not os.path.exists(destfile):
            if no_download:
                continue
            for version in [4,3,2,1]:
                fname = f"AF-{acc}-F1-model_v{version}.pdb"
                url = f"{URL_BASE}{fname}"
                print(f"Trying download of {acc} from {url}...")
                response = requests.get(url)
                if response.status_code == 200:
                    with open(destfile, 'wb') as f:
                        f.write(response.content)
                    break # success, no need to try further
            if not os.path.exists(destfile):
                print(f"Failed to download {acc}")
                continue

        # Read PDB file
        parser = PDBParser()
        protein_count += 1
        used_accessions.append(acc)
        structure = parser.get_structure(acc, destfile)

        # Calculate SASA using FreeSASA
        # structure_sasa = freesasa.Structure(destfile)
        # result = freesasa.calc(structure_sasa)

        sasa_result, dummy = freesasa.calcBioPDB(structure)
                
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

        # Example of extracting CA atoms, requires Bio.PDB and processing
        for model in structure:
            for chain in model:
                chain_id = chain.id
                for residue in chain:
                    chain_length = len(chain)
                    res_name = residue.resname.strip() 
                    if 'CA' in residue:
                        ca_atom = residue['CA']
                        residue_id = residue.get_id()[1]  # Getting residue number
                        rel_residue_id = residue_id/chain_length
                        log_residue_id = math.log10(residue_id)
                        log_rev_residue_id = math.log10(chain_length-residue_id+1)
                        aa_code= seq1(res_name)
                        b_factor = ca_atom.get_bfactor()
                        # print("Accessing chain id, residue id:", chain_id, residue_id, type(residue_id))
                        sasa_residue = sasa_result.residueAreas()[chain_id][str(residue_id)]
                        polar = sasa_residue.polar*sasa_scale
                        rel_polar = sasa_residue.relativePolar
                        apolar = sasa_residue.apolar*sasa_scale
                        rel_apolar = sasa_residue.relativeApolar
                        main_chain = sasa_residue.mainChain*sasa_scale
                        rel_main_chain = sasa_residue.relativeMainChain
                        side_chain = sasa_residue.sideChain*sasa_scale
                        rel_side_chain = sasa_residue.relativeSideChain
                        total = sasa_residue.total*sasa_scale
                        rel_total = sasa_residue.relativeTotal
                        # sasa_residue.residueType
                        has_rel_areas = int(sasa_residue.hasRelativeAreas)

                        structured = 1  # Default to structured
                        # Check if this residue is within any disordered region
                        disregions = disprot23[(disprot23['acc'] == acc) & (disprot23['start'] <= residue_id) & (disprot23['end'] >= residue_id)]
                        if not disregions.empty:
                            structured = 0  # Mark as disordered
                        # sasa = sasa_per_residue[f'residue {chain.id} {residue_id}']
                        # Append dictionary of data to the list
                        data_list.append({'Accession': acc,
                            'AACode':aa_code,
                            'Confidence': b_factor,
                            'ResNo': residue_id/100,
                            'RelResNo':rel_residue_id,
                            'lgResNo':log_residue_id,
                            'lgRevResNo':log_rev_residue_id,
                            'Polar':polar, 'RelPolar':rel_polar,
                            'Apolar':apolar, 'RelApolar':rel_apolar,
                            'MainChain':main_chain, 'RelMainChain':rel_main_chain,
                            'SideChain':side_chain, 'RelSideChain':rel_side_chain,
                            'Total':total, 'RelTotal':rel_total,
                            'structured': structured
                            })

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
    show=True, allowed_regions=['disorder']):
    disprot23 = pd.read_csv(disprot_file, sep='\t')  # Assuming data in CSV for this example
    print(disprot23.columns)
    print(disprot23['term_name'].value_counts())
    disprot23 = disprot23[disprot23['term_name'].isin(allowed_regions)].reset_index(drop=True)
    results = run_disprot_structure(disprot23)
    disorder = results['disorder']
    protein_count = len(results['accessions'])
    if not os.path.exists(TABLE_DIR):
        print("Creating output directory for table data:", TABLE_DIR)
        os.makedirs(TABLE_DIR)
    outfile = os.path.join(TABLE_DIR, f'disorder_{protein_count}.tsv')
    print("Writing to output file", outfile)
    disorder.to_csv(outfile, sep='\t')
    numeric_df = disorder.select_dtypes(include=[np.number])
    corr_mtx = numeric_df.corr()
    if show:
        plot_heatmap(corr_mtx)
        plot_violin(disorder)
        plot_structure_features(disorder,x='Polar', y='Confidence')
        plot_structure_features(disorder,x='Apolar', y='Confidence')
        plot_structure_features(disorder,x='Polar', y='Apolar')
        plot_structure_features(disorder,x='Total', y='Confidence')
        plot_structure_features(disorder,x='RelPolar', y='Confidence')
        plot_structure_features(disorder,x='RelApolar', y='Confidence')
        plot_structure_features(disorder,x='RelPolar', y='RelApolar')
        plot_structure_features(disorder,x='RelTotal', y='Confidence')

if __name__ == '__main__':
    main()