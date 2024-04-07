# Evaluating Folding of Protein Dimers

This toolkit provides utilities for loading protein sequences, predicting protein folding structures using the ESMfold method via the bioml.ai API, and analyzing PDB structures for protein-protein interactions.
Disclaimer: This implementation is NOT finished!

## Dependencies

- Python 3.x
- Biopython
- Pandas
- Biolmai (bioml.ai API)

Make sure you have the above dependencies installed. You can install them using pip:

```bash
pip install pipenv
pipenv install
```

## Configuration

### API Key

To use the bioml.ai API for protein folding predictions, you need to set your API key in a `.env` file:

BIOLMAI_TOKEN=YOUR_API_KEY

You can obtain an API key by registering at bioml.ai and generating a new key under the "Settings" section.

### Input Data

Ensure your protein sequence data is in FASTA format. If your FASTA file is gzip-compressed (`.gz`), the scripts can handle it automatically.

## Usage

The idea is to run two scripts: i) data preparation with data_prepare.py and actual model application for predicting sites with pairfolding.py. The scripts should be run with `pipenv` instead of directly.

### Preparing data

Again, as preliminary create the virtual environment and install needed packages with following 2 lines of code:

```
# cd to project home directory
pip install pipenv
pipenv install
```

Actually running the script that prepares input data:

```
# current working directory is project home
pipenv run python discomplex/data_prepare.py
```

### Perfomring predictions

```
pipenv run python discomplex/pairfolding.py 
```

### Loading Protein Sequences

Use the `load_fasta_sequences` function from the main script to load sequences from a FASTA file. The sequences are returned as a dictionary indexed by UniProt ID.

### Protein Folding Prediction

The `protein_folding_predict` function provides an interface to perform protein folding prediction using the ESMfold method. It requires sequences to be in list format.

### Analyzing Protein-Protein Interactions

The `find_interacting_residues` function analyzes a PDB file to find interacting regions in multi-chain protein complexes. It uses a cutoff distance to identify interactions.

### Running the Main Scripts

The repository contains two main scripts:

1. `data_prepare.py`: 
   ```
   pipenv run python disfold/data_prepare.py
   ```
2. `main_source_file_2.py`: (Replace with actual filename)
   - Used for...

Run these scripts with appropriate arguments to perform the analyses.

### Output

The scripts generate various output files, including PDB structures, interaction predictions, and intermediate data files.

## Additional Information

For more details on the individual functions and their parameters, refer to the docstrings provided within the scripts.
