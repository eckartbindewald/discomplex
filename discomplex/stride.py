
STRIDE_EXAMPLE="""\
REM  --------------------------------------------------------------------      
REM                                                                            
REM  STRIDE: Knowledge-based secondary structure assignment                    
REM  Please cite: D.Frishman & P.Argos, Proteins XX, XXX-XXX, 1995             
REM                                                                            
REM  Residue accessible surface area calculation                               
REM  Please cite: F.Eisenhaber & P.Argos, J.Comp.Chem. 14, 1272-1280, 1993     
REM               F.Eisenhaber et al., J.Comp.Chem., 1994, submitted           
REM                                                                            
REM  ------------------------ General information -----------------------      
REM                                                                            
HDR                                          01-JUN-22                         
CMP  MOL_ID: 1;                                                                
CMP   MOLECULE: GLUTAMINE SYNTHETASE INACTIVATING FACTOR IF7;                  
CMP   CHAIN: A                                                                 
SRC  MOL_ID: 1;                                                                
SRC   ORGANISM_SCIENTIFIC: SYNECHOCYSTIS SP. (STRAIN PCC 6714);                
SRC   ORGANISM_TAXID: 1147                                                     
REM                                                                            
REM  -------------------- Secondary structure summary -------------------      
REM                                                                            
CHN  /Users/eckart/spyder/discomplex/data-raw A                                
REM                                                                            
REM                .         .         .         .         .                   
SEQ  1    MSTQQQARALMMRHHQFIKNRQQSLLSRAAAEIGVQAEKDFWTTVQGKPQ   50              
STR         HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH    TTTT  EETTEE                    
REM                                                                            
REM                .                                                           
SEQ  51   SSFRETYDRSSASLS                                      65              
STR       HHHHHHH  TTTTT                                                       
REM                                                                            
REM                                                                            
REM                                                                            
LOC  AlphaHelix   THR     3 A      ILE     33 A                                
LOC  AlphaHelix   SER    51 A      TYR     57 A                                
LOC  Strand       THR    44 A      VAL     45 A                                
LOC  Strand       LYS    48 A      PRO     49 A                                
LOC  TurnI        GLU    38 A      PHE     41 A                                
LOC  TurnIV       THR    44 A      GLY     47 A                                
LOC  TurnI'       VAL    45 A      LYS     48 A                                
LOC  TurnI        SER    60 A      SER     63 A                                
LOC  TurnI        SER    61 A      LEU     64 A                                
REM                                                                            
REM  --------------- Detailed secondary structure assignment-------------      
REM                                                                            
REM  |---Residue---|    |--Structure--|   |-Phi-|   |-Psi-|  |-Area-|          
ASG  MET A    1    1    C          Coil    360.00    128.38     202.8          
ASG  SER A    2    2    C          Coil    -76.09    159.08      47.6          
ASG  THR A    3    3    H    AlphaHelix    -56.89    -35.61     104.8          
ASG  GLN A    4    4    H    AlphaHelix    -66.68    -40.29     133.0          
ASG  GLN A    5    5    H    AlphaHelix    -64.76    -39.89     112.5          
ASG  GLN A    6    6    H    AlphaHelix    -65.39    -45.47     117.7          
ASG  ALA A    7    7    H    AlphaHelix    -59.13    -43.97      49.6          
ASG  ARG A    8    8    H    AlphaHelix    -61.65    -44.24     190.4          
ASG  ALA A    9    9    H    AlphaHelix    -61.54    -41.97      55.2          
ASG  LEU A   10   10    H    AlphaHelix    -63.49    -44.82      91.2          
ASG  MET A   11   11    H    AlphaHelix    -61.21    -45.32     134.8          
ASG  MET A   12   12    H    AlphaHelix    -65.58    -40.46     120.9          
ASG  ARG A   13   13    H    AlphaHelix    -63.02    -45.44     137.0          
ASG  HIS A   14   14    H    AlphaHelix    -62.22    -42.26      53.7          
ASG  HIS A   15   15    H    AlphaHelix    -56.27    -47.14      48.4          
ASG  GLN A   16   16    H    AlphaHelix    -62.82    -37.33      84.6          
ASG  PHE A   17   17    H    AlphaHelix    -63.30    -46.02      84.8          
ASG  ILE A   18   18    H    AlphaHelix    -59.57    -45.06       1.8          
ASG  LYS A   19   19    H    AlphaHelix    -57.64    -52.14      58.9          
ASG  ASN A   20   20    H    AlphaHelix    -60.13    -41.01      76.1          
ASG  ARG A   21   21    H    AlphaHelix    -66.53    -43.99     110.8          
ASG  GLN A   22   22    H    AlphaHelix    -58.30    -49.56       4.6          
ASG  GLN A   23   23    H    AlphaHelix    -59.47    -44.04      45.0          
ASG  SER A   24   24    H    AlphaHelix    -61.03    -45.13      80.8          
ASG  LEU A   25   25    H    AlphaHelix    -66.97    -43.55      66.9          
ASG  LEU A   26   26    H    AlphaHelix    -67.34    -43.74      69.2          
ASG  SER A   27   27    H    AlphaHelix    -59.49    -40.31      55.6          
ASG  ARG A   28   28    H    AlphaHelix    -64.04    -51.12     169.3          
ASG  ALA A   29   29    H    AlphaHelix    -61.73    -39.53      44.3          
ASG  ALA A   30   30    H    AlphaHelix    -61.98    -44.32      13.4          
ASG  ALA A   31   31    H    AlphaHelix    -61.25    -37.71      73.6          
ASG  GLU A   32   32    H    AlphaHelix    -64.02    -33.31     132.6          
ASG  ILE A   33   33    H    AlphaHelix    -95.14      7.09     128.2          
ASG  GLY A   34   34    C          Coil     71.60     16.93      70.1          
ASG  VAL A   35   35    C          Coil    -93.20    131.23      72.4          
ASG  GLN A   36   36    C          Coil    -72.59    112.72     182.3          
ASG  ALA A   37   37    C          Coil    -67.67    109.50      34.4          
ASG  GLU A   38   38    T          Turn    -52.99    142.01     137.8          
ASG  LYS A   39   39    T          Turn    -58.16    -19.61     187.5          
ASG  ASP A   40   40    T          Turn   -112.60     25.94     130.3          
ASG  PHE A   41   41    T          Turn    -70.19     70.94     103.3          
ASG  TRP A   42   42    C          Coil   -139.94    118.56      46.6          
ASG  THR A   43   43    C          Coil    -99.68    115.81      95.7          
ASG  THR A   44   44    E        Strand    -94.28    142.02      65.6          
ASG  VAL A   45   45    E        Strand    -93.94    113.95      77.1          
ASG  GLN A   46   46    T          Turn     54.06     40.35     174.2          
ASG  GLY A   47   47    T          Turn     73.50      0.77      50.1          
ASG  LYS A   48   48    E        Strand   -102.42    150.79     134.1          
ASG  PRO A   49   49    E        Strand    -56.86    145.80      15.5          
ASG  GLN A   50   50    C          Coil    -53.39    139.02      79.1          
ASG  SER A   51   51    H    AlphaHelix    -59.97    -35.38      82.7          
ASG  SER A   52   52    H    AlphaHelix    -59.79    -31.51      56.3          
ASG  PHE A   53   53    H    AlphaHelix    -74.26    -23.70      37.1          
ASG  ARG A   54   54    H    AlphaHelix    -66.25    -33.86     131.3          
ASG  GLU A   55   55    H    AlphaHelix    -82.99    -38.27     163.7          
ASG  THR A   56   56    H    AlphaHelix    -93.18    -28.26      84.2          
ASG  TYR A   57   57    H    AlphaHelix   -109.41     -3.38      46.1          
ASG  ASP A   58   58    C          Coil    -65.51    158.25      74.6          
ASG  ARG A   59   59    C          Coil    -62.67    127.22     162.0          
ASG  SER A   60   60    T          Turn    -59.32    134.58      75.0          
ASG  SER A   61   61    T          Turn    -56.83    -20.62     110.2          
ASG  ALA A   62   62    T          Turn    -69.47    -10.76      85.1          
ASG  SER A   63   63    T          Turn    -96.95      8.02      74.3          
ASG  LEU A   64   64    T          Turn    -86.71     82.05     158.2          
ASG  SER A   65   65    C          Coil   -129.94    360.00     183.1
"""

import pandas as pd
from io import StringIO
import re
import os
import subprocess

STRUCTURE_CODE = {
    'C':0, # coil (unstructured)
    'T':1, # turn
    'G':2, # 3_10 helix,
    'I':3, # pi-helix
    'H':4, # alpha helix,
    'b':5, # "isolated bridge",
    'B':6, # "bridge"
    'E':7 # beta-sheet
}

def parse_stride(text:str) -> pd.DataFrame:
    """
    Parses output from protein secondary structure prediction
    program called STRIDE
    Returns content as dataframe
    """
    if not isinstance(text, str):
        raise ValueError("parse_stride should be called with str data type but received " + str(type(text)) + " : " + str(text))
    lines = text.split("\n")
    last_index = -1
    for i in range(len(lines)):
        lines[i] = lines[i].strip()
        if lines[i].startswith("REM"):
            last_index = i
    if last_index < 0 or last_index >= len(lines):
        raise ValueError("Could not find start of data section - possibly due to wrong format? Expecting output of program STRIDE")
    region = lines[last_index+1:]
    tsv_data = [re.sub(r' +', '\t', line) for line in region]
    # print("tsv data:")
    # print(tsv_data)
    tsv_string_io = StringIO('\n'.join(tsv_data))
    df = pd.read_csv(tsv_string_io, sep='\t', header=None)
    
    # REM  |---Residue---|    |--Structure--|   |-Phi-|   |-Psi-|  |-Area-|          
    # ASG  MET A    1    1    C          Coil    360.00    128.38     202.8          

    df.columns=['Slug','ResName', 'Chain', 'ResPdbId', 'ResCountId', 'Struct', 'StructName', 'Phi', 'Psi', 'Area']
    df['StructCode'] = 0
    for i in range(len(df)):
        df.at[i, 'StructCode'] = STRUCTURE_CODE[df.at[i,'Struct']] 
    df = df.drop(columns=['Slug'])

    # print(df)
    return df


def test_parse_stride(text=STRIDE_EXAMPLE):
    result = parse_stride(text)
    assert isinstance(result, pd.DataFrame)
    assert len(result) > 0
    print(result.columns)


def run_stride(pdbfile:str, stride_binary='stride', verbose=True) -> pd.DataFrame:
    """
    Runs STRIDE program on a PDB file and returns
    a data frame with the protein secondary structure
    """
    assert ',' not in pdbfile
    assert '(' not in pdbfile
    command = [stride_binary, pdbfile]
    joined_command = ' '.join(command)
    assert ',' not in joined_command
    assert '(' not in joined_command
    if verbose:
        print("Stride command:")
        print(" ".join(command))
    try:
        result_text = subprocess.check_output(command, text=True)
        result = parse_stride(result_text)
        return result
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
        return pd.DataFrame()  # Return empty dataframe on error

def run_all_tests():
    test_parse_stride()


if __name__ == "__main__":
    run_all_tests()