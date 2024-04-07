import os
import pandas as pd

RAW_DIR = 'data/raw'

FILENAMES_RAW = {
    'complexportal':'idmapping_ComplexPortal_slim.tsv',
    'pdb':'idmapping_PDB.tsv',
    'disprot':'idmapping_disprot.tsv'
}
PARAMS = {

}
tables = {

}

def process_data(raw_names:dict=FILENAMES_RAW ):
    for d in raw_names.keys():
        indir = os.path.join(RAW_DIR, d)
        assert os.path.exists(indir)
        infile = os.path.join(indir, raw_names[d])
        assert os.path.exists(indir)
        
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
    