import pandas as pd 
import csv
from tqdm import tqdm
import os

def main(traitfile, speciesfile, outfile):
    os.makedirs(os.path.dirname(outfile), exist_ok=True)

    with open(traitfile, 'r') as tfile:
        dialect = csv.Sniffer().sniff(tfile.readline())
        
    trait_df = pd.read_csv(traitfile, delimiter=dialect.delimiter)
    species_df = pd.read_csv(speciesfile)
    
    specieslist = set(list(species_df.stack().values))
    for sp in tqdm(specieslist, desc="Matching taxa to traits"):
        trait_df.loc[trait_df.iloc[:, 0] == sp, "tree_name"] = sp
                
    trait_df.to_csv(outfile, index=False)