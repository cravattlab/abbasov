import pandas as pd
import numpy as np
import json 

codes = json.loads(open('manuscript_codes.json', 'r').read())
print(codes)

def try_match(s: str, dictionary): 
    if '_' not in s:
        return s
    name = s.split('_')[0]
    conc = s.split('_')[1]
    if name in dictionary:
        return dictionary[name] + '_' + conc
    else:
        for (k, v) in dictionary.items():
            if k in name:
                return v + '_' + conc
    return name

df = pd.read_csv('../../peptides_compounds.tsv', sep='\t')
print(df)

df = df.rename(lambda x: try_match(x, codes), axis="columns")
df['gene_name'] = df.apply(lambda x: x['description'].split(' ')[0], axis=1)
df.to_csv('renamed.csv')
print(df)