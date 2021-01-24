import pandas as pd
import numpy as np
import json 

codes = json.loads(open('manuscript_codes.json', 'r').read())
print(codes)

def try_match(s: str, dictionary): 
    if s in dictionary:
        return dictionary[s]
    else:
        for (k, v) in dictionary.items():
            if k in s:
                return v
    return None

df = pd.read_csv('../../peptides_compounds.tsv', sep='\t')
print(df)
cmpds = list(map(lambda x: try_match(x.split('_')[0], codes), df.columns[7:]))
print(cmpds)
def f(row):
    r = pd.Series()
    r['accession'] = row['accession']
    r['gene_name'] = row['description'].split(' ')[0]
    r['description'] = row['description']
    r['sequence'] = row['sequence']
    r['site'] = row['site']
    r['max_ratio'] = row['max_ratio']
    r['average_ratio'] = row['average_ratio']
    hits = set(map(lambda x: x[1], filter(lambda x: x[0] >= 4, zip(row[7:], cmpds))))
    # map #2 (filter )
    r['hit_compounds'] = ','.join(hits)
    return r

df2 = df.apply(f, axis=1)
df2.to_csv('hit_compounds.csv')
print(df2)
