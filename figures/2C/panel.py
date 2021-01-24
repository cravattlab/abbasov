import pandas as pd 
import numpy as np 

df = pd.read_csv('../../peptides_chemotypes.tsv', sep='\t')
dup = df.copy()
dup.pop('description')
dup.pop('accession')
dup.pop('site')

df = df[~dup.duplicated()]
chemotypes = df.columns[7:]
print(chemotypes)

data = {}
for c in chemotypes:
    data[c] = dict(unique=0, shared=0)

for idx, row in df.iterrows():
    liganding_chemotypes = []
    for c in chemotypes:
        if row[c] >= 4.0:
            liganding_chemotypes.append(c)
    for c in liganding_chemotypes:
        if len(liganding_chemotypes) == 1:
            data[c]['unique'] = data[c]['unique'] + 1
        else:
            data[c]['shared'] = data[c]['shared'] + 1

df = pd.DataFrame(data).T 
df.to_csv('figure 2c.csv')