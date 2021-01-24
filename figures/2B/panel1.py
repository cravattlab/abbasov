# Left panel of figure 2b: pie chart, distribution of competed lysines targeted
# by number of distinct fragments
import pandas as pd
import numpy as np

df = pd.read_csv('../../peptides_compounds.tsv', sep='\t')


dup = df.copy()
dup.pop('description')
dup.pop('accession')
dup.pop('site')

df = df[~dup.duplicated()]

cmpds = df.columns[6:]
print(cmpds)
# Count how many proteins are liganded in each chemotype or compound
data = dict()
for c in cmpds:
    c = c.split('_')[0]
    data[c] = {'liganded': 0}

# Keep track of which compounds (vlaues) target which proteins (key)
liganded_by = {}

for index, row in df.iterrows():
    for c in cmpds:
        if row[c] >= 4.0:
            c = c.split('_')[0]
            print(row['accession'], c)
            data[c]['liganded'] = data[c]['liganded'] + 1
            liganded_by[row['accession']] = liganded_by.get(
                row['accession'], []) + [c]


# Go through each protein, and count how many compounds have ratios >= 4
ligand_counts = [0 for i in range(10)]
for acc, ligands in liganded_by.items():
    ligand_counts[min(10, len(ligands)) - 1] += 1


df = pd.DataFrame({'Distinct chemotypes': [x+1 for x in range(10)], 'Number of liganded lysines': ligand_counts})
df.to_csv('figure 2b left.csv')
