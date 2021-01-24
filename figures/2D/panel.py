import pandas as pd
import numpy as np
import json

def count(row, compounds):
    d = dict()
    for c in compounds:
        d[c] = len(list(filter(lambda x: x >= 4.0, row[c])))
    return d 

df = pd.read_csv('../../peptides_compounds.tsv', sep='\t')
dup = df.copy()
dup.pop('description')
dup.pop('accession')
dup.pop('site')

df = df[~dup.duplicated()]

cmpds = df.columns[7:]
print(cmpds)

liganded_peptides = count(df, cmpds)
liganded_proteins = dict()
liganded_peptides_per_protein = dict()
competitors_proteins = dict()

for acc in set(df['accession']):    
    # counts for this protein
    subset = df[df['accession'] == acc]
    comps = 0
    for c in cmpds:
        if np.max(subset[c]) >= 4:
            c = c.split('_')[0]
            liganded_proteins[c] = liganded_proteins.get(c,0) + 1
            comps += 1
    competitors_proteins[str(comps)] = competitors_proteins.get(str(comps), 0) + 1

    peps = 0
    for idx, row in subset.iterrows():
        if len(list(filter(lambda x: x >= 4.0, row[cmpds]))) > 0:
            peps += 1
    liganded_peptides_per_protein[str(peps)] = liganded_peptides_per_protein.get(str(peps), 0) + 1

competitors_proteins_fmt = {
    '1': competitors_proteins['1'],
    '2': competitors_proteins['2'],
    '3': competitors_proteins['3'],
    '4': competitors_proteins['4'],
    '5': competitors_proteins['5'],
    '>=6': sum([competitors_proteins[x] for x in competitors_proteins.keys() if int(x) >= 6])
}

liganded_peptides_per_protein_fmt = {
    '1': liganded_peptides_per_protein['1'],
    '2': liganded_peptides_per_protein['2'],
    '3': liganded_peptides_per_protein['3'],
    '4': liganded_peptides_per_protein['4'],
    '5': liganded_peptides_per_protein['5'],
    '>=6': sum([liganded_peptides_per_protein[x] for x in liganded_peptides_per_protein.keys() if int(x) >= 6])
}

counts = ['1', '2', '3', '4', '5', '>=6']
df = pd.DataFrame({'Count': counts, 'Number of competitors': [competitors_proteins_fmt[x] for x in counts],
    'Liganded peptides per protein': [liganded_peptides_per_protein_fmt[x] for x in counts]
})
df.to_csv('figure 2d.csv')