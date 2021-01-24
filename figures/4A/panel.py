import pandas as pd
import numpy as np

df = pd.read_csv('../../peptides_chemotypes.tsv', sep='\t')
dup = df.copy()
dup.pop('description')
dup.pop('accession')
dup.pop('site')

df = df[~dup.duplicated()]
yes = 0
no = 0

def any(arr, pred):
    p = False
    for a in arr:
        if pred(a):
            p = True 
    return p

for index, row in df.iterrows():
    if any(row[5:], (lambda x: x >= 4)):
        yes += 1
    else:
        no += 1
            

# df = pd.DataFrame({'Liganded peptides': yes, 'Unliganded peptides': no})
open('figure 4a left.csv', 'w').write('Liganded peptides, {}\nUnliganded peptides, {}\n'.format(yes, no))


df = pd.read_csv('../../proteins_chemotypes.tsv', sep='\t')
yes = 0
no = 0

for index, row in df.iterrows():
    if any(row[3:], (lambda x: x >= 4)):
        yes += 1
    else:
        no += 1


# df = pd.DataFrame({'Liganded peptides': yes, 'Unliganded peptides': no})
open('figure 4a right.csv', 'w').write(
    'Liganded proteins, {}\nUnliganded proteins, {}\n'.format(yes, no))
