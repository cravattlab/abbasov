import pandas as pd
import numpy as np
import json 

codes = json.loads(open('manuscript_codes.json', 'r').read())
# print(codes)

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

cells = set(['231', 'Ramos'])

counts = dict()

with open('expts', 'r') as f:
    # data = f.read()
    data = [x.split('/')[1].strip().replace('MIKA-', '').replace('.combined', '') for x in f.readlines()]

    for x in data:
        xs = x.split('_')
        cmpd = None 
        cell = None
        for val in xs:
            if val in codes:
                cmpd = codes[val]
            elif val in cells:
                cell = val 
            elif cmpd is None:
                for (k, v) in codes.items():
                    if k in val:
                        cmpd=v

            
        if cmpd is not None and cell is not None:
            tmp = counts.get(cmpd, dict())
            tmp[cell] = tmp.get(cell, 0) + 1
            counts[cmpd] = tmp
        else:
            print(x, cmpd, cell)

print(counts)

total = 0
with open('Tab2Experiments.csv', 'w') as f:
    f.write("Compound,231 datasets,Ramos datasets\n")
    for cmpd, v in counts.items():
        if cmpd in ['28q', '28r', '28s', '28t']:
            continue
        f.write('{},{},{}\n'.format(cmpd, v.get('231', 0), v.get('Ramos', 0)))