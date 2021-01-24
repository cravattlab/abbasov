import os
import zipfile
import pandas as pd 

z = zipfile.ZipFile('figures.zip', 'w')

for root, dirs, files in os.walk('.'):
    for name in files:
        if '.csv' in name:
            z.write('{}\\{}'.format(root,name), arcname=name)

z.write('processed_excel.xlsx', arcname='processed_excel.xlsx')