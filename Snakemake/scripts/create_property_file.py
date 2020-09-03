print('skdj')
_files = snakemake.output
'''
with open('Utils/features_config.json') as f:
    config = json.load(f)


#_files = snakemake.output
#name = snakemake.feature
   
#method = 'upload' 
#description = 'Conservation primate phyloP46way downloaded from USCS' 
name = 'priPhyloPXway'
_files = ['input/features/hg38/priPhyloPXway/chr1_KI270707v1_random.phyloP30way.wigFix.gz',
#'input/features/hg38/priPhyloPXway/chr1_KI270706v1_random.phyloP30way.wigFix.gz']

import json
import datetime
from glob import glob
def get_properties(d):
    dd = []
    for key, value in d.items():
        dd.append(key + ' = ' + value)
    return " \n".join(dd[2:])


date = str(datetime.date.today())

files = " \n".join(['file = ' + file for file in _files])
property_list = get_properties(config['hg38'][name])

properties = " \n".join([name,date,property_list,files])

with open('proper', 'w') as outfile:
    for f in properties:
        outfile.write(f)
'''