import pandas as pd
input_file = str(snakemake.input)
df = pd.read_csv(input_file,header = None,sep ='\t')
CHROMS = ['chr'+str(s)  for s in list(range(1,23))+['Y','X']]
for c in CHROMS:
    file = 'input/features/hg38/numTFBSConserved/numTFBSConserved.'+c+'.bed.gz.vartmp'
    d = df[df[0]==c]
    d.to_csv(file, header = None,index = None, sep='\t')