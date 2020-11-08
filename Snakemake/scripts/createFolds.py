import pandas as pd

file = str(snakemake.input.m)
df = pd.read_csv(file, compression = 'gzip',header = None, sep = '\t')
df = df.set_index([1])
df = df[~pd.isna(df.index)]

band = pd.read_csv(str(snakemake.input.c), compression = 'gzip', header = None, sep = '\t')
band  = band.set_index([3], drop = True)
band = band[~pd.isna(band.index)]

d = pd.concat([band,df[3]],axis = 1).reset_index()[[0,1,2,'index',3]]
d.to_csv(str(snakemake.output),  compression="gzip", sep = '\t', index = None,header=None)
