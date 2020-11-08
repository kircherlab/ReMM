import pandas as pd
import numpy as np
import pickle
path = '/fast/work/groups/ag_kircher/ReMM/MA_Lusi/Snakemake/'


data = str(snakemake.output)
input_file =str(snakemake.input.f)
cytoband =str(snakemake.input.cb)

df0 = pd.read_table(input_file)#[:7000000]
cytoband = pd.read_table(cytoband, header =None)

df = df0

print(df)
df['fold']=None
cytoband[[1,2]]=cytoband[[1,2]]+1
for i in range(len(cytoband)):
    ind = df[(df['CHR']==cytoband[0][i])&(cytoband[1][i]<=df['POSITION'])&(cytoband[2][i]>df['POSITION'])].index
    if len(ind>0):
        #print(len(ind))
        df.loc[ind,'fold']=cytoband[4][i]
              
df.to_csv(data, index = None, compression = 'gzip', sep = '\t', na_rep = 0)

