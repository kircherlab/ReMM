import pandas as pd
import numpy as np
import pickle
path = '/fast/work/groups/ag_kircher/ReMM/MA_Lusi/Snakemake/'

folds = str(snakemake.output.f)
labels = str(snakemake.output.l)
data = str(snakemake.output.d)

input_file =str(snakemake.input.f)
cytoband =str(snakemake.input.cb)

df = pd.read_table(input_file,low_memory = False)#[:1000000]
cytoband = pd.read_table(cytoband, header =None)

#if 'hg38' or 'test' in input_file:
#    df = df0.drop(3, axis = 1)
#else:
#    df = df0
#df = df0.drop(3, axis = 1)
df['FOLD']=1
df['LABEL']=1
#cytoband[[1,2]]=cytoband[[1,2]]+1
#for i in range(len(cytoband)):
#    ind = df[(df[1]==cytoband[0][i])&(cytoband[1][i]<=df[2])&(cytoband[2][i]>df[2])].index
#    if len(ind>0):
        #print(len(ind))
#        df.loc[ind,'FOLD']=cytoband[4][i]

df['LABEL'].to_csv(labels, sep = '\t', header = None, index = None)
df['FOLD'].to_csv(folds, sep = '\t', header = None, index = None)
df.iloc[:,np.r_[3:28]].to_csv(data,sep = '\t', header = None, index = None,na_rep = 0)
    
    