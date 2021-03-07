import pandas as pd


remm = pd.read_csv('input/variants/hg19/ReMM.v0.3.1.tsv.gz',sep = '\t',skiprows=2,header = None,compression = 'gzip',low_memory=False )
remm['ind']='chr'+remm[0].astype(str)+'-'+remm[1].astype(str)


var = pd.read_csv('output/features/annotated/hg19/SNVs.hg19.combined.txt.gz',sep = '\t', compression = 'gzip',header = None, usecols = [0,1,2])
var['ind']=var[1].astype(str)+'-'+var[2].astype(str)

file = 'output/predictions/hg19/SNVs.hg19.predictions.txt'
df = pd.read_csv(file, sep = '\t',header = None)

var['prediction'] = df[1]

f = var[[0,1,2,'prediction']].join(remm[2],how = 'inner',lsuffix = '_var',rsuffix = '_remm')
print('sjaks')
f.to_csv('output/predictions/hg19/remm_with_ps_predictions.csv', sep='\t',header=None, index=None)

