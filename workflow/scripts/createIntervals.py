import pandas as pd
import warnings
warnings.filterwarnings('ignore')
#output_file = "input/features/hg38/numTFBSConserved/numTFBSConserved.chr1.bed.gz"
#input_file = "input/features/hg38/numTFBSConserved/numTFBSConserved.bed.gz.vartmp"

input_file = str(snakemake.input)
output_file = str(snakemake.output)
CHROMS = ['chr'+str(s)  for s in list(range(1,23))+['Y','X']]

all_intervals = pd.DataFrame()
all_intervals.to_csv(output_file)


all_intervals = pd.DataFrame()

if 'numTFBSConserved' in output_file:
    df_raw =pd.read_csv(input_file,header = None,sep = '\t',compression='gzip')
    df_raw[1] = df_raw[1].astype(dtype = str, errors = 'ignore')
    CHROMS = [output_file.split('.')[-4]]

    
if 'DGVCount' in output_file:
    df_raw =pd.read_csv(input_file,header = None, compression='gzip',sep = '\t')
    df_raw = df_raw[[0,1,2,5]]
    df_raw[1] = df_raw[1].astype(str)
    CHROMS = [output_file.split('.')[-3]]
    
    
else:
    df_raw =pd.read_csv(input_file,header = None, compression='gzip')[0].str.split('\t',expand =True)
df_raw = df_raw[df_raw[1].str.isnumeric()]

if input_file.split('/')[-1].strip()=='dbVARCount.all.bed.gz.vartmp':
    repl = pd.read_csv('utils/GRCh38_RefSeq2UCSC.txt', sep='\t', header= None)
    repl = dict(repl.values)
    
    df_raw[0] = df_raw[0].replace(repl)
    

for CHR in CHROMS:
    print(CHR)
    df0 = df_raw[df_raw[0]==CHR]
    df0 = df0[df0[1].str.isnumeric()]
    df0[[1,2]] = df0[[1,2]].apply(lambda x: x.astype(int))

    df = df0#.iloc[:1000]
    df =df.reset_index(drop = True)
    borders =pd.Series(pd.Series(df[1].to_list()+df[2].to_list()).unique()).sort_values().reset_index(drop=True)

    intervals =pd.DataFrame()

    start =[]
    end = []

    for i in range(len(borders)-1):
        start.append(borders[i])
        end.append(borders[i+1])

    intervals['start']=start
    intervals['end']=end
    intervals['variants']=''

    for raw in df.itertuples():
        s = raw[2]
        e = raw[3]
        v = raw[-1]
        dd = intervals.query('start >={}  & end <={}'.format(s,e)).loc[:,'variants']+v+ ', '
        intervals.loc[intervals.index.isin(dd.index.tolist()), 'variants']=dd

    intervals['count']=intervals['variants'].str.count(', ')

    intervals['chr'] = CHR

    intervals = intervals[['chr','start','end','count','variants']]
    intervals['variants']=intervals['variants'].str[:-2]

    all_intervals = all_intervals.append(intervals)

    all_intervals = all_intervals[all_intervals['count']!=0]
all_intervals.to_csv(output_file, sep = '\t', index = None, header = None, compression="gzip")