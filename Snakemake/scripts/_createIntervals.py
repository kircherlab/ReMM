import pandas as pd

#output_file = "input/features/hg38/DGVCount/DGVCount.all.bed.gz"
#input_file = "input/features/hg38/DGVCount/DGVCount.all.bed"

input_file = str(snakemake.input)
output_file = str(snakemake.output)

CHROMS = ['chr'+str(s)  for s in list(range(1,23))+['Y','X']]

all_intervals = pd.DataFrame()
if output_file.split('/')[-1]=='numTFBSConserved.all.bed.gz':
    df_raw =pd.read_csv(input_file,header = None,sep = '\t')
    df_raw[1] = df_raw[1].astype(dtype = str, errors = 'ignore')
else:
    df_raw =pd.read_csv(input_file,header = None, compression='gzip')[0].str.split('\t',expand =True)
df_raw = df_raw[df_raw[1].str.isnumeric()]

if input_file.split('/')[-1].strip()=='dbVARCount.all.bed.gz.vartmp':
    repl = pd.read_csv('utils/GRCh38_RefSeq2UCSC.txt', sep='\t', header= None)
    repl = dict(repl.values)
    print(repl)
    print(df_raw)
    df_raw[0] = df_raw[0].replace(repl)
    print(df_raw)
    #df_raw[0] = df_raw[0].str[:-3].replace(repl)


for CHR in CHROMS:
    #print(CHR)
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

    for i in range(len(df)):
        #if i % 20000 == 0:
            #print(i)
        s = df[1][i]
        e = df[2][i]
        intervals.loc[(intervals['start']>=s) & (intervals['end']<=e),'variants'] = intervals.loc[(intervals['start']>=s) & (intervals['end']<=e),'variants'] + df[df.columns[-1]][i]+', '
    intervals['count']=intervals['variants'].str.count(', ')

    intervals['chr'] = CHR

    intervals = intervals[['chr','start','end','count','variants']]
    intervals['variants']=intervals['variants'].str[:-2]

    all_intervals = all_intervals.append(intervals)

all_intervals = all_intervals[all_intervals['count']!=0]
all_intervals.to_csv(output_file, sep = '\t', index = None, header = None, compression="gzip")
