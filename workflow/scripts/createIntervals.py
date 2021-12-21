import pandas as pd
import warnings
warnings.filterwarnings('ignore')

input_file = str(snakemake.input)
output_file = str(snakemake.output)


all_intervals = pd.DataFrame()

df_raw =pd.read_csv(input_file,header = None,sep = '\t',compression='gzip')
df_raw[1] = df_raw[1].astype(dtype = str, errors = 'ignore')

for CHR in df_raw[0].unique():
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
    intervals['id']=''

    for raw in df.itertuples():
        start_tuple = raw[2]
        end_tuple = raw[3]
        id_tuple = raw[4]
        dd = intervals.query('start >={}  & end <={}'.format(start_tuple,end_tuple)).loc[:,'id']+id_tuple+ ','
        intervals.loc[intervals.index.isin(dd.index.tolist()), 'id']=dd

    intervals['count']=intervals['id'].str.count(',')

    intervals['chr'] = CHR

    intervals = intervals[['chr','start','end','count','id']]
    intervals['id']=intervals['id'].str[:-1]

    all_intervals = all_intervals.append(intervals)

    all_intervals = all_intervals[all_intervals['count']!=0]
all_intervals.to_csv(output_file, sep = '\t', index = None, header = None, compression="gzip")
