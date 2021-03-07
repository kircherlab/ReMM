import pandas as pd

df_neg = pd.read_csv(str(snakemake.input.na), compression = 'gzip', sep = '\t')
df_pos = pd.read_csv(str(snakemake.input.pa), skiprows  = 0, sep = '\t', compression = 'gzip')
df = pd.read_csv(str(snakemake.input.p), skiprows  = 7, sep = '\t', compression = 'gzip')

df_pos = df_pos[(df['REF'].str.len()==1)&(df['ALT'].str.len()==1)]

cols = df_neg.columns[~df_neg.columns.str.contains('fantom5MaxExpCell|fantom5MaxExpTissue|fantom5MaxUbiq|ERB91MotifFeatures|jaspar2018Cutoff400Max|jaspar2018Max')]
df_neg[cols].to_csv(str(snakemake.output.n), sep = '\t', compression = 'gzip', index = None,na_rep = 0)
df_pos[cols].to_csv(str(snakemake.output.p),sep = '\t', index = None,na_rep = 0,compression = 'gzip')

