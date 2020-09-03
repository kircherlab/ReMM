from ftplib import FTP
import pandas as pd
import os

#url = 'hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP30way/hg38.30way.phyloP/'
#filter_ ='chr[0-9]|chrX|chrY.*.phyloP30way.wigFix.gz'
#path = 'input/features/hg38/priPhyloP46way/rsync.txt'

CHROMS = ['chr'+str(s)  for s in list(range(1,23))+['Y','X']]

 
link = snakemake.params.url
filter_ = snakemake.params.filter
path = snakemake.output.path

ftp = link.split('/', 1)[0]
directory = '/' + link.split('/', 1)[1]
filter_ = '|'.join([filter_.replace('chr[0-24]',s) for s in CHROMS])

f = FTP(ftp)
f.login()
rsync = 'rsync -a -P rsync://' + ftp
f.cwd(directory)
files = f.nlst()

files =pd.Series(files)[pd.Series(files).str.contains(filter_)]

rsync = 'rsync -a -P rsync://' + ftp + directory + files +' ' + './' + path.replace('rsync.sh','') + ' ;'
p = os.getcwd()+'/' +path
with open(p, 'w') as f:
    f.write("\n".join(str(item) for item in rsync))
    
    
#os.system('bash ' + p)

