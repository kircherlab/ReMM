import pandas as pd
from pyfaidx import Fasta
import difflib


hg38 = Fasta('/fast/projects/cubit/current/static_data/reference/GRCh38/hs38/hs38.fa',
               one_based_attributes=False )
hg37 = Fasta('/fast/projects/cubit/current/static_data/reference/GRCh37/hs37/hs37.fa',
               one_based_attributes=False )

file_input = str(snakemake.input)
file_output = str(snakemake.output)

var = pd.read_csv(file_input, sep = '\t', header = None,compression="gzip")

ref_hg38 = []
diffs = []
for i in range(len(var)):
    ref = hg38[var[0][i]][var[1][i]:var[2][i]].seq
    n = 5
    ref_int_hg38 = hg38[var[0][i]][var[1][i]-n:var[2][i]+n].seq
    p = int(var[3][i].split(':')[1])
    ref_int_hg37 = hg37[var[0][i][3:]][p-1-n:p+n].seq
    diff= sum([i[0] != ' '  for i in difflib.ndiff(ref_int_hg38, ref_int_hg37)]) / 2
    diffs.append(1-diff/1001)
    ref_hg38.append(ref)
    
    
ref_hg37 = var[3].str.split(':', expand = True)[2]
vcf = pd.DataFrame()
vcf[[0,2]]=var[[0,3]]
vcf[1]=var[1]+1
# ID 
vcf[4]='.'
# or remain 
#vcf[4]=var[3].str.split(':', expand = True)[3]
vcf[3]=ref_hg38
vcf[[5,6]]='.'
vcf[7]=var[3].str.split(':', expand = True)[4]

vcf = vcf[[0,1,2,3,4,5,6,7]]
vcf = vcf[vcf[3]!=vcf[4]]

vcf = vcf[vcf[3].str.len()==1]
vcf = vcf[vcf[2].str.split(':',expand = True)[3].str.len()==1]

vcf.columns = ["#CHROM",'POS','ID','REF','ALT','QUAL','FILTER','INFO']

string = '##fileformat=VCFv4.1'
                    
cols = [string,'','','','','','','']
vcf = vcf.T.reset_index().T 
vcf.columns = cols



vcf.to_csv(file_output, sep = '\t', index = None) 


