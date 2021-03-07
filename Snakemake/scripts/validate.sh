#!/bin/bash



# 



a=20284053
c=22756014 
b="15"

#echo $b, $a,$c;
#for file in nstd75.GRCh37.variant_call.tsv.gz nstd45.GRCh37.variant_call.tsv.gz nstd37.GRCh37.variant_call.tsv.gz; do
#    echo delete/$file
#    zcat delete/$file |  awk -v a=$a -v b=$b -v c=$c 'BEGIN{ OFS="\t"}  {if ($7==15) print $0;}'
    #zcat delete/$file |  awk -v a=$a -v b=$b -v c=$c 'NR==2 { print $8, $12, $13;}'


#done;
      
#zcat delete/nstd37.GRCh37.variant_call.tsv.gz | awk 'BEGIN{ OFS="\t"}  {if ($7==1) print $0}'

awk 'NR>1{a[$8]++} END{for(b in a) print b}' delete/nstd37.GRCh37.variant_call.tsv.gz 

# EncH3K27Ac EncH3K4Me1 EncH3K4Me3

#echo $b, $a;
#for l in EncH3K27Ac  EncH3K4Me1 EncH3K4Me3; do
#    echo $l
#    for file in input/features/hg38/$l/*wig.gz; do
#        zcat $file | awk -v a=$a -v b=$b 'BEGIN{ OFS="\t"}  {if (($1==b) && ($2<=a) && ($3>=a)) print $0;}'
#    done;
#        done;




#zcat output/features/combined/hg38/featureSet.hg38.vcf.gz | awk -v a=$a -v b=$b 'BEGIN{ OFS="\t"}  {if (($1==b) && ($2==a)) print $0;}' 


# from combined file for hg19
#zcat /fast/work/groups/ag_kircher/ReMM/ReMM/data/features/combined/ReMM/GRCh37/20180119/features.ReMM.20180119.vcf.gz  | awk -v a=$a -v c=$c 'BEGIN{ OFS="\t"}  {if (($1==b) && ($2==c)) print $0;}' 

#zcat input/variants/hg19/SNVs.hg19.negative.annotated.tsv.gz| awk 'BEGIN{ OFS="\t"}  {if (($1=="chr1") && ($2==10000069	)) print $0;}' 
