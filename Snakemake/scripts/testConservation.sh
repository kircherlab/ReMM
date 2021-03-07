#!/bin/bash

#CHR_38                chr1
#POSITION_38        7961859
#priPhyloP_38          0.45
#priPhastCons_38       0.45
#verPhyloP_38        -0.269
#verPhastCons_38      0.001
#mamPhyloP_38         0.857
#mamPhastCons_38      0.079
#mamPhastCons_19          0
#mamPhyloP_19        -0.032
#priPhastCons_19      0.135
#priPhyloP_19         0.305
#verPhastCons_19          0
#verPhyloP_19        -0.032
#Name: chr1-8021919


#hg19: chr11:17437007-17437008
#hg38: chr11:17415460-17415461


a=23464425
c=23464435

b="chr5"

echo $b, $a,$c;
for file in input/features/hg38/EncH3K4Me3/*.wig.gz; do
echo $file
    zcat $file | awk -v a=$a -v b=$b -v c=$c 'BEGIN{ OFS="\t"}  {if (($1==b) && ($2<=a) && ($3>=c)) print $0;}'
done;
        
#zcat output/features/combined/hg38/featureSet.hg38.vcf.gz| awk 'BEGIN{ OFS="\t"}  {if (($1=="chr5") && ($2<=23464240) && ($3>=23464250)) print $0;}'


