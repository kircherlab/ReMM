#!/bin/bash


## Need: bedtools, attributedb in bin folder

featureFile=$1
N=$2
refG=$3
output=$4
outputTemp=$4.temp
outputTempVcf=$4.temp.vcf.gz




bedtools random -n 120000 -l 1 -seed 42 -g $refG |sort -k1,1 -k2,2n -k3,3n | bgzip -c >  $outputTemp;
(echo -e "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"; zcat $outputTemp | sort -k1,1 -k2,2n | grep -v chrM | awk -v 'OFS=\t' 'length($1) < 6 { print $0 }' | awk -v 'OFS=\t' '{print $1,$2,".","T","A",".","PASS","."}') | bgzip -c > $outputTempVcf;

tabix $outputTempVcf;
java -Xmx2g -jar bin/attributedb-cli-0.0.1-jar-with-dependencies.jar annotate-vcf --annotation-vcf $featureFile --file random/hg38random.vcf.gz | bgzip -c > $output

rm $outputTemp $outputTempVcf $outputTempVcf.tbi

