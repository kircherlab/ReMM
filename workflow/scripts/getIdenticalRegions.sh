#!/bin/bash

region=$1
featureFile=$2
output=$3


(echo -e "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"; tabix $featureFile  $region | awk -v 'OFS=\t' '{print $1,$2,".","T","A",".","PASS","."}') |bgzip -c  > $output.temp.gz;

tabix -f $output.temp.gz;

java -Xmx2g -jar bin/attributedb-cli-0.0.1-jar-with-dependencies.jar annotate-vcf --annotation-vcf $featureFile --file $output.temp.gz  | bgzip -c  > $output

#rm $output.temp.gz $output.temp.gz.tbi
