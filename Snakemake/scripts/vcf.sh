#!/bin/bash


export LC_ALL=C;
java -Xmx5g -jar bin/attributedb-cli-0.0.1-jar-with-dependencies.jar \
vcf -p input/features/hg38/PropertyFiles/priPhyloPXway.properties --output output/features/single/hg38/priPhyloPXway.vcf.gz;
#sort -k1,1 -k2,2n output/features/single/hg38/DnaseClusteredScore.vcf.gz  #output/features/single/hg38/DnaseClusteredScore.vcf.gz #output/features/single/hg38/priPhastConsXway.vcf.gz;
#tabix output/features/single/hg38/priPhastConsXway.vcf.gz;




