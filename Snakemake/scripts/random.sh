#!/bin/bash


##HG38
bedtools random -n 120000 -l 1 -seed 42 -g /fast/projects/cubit/current/static_data/reference/GRCh38/hs38/hs38.fa.genome |sort -k1,1 -k2,2n -k3,3n | bgzip -c > output/predictions/random/hg38random.txt.gz ;

(echo -e "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"; zcat output/predictions/random/hg38random.txt.gz | sort -k1,1 -k2,2n | grep -v chrM | awk -v 'OFS=\t' 'length($1) < 6 { print $0 }' | awk -v 'OFS=\t' '{print $1,$2,".","T","A",".","PASS","."}') | bgzip -c > output/predictions/random/hg38random.vcf.gz;

tabix output/predictions/random/hg38random.vcf.gz;

java -Xmx2g -jar bin/attributedb-cli-0.0.1-jar-with-dependencies.jar annotate-vcf --annotation-vcf output/features/combined/hg38/featureSet.hg38.vcf.gz --file output/predictions/random/hg38random.vcf.gz | bgzip -c > output/predictions/random/hg38random.features.vcf.gz;


##HG19
#lifover requires a bed file so -1
zcat  output/predictions/random/hg38random.txt.gz | sort -k1,1 -k2,2n | grep -v chrM | awk -v 'OFS=\t' 'length($1) < 6 { print $0 }'  | awk 'BEGIN{ OFS="\t"}{ print $1,$2-1,$2,$1":"$2":"$3 }' | gzip -c > output/predictions/random/hg19random.bed.gz;

liftOver output/predictions/random/hg19random.bed.gz input/variants/hg19/hg38ToHg19.over.chain.gz >( gzip -c > output/predictions/random/hg19random.lifted.bed.gz ) >( gzip -c > output/predictions/random/hg19random.failed.bed.gz);

zcat output/predictions/random/hg19random.lifted.bed.gz | awk 'BEGIN{ OFS="\t"}{ print $1,$3,$4 }' | sort -k1,1 -k2,2n  | bgzip -c > output/predictions/random/hg19random.lifted.txt.gz;


(echo -e "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";zcat output/predictions/random/hg19random.lifted.txt.gz | awk -v 'OFS=\t' '{print $1,$2,$3,"T","A",".","PASS","."}') | bgzip -c > output/predictions/random/hg19random.vcf.gz;

tabix output/predictions/random/hg19random.vcf.gz;

java -Xmx2g -jar bin/attributedb-cli-0.0.1-jar-with-dependencies.jar annotate-vcf --annotation-vcf /fast/groups/ag_kircher/ReMM/ReMM/data/features/combined/ReMM/GRCh37/20180119/features.ReMM.20180119.vcf.gz --file output/predictions/random/hg19random.vcf.gz| bgzip -c > output/predictions/random/hg19random.features.vcf.gz;



awk 'BEGIN{ OFS="\t"} NR==FNR{ a[$1 FS $2];next} ("chr"$1 FS $2) in a {print "chr"$1, $2, $3} ' <(zcat output/predictions/random/hg19random.lifted.txt.gz) <(zcat input/variants/hg19/ReMM.v0.3.1.tsv.gz) |gzip -c  > output/predictions/random/hg19random.remm.txt.gz ;



