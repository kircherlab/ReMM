#!/bin/bash

#DOWNLIFTEN
#get the variant positions + predictions

cd /fast/groups/ag_kircher/ReMM/MA_Lusi/Snakemake/;

#originlal: 49468603 und 169482752
# 48965345 und 169764963


paste  -d "\t"  <(zcat output/features/annotated/hg38/SNVs.hg38.combined.txt.gz | cut -f 1-3) output/predictions/hg38/SNVs.hg38.predictions.txt | awk 'BEGIN{ OFS="\t" }{ print $2,$3,$4,$1}' | gzip -c > output/predictions/lifted/new/hg38.predictions.txt.gz;


#lifover requires a bed file so -1
zcat  output/predictions/lifted/new/hg38.predictions.txt.gz  | awk 'BEGIN{ OFS="\t"}{ print $1,$2-1,$2,$1":"$2":"$3":"$4 }' | gzip -c > output/predictions/lifted/new/hg38.positions.bed.gz;


# liftover
liftOver output/predictions/lifted/new/hg38.positions.bed.gz input/variants/hg19/hg38ToHg19.over.chain.gz >( gzip -c > output/predictions/lifted/new/hg38.predictions.lifted.bed.gz ) >( gzip -c > output/predictions/lifted/new/failed.bed.gz);

# after liftover add 1 to start position
# these two positions are read in wrongly (a separator error?) so that not the end of the interval but the beginning 

#zcat  output/predictions/lifted/new/hg38.predictions.lifted.bed.gz | awk 'BEGIN{ OFS="\t"}{ print $1,$2+1,$4}' | sed "s/\b169482751\b/169482752/g" | sed  "s/\b49468602\b/49468603/g" |gzip -c > output/predictions/lifted/new/hg38.predictions.lifted.txt.gz;

zcat  output/predictions/lifted/new/hg38.predictions.lifted.bed.gz | awk 'BEGIN{ OFS="\t"}{ print $1,$2+1,$4}' |gzip -c > output/predictions/lifted/new/hg38.predictions.lifted.txt.gz;


# sort lifted bed 
#zcat output/predictions/lifted/new/hg38.predictions.lifted.txt.gz | sed s/^chr//g | sort -k1,1 -k2,2n | gzip -c > output/predictions/lifted/new/hg38.predictions.lifted.sorted.txt.gz;

# join with ReMM score
awk 'BEGIN{ OFS="\t"} NR==FNR{ a[$1 FS $2];next} ("chr"$1 FS $2) in a {print "chr"$1, $2, $3} ' <(zcat output/predictions/lifted/new/hg38.predictions.lifted.txt.gz) <(zcat input/variants/hg19/ReMM.v0.3.1.tsv.gz) |gzip -c  > output/predictions/lifted/new/SNVs.lifted.remm.txt.gz;

