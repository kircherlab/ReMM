#!/bin/bash
for file in  /fast/groups/ag_kircher/ReMM/MA_Lusi/Snakemake/input/features/hg38/mamPhyloXway/*; 
    do a=$(basename $file .bed.gz); zcat $file | cut -f 1-3,5 | bgzip -c  > _$a.bed.gz;
        done;
        
        
#rename _ '' *; 
zcat GRCh38_hg38_variants_2020-02-25.bed.gz |head -1  | tr '\t' '\n' | cat -n | grep "frequency"


#zcat hg38.phastCons17way.wig.gz | wig2bed | bgzip  > _hg38.phastCons17way.bed.gz


#rm -rf chr*


zcat d.bed.gz | cut -f 1-3,5 | bgzip -c  > _d.bed.gz.bed.gz;