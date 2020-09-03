#!/bin/bash
{params.path} = input/features/hg38/DGVCount/
 for file in  input/features/hg38/DGVCount/*; do 
        if  [[ $file == *txt.gz ]]; then
            a=$(basename $file .txt.gz); 
            zcat $file | cut -f 2- | bgzip >input/features/hg38/DGVCount/$a.bed.gz; rm $file;
            elif  [[ $file == *txt ]]; then
                a=$(basename $file .txt); 
                cat $file | cut -f 2- | gzip > input/features/hg38/DGVCount/$a.bed.gz; rm $file;
                elif [[ $file == *wigFix.gz ]]; then
                    rename .wigFix.gz .wig.gz input/features/hg38/DGVCount/*.wigFix.gz;  
                    elif [[ $file == *gvf.gz ]]; then
                    a=$(basename $file .txt.gz); 
                    gvcf2bed -I $file -O {params.path}$a.bed; rm $file;
                  #  convert2bed --input=gvf < GRCh38.variant_call.all.gvf > GRCh38.variant_call.all.bed
                   #     a=$(basename $file .gvf.gz); 
          fi;  
done;

