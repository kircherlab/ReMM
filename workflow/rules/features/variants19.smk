'''
rule get_feature_ISCApath:
    output:
        "input/features/{genomeBuild}/ISCApath/ISCApath.nstd37_nstd45_nstd75.bed.gz.vartmp"
    params:
        nstd75=lambda wildcards: "%sall_variants_for_nstd37.csv.gz " % (config[wildcards.genomeBuild]['ISCApath']['url']),
        nstd37=lambda wildcards: "%sall_variants_for_nstd45.csv.gz" % (config[wildcards.genomeBuild]['ISCApath']['url']),
        nstd45=lambda wildcards: "%sall_variants_for_nstd75.csv.gz " % (config[wildcards.genomeBuild]['ISCApath']['url']),
    shell:
        """       
    (curl {params.nstd75} | zcat; \
     curl {params.nstd37} | zcat; \
     curl {params.nstd45} | zcat) | egrep -v "^Study" | grep "GRCh37" | cut -f 1,8,12,13| awk -vIFS='\\t' -vOFS='\\t' '{{print "chr"$2,$3-1,$4,$1}}' | egrep -v "^chr\s" | sort -k 1,1 -k2,2n | bgzip -c > {output}
        """
all_variants_for_nstd37.csv.gz 
        
    
    

rule get_feature_ISCApath:
    output:
        "input/features/{genomeBuild}/ISCApath/ISCApath.nstd37_nstd45_nstd75.bed.gz.vartmp"
    params:
        nstd75=lambda wildcards: "%s%s.GRCh37.variant_call.tsv.gz" % (config[wildcards.genomeBuild]['ISCApath']['url'], 'nstd75'),
        nstd37=lambda wildcards: "%s%s.GRCh37.variant_call.tsv.gz" % (config[wildcards.genomeBuild]['ISCApath']['url'], 'nstd37'),
        nstd45=lambda wildcards: "%s%s.GRCh37.variant_call.tsv.gz" % (config[wildcards.genomeBuild]['ISCApath']['url'], 'nstd45'),
    shell:
        """       
    (curl {params.nstd75} | zcat; \
     curl {params.nstd37} | zcat; \
     curl {params.nstd45} | zcat) | egrep -v "^Study" | egrep -v "Pathogenic|pathogenic" | grep "GRCh37" | cut -f 1,8,12,13| awk -vIFS='\\t' -vOFS='\\t' '{{print "chr"$2,$3-1,$4,$1}}' | egrep -v "^chr\s" | sort -k 1,1 -k2,2n | bgzip -c > {output}
        """
'''


rule get_feature_dbVARCount19:
    output:
        "results/features/{genomeBuild}/dbVARCount/dbVARCount.all.bed.gz.vartmp",
    params:
        url=config["hg38"]["dbVARCount"]["url"],
    shell:
        """
        curl {params.url}| zcat  | egrep -v "#" | \
        egrep -v "benign|Benign" | sed  -E "s/[\|\;]/\\t/g" | \
        awk -vIFS='\\t' -vOFS='\\t' '{{print $1,$4-1,$5,$9}}' | \
        sort -k 1,1 -k2,2n | bgzip -c > {output}
        """


'''
rule get_feature_DGVCount19:
    output:
        "input/features/{genomeBuild}/DGVCount/DGVCount.all.bed.gz.vartmp"
    params:
        url=config['hg38']['DGVCount']['url']
    shell:
        """ 
        curl {params.url}  | egrep -v "^chr" | awk -vIFS='\\t' -vOFS='\\t' '{{print "chr"$2,$3-1,$4,$5,$6,$1}}' | sort -k 1,1 -k2,2n | awk '{{ for (i=1; i<NF;i++) {{ printf("%s\t",$i)}}; print $NF }}' | bgzip -c > {output} 
        """
'''


rule get_feature_Intervals19:
    input:
        "results/features/{genomeBuild}/{feature}/{feature}.{files}.bed.gz.vartmp",
    output:
        "results/features/{genomeBuild}/{feature}/{feature}.{files}.bed.gz",
    script:
        "../../scripts/createIntervals.py"
