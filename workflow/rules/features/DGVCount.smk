
rule getDGVCount:
    output:
        "results/features/{genomeBuild}/DGVCount/DGVCount.all.bed.gz.dgv",
    params:
        url=lambda wc: config[wc.genomeBuild]["DGVCount"]["url"],
    shell:
        """ 
        curl {params.url} | \
        awk -vIFS='\\t' -vOFS='\\t' '{{print "chr"$2,$3-1,$4,$5,$6,$1}}' | \
        sort -k 1,1 -k2,2n | \
        bgzip -c > {output} 
        """


#        """
#        curl {params.url}  | egrep -v "^chr" | awk -vIFS='\\t' -vOFS='\\t' '{{print "chr"$2,$3-1,$4,$5,$6,$1}}' | sort -k 1,1 -k2,2n | awk '{{ for (i=1; i<NF;i++) {{ printf("%s\t",$i)}}; print $NF }}' | bgzip -c > {output}
#       """


rule getDGVCountIntervals:
    input:
        "results/features/{genomeBuild}/DGVCount/DGVCount.all.bed.gz.dgv",
    output:
        "results/features/{genomeBuild}/DGVCount/DGVCount.{file}.bed.gz",
    script:
        "../../scripts/createIntervals.py"
