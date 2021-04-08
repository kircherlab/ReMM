## use script getISCApath38 to get ISCApath from the file with all variants


rule getISCApath:
    output:
        "results/features/download/ISCApath/{genomeBuild}/ISCApath.allContigs.bed.gz",
    params:
        nstd75=lambda wc: "%s%s.GRCh38.variant_call.tsv.gz" % (
            features["ISCApath"][wc.genomeBuild]["url"],
            "nstd75",
        ),
        nstd37=lambda wc: "%s%s.GRCh38.variant_call.tsv.gz" % (
            features["ISCApath"][wc.genomeBuild]["url"],
            "nstd37",
        ),
        nstd45=lambda wc: "%s%s.GRCh38.variant_call.tsv.gz" % (
            features["ISCApath"][wc.genomeBuild]["url"],
            "nstd45",
        ),
    shell:
        """       
        (
            curl {params.nstd75} | zcat; \
            curl {params.nstd37} | zcat; \
            curl {params.nstd45} | zcat
        )  | cut -f 1,8,12,13 | \
        awk -vIFS='\\t' -vOFS='\\t' '{{print "chr"$2,$3-1,$4,$1}}' | \
        egrep -v "^chr\s" | \
        sort -k1,1 -k2,2n | \
        bgzip -c > {output}
        """


rule downloadRefSeq2UCSC:
    output:
        "resources/hg38/RefSeq2UCSC.txt",
    params:
        url="https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh38_RefSeq2UCSC.txt",
    shell:
        """
        curl {params.url} | sed 's/\\r//' > {output}
        """


rule getDbVARCount:
    input:
        "resources/{genomeBuild}/RefSeq2UCSC.txt",
    output:
        "results/features/download/dbVARCount/{genomeBuild}/dbVARCount.allContigs.bed.gz",
    params:
        url=features["dbVARCount"]["hg38"]["url"],
    shell:
        """
        join -t $'\\t' <(cat {input} | sort -k1,1) \
        <(
            curl {params.url}| zcat  | egrep -v "#" | \
            egrep -v "benign|Benign" | sed  -E "s/[\|\;]/\\t/g" | \
            awk -vIFS='\\t' -vOFS='\\t' '{{print $1,$4-1,$5,$9}}' | \
            sort -k1,1
        )  | cut -f 2- | sort -k1,1 -k2,2n | bgzip -c > {output}
        """


rule getDGVCount:
    output:
        "results/features/download/DGVCount/{genomeBuild}/DGVCount.allContigs.bed.gz",
    params:
        url=lambda wc: features["DGVCount"][wc.genomeBuild]["url"],
    shell:
        """ 
        curl {params.url}   | \
        awk -vIFS='\\t' -vOFS='\\t' '{{print "chr"$2,$3-1,$4,$1,$5,$6}}' | \
        sort -k1,1 -k2,2n | bgzip -c > {output} 
        """


rule getNumTFBSConserved:
    output:
        "results/features/download/numTFBSConserved/{genomeBuild}/numTFBSConserved.allContigs.bed.gz",
    params:
        url=lambda wc: features["numTFBSConserved"][wc.genomeBuild]["url"],
    shell:
        """
        curl {params.url}  > {output}
        """


rule getChromosomes:
    input:
        "results/features/download/{feature}/{genomeBuild}/{feature}.allContigs.bed.gz",
    output:
        temp("results/features/download/{feature}/{genomeBuild}/{feature}.split_{chr}.bed.gz"),
    params:
        chr="{chr}",
    wildcard_constraints:
        chr="|".join(["(chr%s)" % str(c) for c in list(range(1, 23)) + ["Y", "X"]]),
    shell:
        """
        zcat {input} | grep -E "^{params.chr}\\s" | sort -k1,1 -k2,2n | gzip -c > {output}
        """


rule features_getIntervals:
    input:
        "results/features/download/{feature}/{genomeBuild}/{feature}.split_{chr}.bed.gz",
    output:
        "results/features/download/{feature}/{genomeBuild}/{feature}.{chr}.bed.gz",
    wildcard_constraints:
        feature="(DGVCount)|(numTFBSConserved)|(dbVARCount)|(ISCApath)",
        chr="|".join(["(chr%s)" % str(c) for c in list(range(1, 23)) + ["Y", "X"]]),
    script:
        "../../scripts/createIntervals.py"
