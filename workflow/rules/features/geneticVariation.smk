## use script getISCApath38 to get ISCApath from the file with all variants


rule getISCApath:
    conda:
        "../../envs/default.yml"
    output:
        "results/features/download/ISCApath/{genomeBuild}/ISCApath.allContigs.bed.gz",
    params:
        nstd75=lambda wc: "%s%s.%s.variant_call.tsv.gz"
        % (
            features["ISCApath"][wc.genomeBuild]["url"],
            "nstd75",
            "GRCh37" if wc.genomeBuild == "hg19" else "GRCh38",
        ),
        nstd37=lambda wc: "%s%s.%s.variant_call.tsv.gz"
        % (
            features["ISCApath"][wc.genomeBuild]["url"],
            "nstd37",
            "GRCh37" if wc.genomeBuild == "hg19" else "GRCh38",
        ),
        nstd45=lambda wc: "%s%s.%s.variant_call.tsv.gz"
        % (
            features["ISCApath"][wc.genomeBuild]["url"],
            "nstd45",
            "GRCh37" if wc.genomeBuild == "hg19" else "GRCh38",
        ),
    log:
        "results/logs/features/genetic_variation/getISCApath.{genomeBuild}.log",
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
        bgzip -c > {output} 2> {log}
        """


rule ISCApath_20211103:
    conda:
        "../../envs/default.yml"
    output:
        "results/features/download/ISCApath_20211103/{genomeBuild}/ISCApath_20211103.allContigs.bed.gz",
    params:
        nstd75=lambda wc: "%s%s.%s.variant_call.tsv.gz"
        % (
            features["ISCApath"][wc.genomeBuild]["url"],
            "nstd75",
            "GRCh37" if wc.genomeBuild == "hg19" else "GRCh38",
        ),
        nstd102=lambda wc: "%s%s.%s.variant_call.tsv.gz"
        % (
            features["ISCApath"][wc.genomeBuild]["url"],
            "nstd102",
            "GRCh37" if wc.genomeBuild == "hg19" else "GRCh38",
        ),
        nstd45=lambda wc: "%s%s.%s.variant_call.tsv.gz"
        % (
            features["ISCApath"][wc.genomeBuild]["url"],
            "nstd45",
            "GRCh37" if wc.genomeBuild == "hg19" else "GRCh38",
        ),
    log:
        "results/logs/features/genetic_variation/ISCApath_20211103.{genomeBuild}.log",
    shell:
        """       
        (
            curl {params.nstd75} | zcat; \
            curl {params.nstd102} | zcat; \
            curl {params.nstd45} | zcat
        )  | \
        cut -f 1,8,12,13 | \
        awk -vIFS='\\t' -vOFS='\\t' '{{print "chr"$2,$3-1,$4,$1}}' | \
        egrep -v "^chr\s" | \
        egrep -v "\s-1\s" | \
        sort -k1,1 -k2,2n | \
        bgzip -c > {output} 2> {log}
        """


rule downloadRefSeq2UCSC:
    conda:
        "../../envs/default.yml"
    output:
        "resources/{genomeBuild}/RefSeq2UCSC.txt",
    params:
        url=lambda wc: "https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh38_RefSeq2UCSC.txt"
        if wc.genomeBuild == "hg38"
        else "https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh37_NCBI2UCSC.txt",
    log:
        "results/logs/features/genetic_variation/downloadRefSeq2UCSC.{genomeBuild}.log",
    shell:
        """
        curl {params.url} | sed 's/\\r//' > {output} 2> {log}
        """


rule getDbVARCount:
    conda:
        "../../envs/default.yml"
    input:
        "resources/{genomeBuild}/RefSeq2UCSC.txt",
    output:
        "results/features/download/{dbVARCount}/{genomeBuild}/{dbVARCount}.allContigs.bed.gz",
    params:
        url=lambda wc: features[wc.dbVARCount][wc.genomeBuild]["url"],
    wildcard_constraints:
        dbVARCount="dbVARCount.*",
    log:
        "results/logs/features/genetic_variation/getDbVARCount.{genomeBuild}.{dbVARCount}.log",
    shell:
        """
        join -t $'\\t' <(cat {input} | sort -k1,1) \
        <(
            curl {params.url}| zcat  | egrep -v "#" | \
            egrep -v "benign|Benign" | sed  -E "s/[\|\;]/\\t/g" | \
            awk -vIFS='\\t' -vOFS='\\t' '{{print $1,$4-1,$5,$9}}' | \
            sort -k1,1
        )  | cut -f 2- | sort -k1,1 -k2,2n | bgzip -c > {output} 2> {log}
        """


rule getDGVCount:
    conda:
        "../../envs/default.yml"
    output:
        "results/features/download/{DGVCount}/{genomeBuild}/{DGVCount}.allContigs.bed.gz",
    params:
        url=lambda wc: features[wc.DGVCount][wc.genomeBuild]["url"],
    wildcard_constraints:
        DGVCount="DGVCount.*",
    log:
        "results/logs/features/genetic_variation/getDGVCount.{genomeBuild}.{DGVCount}.log",
    shell:
        """ 
        curl {params.url}   | \
        awk -vIFS='\\t' -vOFS='\\t' '{{print "chr"$2,$3-1,$4,$1,$5,$6}}' | \
        sort -k1,1 -k2,2n | bgzip -c > {output}  2> {log}
        """


rule getNumTFBSConserved:
    conda:
        "../../envs/default.yml"
    output:
        "results/features/download/numTFBSConserved/{genomeBuild}/numTFBSConserved.allContigs.bed.gz",
    params:
        url=lambda wc: features["numTFBSConserved"][wc.genomeBuild]["url"],
    log:
        "results/logs/features/genetic_variation/getNumTFBSConserved.{genomeBuild}.log",
    shell:
        """
        curl {params.url}  > {output} 2> {log}
        """


rule getEncRegTfbsClustered:
    conda:
        "../../envs/default.yml"
    output:
        "results/features/download/encRegTfbsClustered/{genomeBuild}/encRegTfbsClustered.allContigs.bed.gz",
    params:
        url=lambda wc: features["encRegTfbsClustered"][wc.genomeBuild]["url"],
    log:
        "results/logs/features/genetic_variation/getEncRegTfbsClustered.{genomeBuild}.log",
    shell:
        """
        curl {params.url}  > {output} 2> {log}
        """


rule getChromosomes:
    conda:
        "../../envs/default.yml"
    input:
        "results/features/download/{feature}/{genomeBuild}/{feature}.allContigs.bed.gz",
    output:
        temp(
            "results/features/download/{feature}/{genomeBuild}/{feature}.split_{chr}.bed.gz"
        ),
    params:
        chr="{chr}",
    wildcard_constraints:
        chr="|".join(["(chr%s)" % str(c) for c in list(range(1, 23)) + ["Y", "X"]]),
    log:
        "results/logs/features/genetic_variation/getChromosomes.{genomeBuild}.{feature}..{chr}.log",
    shell:
        """
        zcat {input} | grep -E "^{params.chr}\\s[0-9]+\\s[0-9]+\\s" | sort -k1,1 -k2,2n | gzip -c > {output} 2> {log}
        """


rule features_getIntervals:
    conda:
        "../../envs/default.yml"
    input:
        "results/features/download/{feature}/{genomeBuild}/{feature}.split_{chr}.bed.gz",
    output:
        "results/features/download/{feature}/{genomeBuild}/{feature}.{chr}.bed.gz",
    wildcard_constraints:
        feature="(DGVCount.*)|(numTFBSConserved)|(dbVARCount.*)|(ISCApath.*)|(encRegTfbsClustered)",
        chr="|".join(["(chr%s)" % str(c) for c in list(range(1, 23)) + ["Y", "X"]]),
    log:
        "results/logs/features/genetic_variation/getIntervals.{genomeBuild}.{feature}.{chr}log",
    script:
        "../../scripts/createIntervals.py"
