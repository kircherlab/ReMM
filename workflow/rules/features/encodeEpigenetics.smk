rule getEncH3K27Ac:
    output:
        "results/features/download/EncH3K27Ac/{genomeBuild}/EncH3K27Ac.{file}.bigWig",
    params:
        url=lambda wildcards: "%s%s.bigWig" % (
            features["EncH3K27Ac"][wildcards.genomeBuild]["url"],
            wildcards.file,
        ),
    shell:
        """
        curl {params.url} > {output}
        """


rule getEncH3K4Me1:
    output:
        "results/features/download/EncH3K4Me1/{genomeBuild}/EncH3K4Me1.{file}.bigWig",
    params:
        url=lambda wildcards: "%s%s.bigWig" % (
            features["EncH3K4Me1"][wildcards.genomeBuild]["url"],
            wildcards.file,
        ),
    shell:
        """
        curl {params.url} > {output}
        """


rule getEncH3K4Me3:
    output:
        "results/features/download/EncH3K4Me3/{genomeBuild}/EncH3K4Me3.{file}.bigWig",
    params:
        url=lambda wildcards: "%s%s.bigWig" % (
            features["EncH3K4Me3"][wildcards.genomeBuild]["url"],
            wildcards.file,
        ),
    shell:
        """
        curl {params.url} > {output}
        """


## ancient is there to not rerun, delete sometime
rule convertBigWigToBedGraph:
    input:
        "results/features/download/{file}/{genomeBuild}/{file}.{files}.bigWig",
    output:
        "results/features/download/{file}/{genomeBuild}/{file}.{files}.encode.bed.gz",
    shell:
        """
        bigWigToBedGraph {input} >(bgzip -c > {output})
        """


rule getDnaseClusteredHyp:
    output:
        "results/features/download/DnaseClusteredHyp/{genomeBuild}/DnaseClusteredHyp.all.bed.gz",
    params:
        url=lambda wildcards: features["DnaseClusteredHyp"][wildcards.genomeBuild][
            "url"
        ],
    shell:
        """
        curl {params.url} | zcat | cut -f 2- | bgzip -c > {output}
        """


rule getDnaseClusteredScore:
    output:
        "results/features/download/DnaseClusteredScore/{genomeBuild}/DnaseClusteredScore.all.bed.gz",
    params:
        url=lambda wc: features["DnaseClusteredScore"][wc.genomeBuild]["url"],
    shell:
        """
        curl {params.url} | zcat | cut -f 2- | bgzip -c > {output}
        """
