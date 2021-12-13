rule features_getEncodeEpigenetics:
    output:
        "results/features/download/{encodeEpigenetic}/{genomeBuild}/{encodeEpigenetic}.{file}.bigWig",
    params:
        url=lambda wc: "%s%s.bigWig" % (
            features[wc.encodeEpigenetic][wc.genomeBuild]["url"],
            wc.file,
        ),
    wildcard_constraints:
        encodeEpigenetic="(EncH3K27Ac)|(EncH3K27Ac_v1_4)|(EncH3K4Me1)|(EncH3K4Me1_v1_4)|(EncH3K4Me3)|(EncH3K4Me3_v1_4)",
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
