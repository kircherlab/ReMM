rule features_encodeEpigenetics_get:
    conda:
        "../../envs/default.yml"
    output:
        "results/features/download/{encodeEpigenetic}/{genomeBuild}/{encodeEpigenetic}.{file}.bigWig",
    params:
        url=lambda wc: "%s%s.bigWig" % (
            features[wc.encodeEpigenetic][wc.genomeBuild]["url"],
            wc.file,
        ),
    wildcard_constraints:
        encodeEpigenetic="(EncH3K27Ac)|(EncH3K27Ac_v1_4)|(EncH3K4Me1)|(EncH3K4Me1_v1_4)|(EncH3K4Me3)|(EncH3K4Me3_v1_4)",
    log:
        temp(
            "logs/features/encodeEpigenetics/get.{encodeEpigenetic}.{genomeBuild}.{file}.log"
        ),
    shell:
        """
        curl {params.url} > {output}
        """


rule features_encodeEpigenetics_mergeBigWigUsingMax:
    conda:
        "../../envs/ucsc_tools.yml"
    input:
        lambda wc: expand(
            "results/features/download/{{encodeEpigenetic}}/{{genomeBuild}}/{{encodeEpigenetic}}.{files}.bigWig",
            files=features[wc.encodeEpigenetic][wc.genomeBuild]["merge"],
        ),
    output:
        "results/features/download/{encodeEpigenetic}/{genomeBuild}/{encodeEpigenetic}.merged.bed.gz",
    log:
        temp(
            "logs/features/encodeEpigenetics/mergeBigWigUsingMax.{encodeEpigenetic}.{genomeBuild}.log"
        ),
    shell:
        """
        bigWigMerge -max {input} >(bgzip -c > {output}) &> {log}
        """


rule features_encodeEpigenetics_convertBigWigToBedGraph:
    conda:
        "../../envs/ucsc_tools.yml"
    input:
        "results/features/download/{file}/{genomeBuild}/{file}.{files}.bigWig",
    output:
        "results/features/download/{file}/{genomeBuild}/{file}.{files}.encode.bed.gz",
    log:
        temp(
            "logs/features/encodeEpigenetics/convertBigWigToBedGraph.{file}.{genomeBuild}.{files}.log"
        ),
    shell:
        """
        bigWigToBedGraph {input} >(bgzip -c > {output})
        """


rule features_encodeEpigenetics_convertBigWigToWig:
    conda:
        "../../envs/ucsc_tools.yml"
    input:
        "results/features/download/{file}/{genomeBuild}/{file}.{files}.bigWig",
    output:
        "results/features/download/{file}/{genomeBuild}/{file}.{files}.encode.wig.gz",
    log:
        temp(
            "logs/features/encodeEpigenetics/convertBigWigToWig.{file}.{genomeBuild}.{files}.log"
        ),
    shell:
        """
        bigWigToWig {input} >(bgzip -c > {output}) &> {log}
        """


rule features_encodeEpigenetics_getDnaseClusteredHyp:
    conda:
        "../../envs/default.yml"
    output:
        "results/features/download/DnaseClusteredHyp/{genomeBuild}/DnaseClusteredHyp.all.bed.gz",
    params:
        url=lambda wildcards: features["DnaseClusteredHyp"][wildcards.genomeBuild][
            "url"
        ],
    log:
        temp("logs/features/encodeEpigenetics/getDnaseClusteredHyp.{genomeBuild}.log"),
    shell:
        """
        curl {params.url} | zcat | cut -f 2- | bgzip -c > {output}
        """


rule features_encodeEpigenetics_getDnaseClusteredScore:
    conda:
        "../../envs/default.yml"
    output:
        "results/features/download/DnaseClusteredScore/{genomeBuild}/DnaseClusteredScore.all.bed.gz",
    params:
        url=lambda wc: features["DnaseClusteredScore"][wc.genomeBuild]["url"],
    log:
        temp("logs/features/encodeEpigenetics/getDnaseClusteredScore.{genomeBuild}.log"),
    shell:
        """
        curl {params.url} | zcat | cut -f 2- | bgzip -c > {output}
        """
