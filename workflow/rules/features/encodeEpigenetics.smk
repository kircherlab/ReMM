rule getEncH3K27Ac:
    output:
        "results/features/{genomeBuild}/EncH3K27Ac/EncH3K27Ac.{file}.bigWig",
    params:
        url=lambda wildcards: "%s%s.bigWig" % (
            config[wildcards.genomeBuild]["EncH3K27Ac"]["url"],
            wildcards.file,
        ),
    shell:
        """
        curl {params.url} > {output}
        """


rule getEncH3K4Me1:
    output:
        "results/features/{genomeBuild}/EncH3K4Me1/EncH3K4Me1.{file}.bigWig",
    params:
        url=lambda wildcards: "%s%s.bigWig" % (
            config[wildcards.genomeBuild]["EncH3K4Me1"]["url"],
            wildcards.file,
        ),
    shell:
        """
        curl {params.url} > {output}
        """


rule getEncH3K4Me3:
    output:
        "results/features/{genomeBuild}/EncH3K4Me3/EncH3K4Me3.{file}.bigWig",
    params:
        url=lambda wildcards: "%s%s.bigWig" % (
            config[wildcards.genomeBuild]["EncH3K4Me3"]["url"],
            wildcards.file,
        ),
    shell:
        """
        curl {params.url} > {output}
        """


## ancient is there to not rerun, delete sometime
rule convertBigWigToBedGraph:
    input:
        "results/features/{genomeBuild}/{file}/{file}.{files}.bigWig",
    output:
        "results/features/{genomeBuild}/{file}/{file}.{files}.encode.bed.gz",
    shell:
        """
        bigWigToBedGraph {input} >(bgzip -c > {output})
        """


rule getDnaseClusteredHyp:
    output:
        "results/features/{genomeBuild}/DnaseClusteredHyp/DnaseClusteredHyp.all.bed.gz",
    params:
        url=lambda wildcards: config[wildcards.genomeBuild]["DnaseClusteredHyp"][
            "url"
        ],
    shell:
        """
        curl {params.url} | zcat | cut -f 2- | bgzip -c > {output}
        """


rule getDnaseClusteredScore:
    output:
        "results/features/{genomeBuild}/DnaseClusteredScore/DnaseClusteredScore.all.bed.gz",
    params:
        url=lambda wildcards: config[wildcards.genomeBuild]["DnaseClusteredScore"][
            "url"
        ],
    shell:
        """
        curl {params.url} | zcat | cut -f 2- | bgzip -c > {output}
        """
