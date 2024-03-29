rule getConservation:
    output:
        "results/features/{genomeBuild}/{feature}/{chr}_{extension}.wigFix.gz",
    params:
        url=(
            lambda wildcards: config[wildcards.genomeBuild][wildcards.feature]["url"]
            + [wildcards.chr][0]
            + "."
            + [wildcards.extension][0]
            + ".wigFix.gz"
        ),
    shell:
        """
        echo {params.url}; curl {params.url} > {output}
        """


rule convertToBigWig:
    input:
        "results/features/{genomeBuild}/{feature}/{chr}.{extension}.wig.gz",
    output:
        "results/features/{genomeBuild}/{feature}/{chr}.{extension}.bw",
    shell:
        """
        bigWigToBedGraph {input} utils/{wildcards.genomeBuild}.chrom.sizes {output}
        """


rule getFantom19:
    output:
        "results/features/hg19/{feature}/{feature}.all.bed.gz",
    params:
        url=lambda wildcards: config["hg19"][wildcards.feature]["url"],
    shell:
        """
        curl {params.url}| sed '1d'  | bgzip > {output}
        """


'''        
rule get_DnaseClusteredHyp:
    output:
        "input/features/hg19/{feature}/{feature}.all.bed.gz"
    params:
        url=lambda wildcards: config['hg19'][wildcards.feature]['url']
    shell:
        """
        curl {params.url}  > {output}
        """        
'''
