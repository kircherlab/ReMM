rule get_cpgIslandExt_file:
    output:
        "results/features/{genome}/cpgIslandExt/cpgIslandExt.all.bed.gz",
    params:
        url="https://hgdownload.soe.ucsc.edu/goldenPath/{genome}/database/cpgIslandExt.txt.gz",
    shell:
        """
        curl {params.url} | zcat | cut -f 2- | bgzip -c > {output}
        """


rule get_GCContent:
    input:
        lambda wc: config["global_files"]["genome_builds"][wc.genomeBuild]["reference"],
    output:
        "results/features/{genomeBuild}/GCContent/GCContent.all.bed.gz",
    params:
        width=lambda wc: config[wc.genomeBuild]["GCContent"]["window_size"],
    conda:
        "../../envs/GCContent.yml"
    shell:
        """
        ruby scripts/getGCContentGenomeWide.rb {input} {params.width} | \
        bgzip -c > {output}
        """
