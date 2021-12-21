rule get_cpgIslandExt_file:
    output:
        "results/features/download/cpgIslandExt/{genome}/cpgIslandExt.all.bed.gz",
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
        "results/features/download/GCContent/{genomeBuild}/GCContent.all.bed.gz",
    params:
        width=lambda wc: features["GCContent"][wc.genomeBuild]["window_size"],
    conda:
        "../../envs/ruby.yaml"
    shell:
        """
        ruby workflow/scripts/getGCContentGenomeWide.rb {input} {params.width} | \
        bgzip -c > {output}
        """
