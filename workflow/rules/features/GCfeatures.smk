rule features_GCfeatures_get_cpgIslandExt_file:
    conda:
        "../../envs/default.yaml"
    output:
        "results/features/download/cpgIslandExt/{genomeBuild}/cpgIslandExt.all.bed.gz",
    params:
        url="https://hgdownload.soe.ucsc.edu/goldenPath/{genomeBuild}/database/cpgIslandExt.txt.gz",
    log:
        temp("results/logs/features/GCfeatures/get_cpgIslandExt_file.{genomeBuild}.log"),
    shell:
        """
        curl {params.url} | zcat | cut -f 2- | bgzip -c > {output} 2> {log}
        """


rule features_GCfeatures_get_GCContent:
    conda:
        "../../envs/ruby.yaml"
    input:
        lambda wc: config["global_files"]["genome_builds"][wc.genomeBuild]["reference"],
    output:
        "results/features/download/GCContent/{genomeBuild}/GCContent.all.bed.gz",
    params:
        width=lambda wc: features["GCContent"][wc.genomeBuild]["window_size"],
    log:
        temp("results/logs/features/GCfeatures/get_GCContent.{genomeBuild}.log"),
    shell:
        """
        ruby workflow/scripts/getGCContentGenomeWide.rb {input} {params.width} | \
        bgzip -c > {output} 2> {log}
        """
