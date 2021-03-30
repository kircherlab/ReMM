rule getCytoband:
    output:
        "resources/cytoBand.{genomeBuild}.txt.gz",
    params:
        url="http://hgdownload.cse.ucsc.edu/goldenpath/{genomeBuild}/database/cytoBand.txt.gz",
    shell:
        """
        curl {params.url}  > {output}
        """


rule getFolds:
    input:
        m="resources/ReMM.20171122.partition.mapping.tsv.gz",
        c="resources/cytoBand.{genomeBuild}.txt.gz",
    output:
        "results/folds/{genomeBuild}/folds.{genomeBuild}.txt.gz",
    script:
        "../../scripts/createFolds.py"
