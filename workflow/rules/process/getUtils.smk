rule process_getUtils_getCytoband:
    conda:
        "../../envs/default.yml"
    output:
        "resources/{genomeBuild}/cytoBand.txt.gz",
    params:
        url="http://hgdownload.cse.ucsc.edu/goldenpath/{genomeBuild}/database/cytoBand.txt.gz",
    log:
        temp("results/logs/process/getUtils/getCytoband.{genomeBuild}.log"),
    shell:
        """
        curl {params.url}  > {output} 2> {log}
        """


rule process_getUtils_getChainFile:
    conda:
        "../../envs/default.yml"
    output:
        "resources/{genomeBuild}/{chain}.over.chain.gz",
    params:
        url=(
            lambda wc: "https://hgdownload.soe.ucsc.edu/gbdb/%s/liftOver/%s.over.chain.gz"
            % (wc.genomeBuild, wc.chain)
        ),
    log:
        temp("results/logs/process/getUtils/getChainFile.{genomeBuild}.{chain}.log"),
    shell:
        """
        curl {params.url}  > {output} 2> {log}
        """
