rule getCytoband:
    output:
        "resources/{genomeBuild}/cytoBand.txt.gz",
    params:
        url="http://hgdownload.cse.ucsc.edu/goldenpath/{genomeBuild}/database/cytoBand.txt.gz",
    shell:
        """
        curl {params.url}  > {output}
        """

rule getChainFile:
    output:
        "resources/{genomeBuild}/{chain}.over.chain.gz",
    params:
        url=(
            lambda wc: "https://hgdownload.soe.ucsc.edu/gbdb/%s/liftOver/%s.over.chain.gz"
            % (wc.genomeBuild, wc.chain)
        ),
    shell:
        """
        curl {params.url}  > {output}
        """
