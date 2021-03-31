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
