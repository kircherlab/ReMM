rule getChainFile:
    output:
        "resources/hg19ToHg38.over.chain.gz",
    params:
        url="https://hgdownload.soe.ucsc.edu/gbdb/{genomebuil}/liftOver/hg19ToHg38.over.chain.gz",
    shell:
        """
        curl {params.url}  > {output}
        """
