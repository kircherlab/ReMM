rule getChainFile:
    output:
        "resources/hg19ToHg38.over.chain.gz",
    params:
        url="https://hgdownload.soe.ucsc.edu/gbdb/{genomebuil}/liftOver/hg19ToHg38.over.chain.gz",
    shell:
        """
        curl {params.url}  > {output}
        """


## anpassen
rule getRefseqFile:
    output:
        "results/variants/{genomebuil}/data/hg38_refseq.ser",
    conda:
        "../../envs/jannovar.yml"
    params:
        path="results/variants/{genomebuil}/",
    shell:
        """
        java -jar jannovar-cli/target/jannovar-cli-0.33-SNAPSHOT.jar download -d hg38/refseq
        """
