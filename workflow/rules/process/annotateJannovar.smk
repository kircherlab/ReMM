rule annotateJannovarNegative:
    input:
        f="results/variants/hg38/SNVs.hg38.negative.vcf.gz",
        d="resources/hg38_refseq.ser",
    output:
        "results/variants/hg38/SNVs.hg38.negative.refseq.vcf.gz",
    conda:
        "../../envs/jannovar.yaml"
    shell:
        """
        jannovar -Xmx4g annotate-vcf -d {input.d} -i {input.f} -o {output};
        tabix {output}
        """

# TODO why implemented if not needed? also this can be used with a wildcard fpr positive or neggative
## not needed for the workflow
rule annotateJannovarPositive:
    input:
        f="resources/variants/hg38/SNVs.hg38.positive.vcf.gz",
        d="resources/data/hg38_refseq.ser",
    output:
        "results/variants/hg38/SNVs.hg38.positive.refseq.vcf.gz",
    conda:
        "../../envs/jannovar.yaml"
    shell:
        """
        jannovar annotate-vcf -d {input.d} -i {input.f} -o {output};
        tabix {output}
        """
