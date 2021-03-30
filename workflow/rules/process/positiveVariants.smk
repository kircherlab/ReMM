rule liftOverPositive:
    input:
        c="resources/hg19ToHg38.over.chain.gz",
        o="resources/SNVs.all.20160109.vcf.gz",
    output:
        t="results/variants/hg38/SNVs.all.20160109.bed.gz",
        f="results/variants/hg38/SNVs.hg38.positive.failed.bed.gz",
        o="results/variants/hg38/SNVs.hg38.positive.bed.gz.pos",
    shell:
        """
        zcat {input.o} | grep -v "^#" | \
        awk 'BEGIN{{ OFS="\\t" }}{{ print $1,$2-1,$2,$1":"$2":"$4":"$5":"$8 }}' | \
        gzip -c > {output.t};
        liftOver {output.t} {input.c} >( gzip -c > {output.o} ) >( gzip -c > {output.f});
        """


rule getPositiveVariants:
    input:
        "results/variants/{genomeBuild}/SNVs.{genomeBuild}.positive.bed.gz.pos",
    output:
        "results/variants/{genomeBuild}/SNVs.{genomeBuild}.positive.vcf",
    script:
        "../../scripts/filter.py"


'''
rule annotateJannovarPositive:
    input:
        f = "input/variants/hg38/SNVs.hg38.positive.vcf",
        d = "input/variants/hg38/data/hg38_refseq.ser"
    output:
        "input/variants/hg38/SNVs.hg38.positive.refseqannotated.vcf.gz"
    conda:
        "../../envs/jannovar.yml"
    shell:
        """
        jannovar -Xmx10G annotate-vcf -d {input.d} -i {input.f} -o {output};
        tabix {output};
        """
'''

# TODO Again a dummy rule???


rule dummyChangePositiveFileName:
    input:
        "results/variants/hg38/SNVs.hg38.positive.vcf",
    output:
        "results/variants/hg38/SNVs.hg38.positive.refseq.filtered.vcf.gz",
    shell:
        """
        bgzip < {input} > {output};
        tabix {output};
        """
