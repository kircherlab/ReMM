# input variants are different for each process step becaus ethey are flexible.
# this trys to collect the correct input. The following steps are possible only in that order:
# liftover > jannovar > bcftools > annotate
def getVariantsInput(variant_set, step):
    variant_set_config = config["variants"][variant_set]
    output = variant_set_config["file"]
    if step == "liftover":
        return output
    if "liftover" in variant_set_config:
        output = expand(
            "results/variants/{variant_set}/liftover/{variant_set}.vcf.gz",
            variant_set=variant_set,
        )
    if step == "jannovar":
        return output
    if "jannovar" in variant_set_config:
        output = expand(
            "results/variants/{variant_set}/jannovar/{variant_set}.vcf.gz",
            variant_set=variant_set,
        )
    if step == "bcftools":
        return output
    if "filters" in variant_set_config:
        if "bcftools" in variant_set_config['filters']:
            output = expand(
                "results/variants/{variant_set}/bcftools/{variant_set}.vcf.gz",
                variant_set=variant_set,
            )
    if step == "annotate":
        return output

    raise Exception("Unknown variant processing step %s" % step)


rule variants_liftover:
    conda:
        "../envs/ReMM.yaml"
    input:
        liftover_config=lambda wc: expand(
            "resources/{genomeBuild}/{liftover}.over.chain.gz",
            genomeBuild=config["variants"][wc.variant_set]["genome_build"],
            liftover=config["variants"][wc.variant_set]["liftover"],
        ),
        variants=lambda wc: getVariantsInput(wc.variant_set, "liftover"),
    output:
        f="results/variants/{variant_set}/liftover/{variant_set}.failed.bed.gz",
        o="results/variants/{variant_set}/liftover/{variant_set}.success.pos.gz",
    shell:
        """
        liftOver <(
            zcat {input.variants} | \
            grep -v "^#" | \
            awk 'BEGIN{{ OFS="\\t" }}{{ print $1,$2-1,$2,$1":"$2":"$4":"$5":"$8 }}'\
        ) {input.liftover_config} >( gzip -c > {output.o} ) >( gzip -c > {output.f});
        """


rule variant_liftover_filter:
    conda:
        "../envs/ReMM.yaml"
    input:
        variants=(
            "results/variants/{variant_set}/liftover/{variant_set}.success.pos.gz"
        ),
        ref_old=lambda wc: config["global_files"]["genome_builds"][
            config["variants"][wc.variant_set]["genome_build"]
        ]["reference"],
        ref_new=lambda wc: config["global_files"]["genome_builds"][
            "hg38"
            if config["variants"][wc.variant_set]["genome_build"] == "hg19"
            else "hg19"
        ]["reference"],
    output:
        "results/variants/{variant_set}/liftover/{variant_set}.vcf.gz",
    shell:
        """
        python workflow/scripts/filterLiftoverVariants.py \
        --input {input.variants} \
        --reference-old {input.ref_old} \
        --reference-new {input.ref_new} \
        --output >(bgzip -c > {output})
        """


rule variants_annotateJannovar:
    conda:
        "../envs/jannovar.yaml"
    input:
        variants=lambda wc: getVariantsInput(wc.variant_set, "jannovar"),
        database=(
            lambda wc: "resources/%s.ser"
            % config["variants"][wc.variant_set]["jannovar"]
        ),
    output:
        variants="results/variants/{variant_set}/jannovar/{variant_set}.vcf.gz",
        idx="results/variants/{variant_set}/jannovar/{variant_set}.vcf.gz.tbi",
    shell:
        """
        jannovar -Xmx4g annotate-vcf -d {input.database} -i {input.variants} -o {output.variants};
        tabix {output.variants};
        """


rule variants_filter_bcftools:
    conda:
        "../envs/ReMM.yaml"
    input:
        lambda wc: getVariantsInput(wc.variant_set, "bcftools"),
    output:
        "results/variants/{variant_set}/bcftools/{variant_set}.vcf.gz",
    params:
        bcftools_filter=lambda wc: config["variants"][wc.variant_set]["filters"][
            "bcftools"
        ]["filter"],
    shell:
        """
        bcftools view {params.bcftools_filter} {input})| bgzip -c > {output};
        tabix {output.o};
        """
