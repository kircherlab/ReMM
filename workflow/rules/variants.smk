#### Variants subworkflow ####

"""
Results will be saved in results/variants/
"""


def getVariantSetGenomeBuild(variant_set):
    """
    Getter of genomeBuild taking liftover into account.
    """
    variantSet_conf = config["variants"][variant_set]
    genomeBuild = variantSet_conf["genome_build"]
    switcher = {"hg38": "hg19", "hg19": "hg38"}
    if "processing" in variantSet_conf and "liftover" in variantSet_conf["processing"]:
        return switcher[genomeBuild]
    else:
        return genomeBuild


def getVariantsInput(variant_set, step, idx=False):
    """
    input variants are different for each process step becaus ethey are flexible.
    this trys to collect the correct input. The following steps are possible only in that order:
    liftover > jannovar > filters (bcftools > downsample) > annotate
    """

    add = ".tbi" if idx else ""
    variant_set_config = config["variants"][variant_set]
    if variant_set_config["type"] == "file":
        output = expand(
            "{file}{add}", file=variant_set_config["properties"]["file"], add=add
        )
    elif variant_set_config["type"] == "generation":
        output = expand(
            "results/variant_generation/{name}/{genomeBuild}/{name}.vcf.gz{add}",
            name=variant_set_config["properties"]["name"],
            genomeBuild=variant_set_config["genome_build"],
            add=add,
        )
    else:
        raise Exception("Unknown variant type %s" % variant_set_config["type"])
    if step == "liftover":
        return output
    # liftover
    if (
        "processing" in variant_set_config
        and "liftover" in variant_set_config["processing"]
    ):
        output = expand(
            "results/variants/{variant_set}/liftover/{variant_set}.vcf.gz{add}",
            variant_set=variant_set,
            add=add,
        )
    # jannovar
    if step == "jannovar":
        return output
    if (
        "processing" in variant_set_config
        and "jannovar" in variant_set_config["processing"]
    ):
        output = expand(
            "results/variants/{variant_set}/jannovar/{variant_set}.vcf.gz{add}",
            variant_set=variant_set,
            add=add,
        )
    # bcftools
    if step == "bcftools":
        return output
    if (
        "processing" in variant_set_config
        and "filters" in variant_set_config["processing"]
    ):
        if "bcftools" in variant_set_config["processing"]["filters"]:
            output = expand(
                "results/variants/{variant_set}/bcftools/{variant_set}.vcf.gz{add}",
                variant_set=variant_set,
                add=add,
            )
    # downsample
    if step == "downsample":
        return output
    if (
        "processing" in variant_set_config
        and "filters" in variant_set_config["processing"]
    ):
        if "downsample" in variant_set_config["processing"]["filters"]:
            output = expand(
                "results/variants/{variant_set}/downsample/{variant_set}.vcf.gz{add}",
                variant_set=variant_set,
                add=add,
            )
    if step == "annotate":
        return output

    raise Exception("Unknown variant processing step %s" % step)


# liftover variants
rule variants_liftover:
    conda:
        "../envs/ReMM.yaml"
    input:
        liftover_config=lambda wc: expand(
            "resources/{genomeBuild}/{liftover}.over.chain.gz",
            genomeBuild=config["variants"][wc.variant_set]["genome_build"],
            liftover=config["variants"][wc.variant_set]["processing"]["liftover"],
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


# filter lifover variants. check if reference is still the same
rule variants_liftover_filter:
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
        vcf="results/variants/{variant_set}/liftover/{variant_set}.vcf.gz",
        idx="results/variants/{variant_set}/liftover/{variant_set}.vcf.gz.tbi",
    shell:
        """
        python workflow/scripts/filterLiftoverVariants.py \
        --input {input.variants} \
        --reference-old {input.ref_old} \
        --reference-new {input.ref_new} \
        --output >(bgzip -c > {output.vcf});
        tabix {output.vcf};
        """


# annotate variants with jannovar
rule variants_annotateJannovar:
    conda:
        "../envs/jannovar.yaml"
    input:
        variants=lambda wc: getVariantsInput(wc.variant_set, "jannovar"),
        database=(
            lambda wc: "resources/%s.ser"
            % config["variants"][wc.variant_set]["processing"]["jannovar"]
        ),
    output:
        variants="results/variants/{variant_set}/jannovar/{variant_set}.vcf.gz",
        idx="results/variants/{variant_set}/jannovar/{variant_set}.vcf.gz.tbi",
    shell:
        """
        jannovar -Xmx4g annotate-vcf -d {input.database} -i {input.variants} -o {output.variants};
        tabix {output.variants};
        """


# filter variants with a bcftools filter command set in the config file
rule variants_filter_bcftools:
    conda:
        "../envs/ReMM.yaml"
    input:
        lambda wc: getVariantsInput(wc.variant_set, "bcftools"),
    output:
        vcf="results/variants/{variant_set}/bcftools/{variant_set}.vcf.gz",
        idx="results/variants/{variant_set}/bcftools/{variant_set}.vcf.gz.tbi",
    params:
        bcftools_filter=lambda wc: config["variants"][wc.variant_set]["processing"][
            "filters"
        ]["bcftools"]["filter"],
    shell:
        """
        (
            echo -e "##fileformat=VCFv4.1\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO";
            bcftools view -H {params.bcftools_filter} {input};
        ) | bgzip -c > {output.vcf};
        tabix {output.vcf};
        """


# filter variants with a bcftools filter command set in the config file
rule variants_filter_downsample:
    conda:
        "../envs/ReMM.yaml"
    input:
        lambda wc: getVariantsInput(wc.variant_set, "downsample"),
    output:
        vcf="results/variants/{variant_set}/downsample/{variant_set}.vcf.gz",
        idx="results/variants/{variant_set}/downsample/{variant_set}.vcf.gz.tbi",
    params:
        sample_size=lambda wc: config["variants"][wc.variant_set]["processing"][
            "filters"
        ]["downsample"],
    shell:
        """
        (
            zcat {input} | egrep "^#";
            zcat {input} | egrep -v "^#" | sort -R | awk 'NR <= {params.sample_size} {{ print }}' | sort -k1,1 -k2,2n;
        ) | bgzip -c > {output.vcf};
        tabix {output.vcf};
        """
