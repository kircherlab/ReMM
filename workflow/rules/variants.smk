############################
### Variants subworkflow ###
############################

"""
Results will be saved in results/variants/<variant_set>
"""


# definitions and functions
include: "variants_defs.smk"


# liftover variants
rule variants_liftover:
    conda:
        "../envs/ucsc_tools.yml"
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
    log:
        temp("results/logs/variants/liftover.{variant_set}.log"),
    shell:
        """
        liftOver <(
            zcat {input.variants} | \
            grep -v "^#" | \
            awk 'BEGIN{{ OFS="\\t" }}{{ print $1,$2-1,$2,$1":"$2":"$4":"$5":"$8 }}'\
        ) {input.liftover_config} >( gzip -c > {output.o} ) >( gzip -c > {output.f}) &> {log};
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
    log:
        temp("results/logs/variants/liftover_filter.{variant_set}.log"),
    shell:
        """
        python workflow/scripts/filterLiftoverVariants.py \
        --input {input.variants} \
        --reference-old {input.ref_old} \
        --reference-new {input.ref_new} \
        --output >(bgzip -c > {output.vcf}) &> {log};
        tabix {output.vcf} &>> {log};
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
    log:
        temp("results/logs/variants/annotateJannovar.{variant_set}.log"),
    shell:
        """
        jannovar -Xmx4g annotate-vcf -d {input.database} -i {input.variants} -o {output.variants} &> {log};
        tabix {output.variants} &>> {log};
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
    log:
        temp("results/logs/variants/filter_bcftools.{variant_set}.log"),
    shell:
        """
        (
            echo -e "##fileformat=VCFv4.1\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO";
            bcftools view -H {params.bcftools_filter} {input};
        ) | bgzip -c > {output.vcf} 2> {log};
        tabix {output.vcf} &>> {log};
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
    log:
        temp("results/logs/variants/filter_downsample.{variant_set}.log"),
    shell:
        """
        (
            zcat {input} | egrep "^#";
            zcat {input} | egrep -v "^#" | sort -R | awk 'NR <= {params.sample_size} {{ print }}' | sort -k1,1 -k2,2n;
        ) | bgzip -c > {output.vcf} 2> {log};
        tabix {output.vcf} &>> {log};
        """
