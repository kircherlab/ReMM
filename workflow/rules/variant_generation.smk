##########################################
#### variant generation sub workflow  ####
##########################################


"""
Results will be saved in `results/variant_generation/<name>`

Generates variants but the same variant set for hg19 and hg38. Uses single positions from a bed file (config type regions) or generates random variants (config type random) of a defined length.

Final output is `results/variant_generation/<name>/hg19/<name>.vcf.gz"`
"""


# definitions and functions
include: "variant_generation_defs.smk"


# for variants make a bed file out of the vcf file
rule variant_generation_variantsToRegions:
    conda:
        "../envs/default.yml"
    input:
        lambda wc: config["variant_generation"][wc.name]["properties"]["file"],
    output:
        "results/variant_generation/{name}/hg38/variants_full.variants.bed.gz",
    log:
        "results/logs/variant_generation/variantsToRegions.{name}.log",
    shell:
        """
        zcat {input} | egrep -v "^#" | awk -v "OFS=\\t" '{{print $1,$2-1,$2,NR}}' | bgzip -c > {output} 2> {log}
        """


# just copy the bed file specified in the config
rule variant_generation_regionsToVariants:
    conda:
        "../envs/default.yml"
    input:
        lambda wc: config["variant_generation"][wc.name]["properties"]["file"],
    output:
        "results/variant_generation/{name}/hg38/regions_full.regions.bed.gz",
    log:
        "results/logs/variant_generation/regionsToVariants.{name}.log",
    shell:
        """
        cp {input} {output} &> {log}
        """


# generate random positions unsing bedtools random
rule variant_generation_random:
    conda:
        "../envs/default.yml"
    input:
        genome=lambda wc: config["global_files"]["genome_builds"]["hg38"]["genome"],
    output:
        "results/variant_generation/{name}/hg38/random_full.variants.bed.gz",
    params:
        seed=lambda wc: config["variant_generation"][wc.name]["properties"]["seed"],
        n=lambda wc: config["variant_generation"][wc.name]["properties"]["n"],
    log:
        "results/logs/variant_generation/random.{name}.log",
    shell:
        """
        bedtools random -n {params.n} -l 1 -seed {params.seed} -g {input.genome} | \
        egrep "^(chr[0-9]+\\s)|^(chrX\\s)|^(chrY\\s)" | \
        sort -k1,1 -k2,2n -k3,3n | \
        bgzip -c > {output} 2> {log};
        """


# liftover the variants to hg19
rule variant_generation_liftover:
    conda:
        "../envs/ucsc_tools.yml"
    input:
        liftover_config="resources/hg38/hg38ToHg19.over.chain.gz",
        variants=lambda wc: getRegionFileForVariantGeneration(wc.name),
    output:
        f="results/variant_generation/{name}/hg19/{name}.{regions_or_variants}.liftover.failed.bed.gz",
        o="results/variant_generation/{name}/hg19/{name}.{regions_or_variants}.liftover.success.pos.gz",
    log:
        "results/logs/variant_generation/liftover.{name}.{regions_or_variants}.log",
    shell:
        """
        liftOver <(
            zcat {input.variants} | \
            awk 'BEGIN{{ OFS="\\t" }}{{ print $1,$2,$3,$1":"$2":"$3":"$4}}'\q
        ) {input.liftover_config} >( gzip -c > {output.o} ) >( gzip -c > {output.f}) 2> {log}
        """


# split regions to window size 1
rule variant_generation_variantsFromRegions:
    """
    Create single position bed from regions. This has to be done for liftover hg19 
    as well as for the hg38 posituions in the idetifier. Otherwise variant_generation_getFinalHg38 will 
    have the region instead of positions.
    """
    conda:
        "../envs/default.yml"
    input:
        "results/variant_generation/{name}/hg19/{name}.regions.liftover.success.pos.gz",
    output:
        "results/variant_generation/{name}/hg19/{name}.variants.liftover.success.pos.gz",
    log:
        "results/logs/variant_generation/variantsFromRegions.{name}.log",
    shell:
        """
        zcat {input} | \
        egrep "^(chr[0-9]+\\s)|^(chrX\\s)|^(chrY\\s)" | \
        bedtools makewindows -w 1 -s 1 -i srcwinnum -b - | \
        sed 's/[:_]/\\t/g' | awk -v "OFS=\\t" '{{print $1,$2,$3,$4":"$5+($8-1)":"$5+$8":"$7"_"$8}}' | \
        bgzip -c > {output} 2> {log}
        """


# get the final hg38 variants (also defined for hg19)
rule variant_generation_getFinalHg38:
    conda:
        "../envs/default.yml"
    input:
        "results/variant_generation/{name}/hg19/{name}.variants.liftover.success.pos.gz",
    output:
        vcf="results/variant_generation/{name}/hg38/{name}.vcf.gz",
        idx="results/variant_generation/{name}/hg38/{name}.vcf.gz.tbi",
    log:
        "results/logs/variant_generation/getFinalHg38.{name}.log",
    shell:
        """
        (
            echo -e "##fileformat=VCFv4.1\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO";
            zcat {input} | \
            egrep "^(chr[0-9]+\\s)|^(chrX\\s)|^(chrY\\s)" | \
            cut -f 4 | \
            awk -F\: -v "OFS=\\t" '{{print $1,$2+1,"variant_"$4,"N",".",".","PASS","."}}' | \
            sort -k1,1 -k2,2n;
        ) | bgzip -c > {output.vcf} 2> {log};
        tabix {output.vcf} &>> {log};
        """


# get the final hg19 variants (also defined for hg38)
rule variant_generation_getFinalHg19:
    conda:
        "../envs/default.yml"
    input:
        "results/variant_generation/{name}/hg19/{name}.variants.liftover.success.pos.gz",
    output:
        vcf="results/variant_generation/{name}/hg19/{name}.vcf.gz",
        idx="results/variant_generation/{name}/hg19/{name}.vcf.gz.tbi",
    log:
        "results/logs/variant_generation/getFinalHg19.{name}.log",
    shell:
        """
        (
            echo -e "##fileformat=VCFv4.1\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO";
            zcat {input} | \
            egrep "^(chr[0-9]+\\s)|^(chrX\\s)|^(chrY\\s)" | \
            sed 's/:/\\t/g' | \
            awk -v "OFS=\\t" '{{print $1,$2+1,"variant_"$7,"N",".",".","PASS","."}}' | \
            sort -k1,1 -k2,2n;
        ) | bgzip -c > {output.vcf} 2> {log};
        tabix {output.vcf} &>> {log};
        """
