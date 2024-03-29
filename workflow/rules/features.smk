##############################
#### Features subworkflow ####
##############################

"""
Results will be saved in results/features/
where single features are under results/features/single_vcf/<feature_name> 
and feature sets under results/features/feature_sets/<feature_set_name>
"""

from snakemake.utils import validate
import yaml

# general feature rule including all features


with open("config/features.yaml", "r") as file:
    # The FullLoader parameter handles the conversion from YAML
    # scalar values to Python the dictionary format
    features = yaml.load(file, Loader=yaml.FullLoader)

validate(features, schema="../schemas/features.schema.yaml")


# include multiple subwokflows for feature groups


include: "features/GCfeatures.smk"
include: "features/conservation.smk"
include: "features/encodeEpigenetics.smk"
include: "features/population.smk"
include: "features/fantom.smk"
include: "features/geneticVariation.smk"


# create a property file (before generating a sigle vcf file)
rule features_createPropertyFile:
    input:
        lambda wc: ancient(
            expand(
                "results/features/download/{file}/{{genomeBuild}}/{file}.{files}.{extension}.gz",
                file=features[wc.feature][wc.genomeBuild]["file"],
                files=features[wc.feature][wc.genomeBuild]["files"],
                extension=features[wc.feature][wc.genomeBuild]["type"],
            )
        ),
    output:
        config=(
            "results/features/single_vcf/{feature}/{genomeBuild}/{feature}.properties"
        ),
    params:
        file_type=lambda wc: features[wc.feature][wc.genomeBuild]["type"].split(".")[-1],
        column=(
            lambda wc: "column=%d" % features[wc.feature][wc.genomeBuild]["column"]
            if features[wc.feature][wc.genomeBuild]["type"].split(".")[-1] == "bed"
            else ""
        ),
        method=lambda wc: features[wc.feature][wc.genomeBuild]["method"],
        description=lambda wc: features[wc.feature]["description"],
    log:
        "results/logs/features/createPropertyFile.{feature}.{genomeBuild}.log",
    run:
        files = " \n".join(["file = %s" % file for file in input])
        shell(
            "echo -e 'name = {wildcards.feature} \n{files} \ntype = {params.file_type} \nmethod = {params.method} \ndescription = {params.description} \n{params.column}' > {output}"
        )


        # create single VCF file of features



rule features_createSingleFeatureVCF:
    input:
        config=ancient(
            "results/features/single_vcf/{feature}/{genomeBuild}/{feature}.properties"
        ),
        files=lambda wc: ancient(
            expand(
                "results/features/download/{file}/{{genomeBuild}}/{file}.{files}.{extension}.gz",
                file=features[wc.feature][wc.genomeBuild]["file"],
                files=features[wc.feature][wc.genomeBuild]["files"],
                extension=features[wc.feature][wc.genomeBuild]["type"],
            )
        ),
    output:
        temp=temp(
            "results/features/single_vcf/{feature}/{genomeBuild}/single/{feature}.temp.vcf.gz"
        ),
        vcf="results/features/single_vcf/{feature}/{genomeBuild}/single/{feature}.vcf.gz",
    params:
        mem="5g",
    conda:
        "../envs/jdk11.yaml"
    log:
        "results/logs/features/createSingleFeatureVCF.{feature}.{genomeBuild}.log",
    shell:
        """
        export LC_ALL=C;
        java -Xmx{params.mem} -jar workflow/bin/attributedb-cli-0.0.1-jar-with-dependencies.jar \
        vcf -p {input.config} --output {output.temp};
        (
            zcat {output.temp}| grep  "#";
            zcat {output.temp} | \
            grep  -v "#" | \
            sort -k1,1 -k2,2n
        ) | bgzip -c > {output.vcf} 2> {log}
        """


# index single feature vcf
rule features_indexSingleFeatureVCF:
    conda:
        "../envs/default.yml"
    input:
        "results/features/single_vcf/{feature}/{genomeBuild}/single/{feature}.vcf.gz",
    output:
        "results/features/single_vcf/{feature}/{genomeBuild}/single/{feature}.vcf.gz.tbi",
    log:
        "results/logs/features/indexSingleFeatureVCF.{feature}.{genomeBuild}.log",
    shell:
        """
        tabix {input} &> {log}
        """


# Average feature of defined position
rule features_average:
    conda:
        "../envs/default.yml"
    input:
        "results/features/single_vcf/{feature}/{genomeBuild}/single/{feature}.vcf.gz",
    output:
        "results/features/single_vcf/{feature}/{genomeBuild}/single/{feature}.avg.tsv.gz",
    log:
        "results/logs/features/average.{feature}.{genomeBuild}.log",
    shell:
        """
        zcat {input} | egrep -v "#" | awk -F'=' '{{ sum+=$2 }} END {{ print sum / NR }}' |  gzip -c > {output} 2> {log}
        """


# Create a feature set defined in config file
rule features_mergeSingleFeatureVCF:
    conda:
        "../envs/default.yml"
    input:
        files=lambda wc: expand(
            "results/features/single_vcf/{feature}/{genomeBuild}/single/{feature}.vcf.gz",
            feature=config["feature_sets"][wc.feature_set]["features"],
            genomeBuild=config["feature_sets"][wc.feature_set]["genome_build"],
        ),
        idx=lambda wc: expand(
            "results/features/single_vcf/{feature}/{genomeBuild}/single/{feature}.vcf.gz.tbi",
            feature=config["feature_sets"][wc.feature_set]["features"],
            genomeBuild=config["feature_sets"][wc.feature_set]["genome_build"],
        ),
    output:
        vcf="results/features/feature_sets/{feature_set}.vcf.gz",
        idx="results/features/feature_sets/{feature_set}.vcf.gz.tbi",
    log:
        "results/logs/features/mergeSingleFeatureVCF.{feature_set}.log",
    shell:
        """
        bcftools merge {input.files} | bgzip -c > {output.vcf} 2> {log};
        tabix {output.vcf} 2>> {log};
        """
