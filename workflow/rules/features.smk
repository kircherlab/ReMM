from snakemake.utils import validate
import yaml

# general feature rule including all features


with open("config/features.yaml", "r") as file:
    # The FullLoader parameter handles the conversion from YAML
    # scalar values to Python the dictionary format
    features = yaml.load(file, Loader=yaml.FullLoader)

validate(features, schema="../schemas/features.schema.yaml")


include: "features/GCfeatures.smk"
include: "features/conservation.smk"
include: "features/encodeEpigenetics.smk"
include: "features/population.smk"
include: "features/fantom.smk"
include: "features/geneticVariation.smk"


rule createPropertyFile:
    input:
        lambda wc: ancient(
            expand(
                "results/features/{file}/{{genomeBuild}}/{file}.{files}.{extension}.gz",
                file=features[wc.feature][wc.genomeBuild]["file"],
                files=features[wc.feature][wc.genomeBuild]["files"],
                extension=features[wc.feature][wc.genomeBuild]["type"],
            )
        ),
    output:
        config="results/features_single/{feature}/{genomeBuild}/{feature}.properties",
    params:
        file_type=lambda wc: features[wc.feature][wc.genomeBuild]["type"].split(".")[
            -1
        ],
        column=(
            lambda wc: "column=%d" % features[wc.feature][wc.genomeBuild]["column"]
            if features[wc.feature][wc.genomeBuild]["type"].split(".")[-1] == "bed"
            else ""
        ),
        method=lambda wc: features[wc.feature][wc.genomeBuild]["method"],
        description=lambda wc: features[wc.feature]["description"],
    run:
        files = " \n".join(["file = " + file for file in input])
        shell(
            """
        echo -e 'name = {wildcards.feature} \n{files} \ntype = {params.file_type}
        \nmethod = {params.method} \ndescription = {params.description} \n{params.column}' > {output}
        """
        )


rule createSingleFeatureVCF:
    input:
        config=ancient(
            "results/features_single/{feature}/{genomeBuild}/{feature}.properties"
        ),
        files=lambda wc: ancient(
            expand(
                "results/features/{file}/{{genomeBuild}}/{file}.{files}.{extension}.gz",
                file=features[wc.feature][wc.genomeBuild]["file"],
                files=features[wc.feature][wc.genomeBuild]["files"],
                extension=features[wc.feature][wc.genomeBuild]["type"],
            )
        ),
    output:
        temp=temp(
            "results/features_single/{feature}/{genomeBuild}/single/{feature}.temp.vcf.gz"
        ),
        vcf="results/features_single/{feature}/{genomeBuild}/single/{feature}.vcf.gz",
    params:
        mem="5g",
    conda:
        "../envs/jdk11.yaml"
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
        ) | bgzip -c > {output.vcf}
        """


rule indexSingleFeatreVCF:
    input:
        "results/features_single/{feature}/{genomeBuild}/single/{feature}.vcf.gz",
    output:
        "results/features_single/{feature}/{genomeBuild}/single/{feature}.vcf.gz.tbi",
    shell:
        """
        tabix {input}
        """


rule mergeSingleFeatureVCF:
    input:
        files=lambda wc: expand(
            "results/features_single/{feature}/{genomeBuild}/single/{feature}.vcf.gz",
            feature=config["feature_sets"][wc.feature_set]['features'],
            genomeBuild=config["feature_sets"][wc.feature_set]['genome_build']
        ),
    output:
        vcf="results/feature_sets/{feature_set}.vcf.gz",
        idx="results/feature_sets/{feature_set}.vcf.gz.tbi",
    shell:
        """
        bcftools merge {input.files} | bgzip -c > {output.vcf};
        tabix {output.vcf};
        """
