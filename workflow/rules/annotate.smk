##################################
#### annotation sub workflow  ####
##################################


"""
Results will be saved in results/annotation/<variant_set>

Annotates variants from the variants subwokflow with features from a feature set.
It replaces the NaN value with the default value given in the config.
Output is a TSV file with headers.
"""


def getFeaturesOfFeatureSet(feature_set):
    """
    Return the features of a feature set
    """
    return config["feature_sets"][feature_set]["features"]


def getFeatureDefaultValue(feature, genome_build):
    """
    return the default value of a feature (and its genome build)
    """
    return features[feature][genome_build]["default_value"]


# annotate variants with features
rule annotate_features:
    conda:
        "../envs/jdk11.yaml"
    input:
        feature_set="results/features/feature_sets/{feature_set}.vcf.gz",
        feature_set_idx="results/features/feature_sets/{feature_set}.vcf.gz.tbi",
        variants_file=lambda wc: getVariantsInput(wc.variant_set, "annotate"),
        variants_file_idx=lambda wc: getVariantsInput(
            wc.variant_set, "annotate", idx=True
        ),
    output:
        "results/annotation/{variant_set}/{variant_set}.{feature_set}.unsorted.tsv.gz",
    shell:
        """
        java -Xmx2g -jar workflow/bin/attributedb-cli-0.0.1-jar-with-dependencies.jar annotate-vcf \
        --annotation-vcf {input.feature_set} --file {input.variants_file} | bgzip -c > {output}
        """


# Sort features and replace NaN value with default of feature defined in config
rule annotate_sort_features:
    input:
        "results/annotation/{variant_set}/{variant_set}.{feature_set}.unsorted.tsv.gz",
    output:
        "results/annotation/{variant_set}/{variant_set}.{feature_set}.sorted.tsv.gz",
    params:
        features=lambda wc: " ".join(
            [
                "--feature %s %f"
                % (
                    feature,
                    getFeatureDefaultValue(
                        feature, getVariantSetGenomeBuild(wc.variant_set)
                    ),
                )
                for feature in getFeaturesOfFeatureSet(wc.feature_set)
            ]
        ),
    shell:
        """
        python workflow/scripts/sortAnnotationFile.py --input {input} --output {output} {params.features}
        """
