def getFeaturesOfFeatuteSet(feature_set):
    return config["feature_sets"][feature_set]["features"]


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

# TODO add default value for feature
rule annotate_sort_features:
    input:
        "results/annotation/{variant_set}/{variant_set}.{feature_set}.unsorted.tsv.gz",
    output:
        "results/annotation/{variant_set}/{variant_set}.{feature_set}.sorted.tsv.gz",
    params:
        features=lambda wc: " ".join(
            ["--feature %s" % f for f in getFeaturesOfFeatuteSet(wc.feature_set)]
        ),
    shell:
        """
        python workflow/scripts/sortAnnotationFile.py --input {input} --output {output} {params.features}
        """
