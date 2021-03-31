

rule annotateFeatures:
    conda:
        "../envs/jdk11.yaml"
    input:
        feature_set="results/features/feature_sets/{feature_set}.vcf.gz",
        feature_set_idx="results/features/feature_sets/{feature_set}.vcf.gz.tbi",
        variants_file=lambda wc: getVariantsInput(wc.variant_set, "annotate"),
    output:
        "results/annotation/{variant_set}/{variant_set}.{feature_set}.tsv.gz",
    shell:
        """
        java -Xmx2g -jar workflow/bin/attributedb-cli-0.0.1-jar-with-dependencies.jar annotate-vcf \
        --annotation-vcf {input.feature_set} --file {input.variants_file} | bgzip -c > {output}
        """
