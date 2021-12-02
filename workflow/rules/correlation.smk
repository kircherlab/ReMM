###################################
#### correlation sub workflow  ####
###################################


"""
Results will be saved in results/correlation/<correlation_set>

Correlation of features and predicted scores through genome builds.
"""


rule correlation_correlate_features:
    input:
        a=lambda wc: expand(
            "results/annotation/{variant_set}/{variant_set}.{feature_set}.sorted.tsv.gz",
            variant_set=config["correlation"][wc.correlation]["A"]["variants"],
            feature_set=config["correlation"][wc.correlation]["A"]["feature_set"],
        ),
        b=lambda wc: expand(
            "results/annotation/{variant_set}/{variant_set}.{feature_set}.sorted.tsv.gz",
            variant_set=config["correlation"][wc.correlation]["B"]["variants"],
            feature_set=config["correlation"][wc.correlation]["B"]["feature_set"],
        ),
    output:
        "results/correlation/{correlation}/feature.correlate.tsv.gz",
    params:
        value_a="GCContent",
        value_b="GCContent",
    wrapper:
        getWrapper("evaluate/correlate")
