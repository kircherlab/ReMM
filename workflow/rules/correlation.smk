###################################
#### correlation sub workflow  ####
###################################


"""
Results will be saved in results/correlation/<correlation_set>

Correlation of features and predicted scores through genome builds.
"""


rule correlation_correlate_feature:
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
        "results/correlation/{correlation}/features/feature.{featureA}.{featureB}.correlate.tsv.gz",
    params:
        value_a="{featureA}",
        value_b="{featureB}",
    wrapper:
        getWrapper("evaluate/correlation")


def correlation_getCorrelateFeatures(correlation):
    features = [
        i.split("=")
        for i in config["correlation"][correlation]["correlate"]["features"]
    ]
    featureSetA = config["feature_sets"][
        config["correlation"][correlation]["A"]["feature_set"]
    ]["features"]
    featureSetB = config["feature_sets"][
        config["correlation"][correlation]["B"]["feature_set"]
    ]["features"]
    for feature in features:
        if feature[0] not in featureSetA:
            raise Exception(
                "Workflow correlation: cannot find feature in %s in feature set for set A!"
                % feature[0]
            )
        if feature[1] not in featureSetB:
            raise Exception(
                "Workflow correlation: cannot find feature in %s in feature set for set A!"
                % feature[0]
            )
    return features


rule correlation_combine_correlate_feature:
    input:
        lambda wc: expand(
            "results/correlation/{{correlation}}/features/feature.{featureA}.{featureB}.correlate.tsv.gz",
            zip,
            featureA=[f[0] for f in correlation_getCorrelateFeatures(wc.correlation)],
            featureB=[f[1] for f in correlation_getCorrelateFeatures(wc.correlation)],
        ),
    output:
        "results/correlation/{correlation}/feature.correlate.tsv.gz",
    params:
        columns=lambda wc: [
            "Feature A=%s" % f[0]
            for f in correlation_getCorrelateFeatures(wc.correlation)
        ],
    wrapper:
        getWrapper("file_manipulation/merge")
