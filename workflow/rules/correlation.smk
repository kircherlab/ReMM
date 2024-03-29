###################################
#### correlation sub workflow  ####
###################################


"""
Results will be saved in `results/correlation/<correlation_set>`

Correlation of features and predicted scores through genome builds.

Final output is `results/correlation/<correlation>/feature.correlate.tsv.gz`
"""


include: "correlation_defs.smk"


### SCORES ###


# join two prediction files by the variant ID
rule correlation_scoreJoin:
    input:
        left=lambda wc: expand(
            "results/predictions/{training}/{variant_set}/predictions/predictions.with_IDs.tsv.gz",
            training=config["correlation"][wc.correlation]["A"]["training"],
            variant_set=config["correlation"][wc.correlation]["A"]["variants"],
        ),
        right=lambda wc: expand(
            "results/predictions/{training}/{variant_set}/predictions/predictions.with_IDs.tsv.gz",
            training=config["correlation"][wc.correlation]["B"]["training"],
            variant_set=config["correlation"][wc.correlation]["B"]["variants"],
        ),
    output:
        "results/correlation/{correlation}/prediction/join.tsv.gz",
    params:
        how="inner",
        suffixes="_A _B",
        left_on="ID",
        right_on="ID",
    log:
        temp("results/logs/correlation/correlation_scoreJoin.{correlation}.log"),
    wrapper:
        getWrapper("file_manipulation/merge")


# correlate the predicted score from both predictions
rule correlation_correlate_score:
    input:
        a="results/correlation/{correlation}/prediction/join.tsv.gz",
    output:
        "results/correlation/{correlation}/score.correlate.tsv.gz",
    params:
        value_a="SCORE_A",
        value_b="SCORE_B",
    log:
        temp(
            "results/logs/correlation/correlation_featcorrelation_correlate_scoreureJoin.{correlation}.log"
        ),
    wrapper:
        getWrapper("evaluate/correlation")


## FEATURES ###


# join two annotation files on the variant ID
rule correlation_featureJoin:
    input:
        left=lambda wc: expand(
            "results/annotation/{variant_set}/{variant_set}.{feature_set}.{missing_value}.sorted.tsv.gz",
            variant_set=config["correlation"][wc.correlation]["A"]["variants"],
            feature_set=config["correlation"][wc.correlation]["A"]["feature_set"],
            missing_value=config["correlation"][wc.correlation]["B"]["missing_value"],
        ),
        right=lambda wc: expand(
            "results/annotation/{variant_set}/{variant_set}.{feature_set}.{missing_value}.sorted.tsv.gz",
            variant_set=config["correlation"][wc.correlation]["B"]["variants"],
            feature_set=config["correlation"][wc.correlation]["B"]["feature_set"],
            missing_value=config["correlation"][wc.correlation]["B"]["missing_value"],
        ),
    output:
        "results/correlation/{correlation}/features/join.tsv.gz",
    params:
        how="inner",
        suffixes="_A _B",
        left_on="ID",
        right_on="ID",
    log:
        temp("results/logs/correlation/correlation_featureJoin.{correlation}.log"),
    wrapper:
        getWrapper("file_manipulation/merge")


# correlate two features defined in the config file
rule correlation_correlate_feature:
    input:
        a="results/correlation/{correlation}/features/join.tsv.gz",
    output:
        "results/correlation/{correlation}/features/feature.{featureA}.{featureB}.correlate.tsv.gz",
    params:
        value_a="{featureA}_A",
        value_b="{featureB}_B",
    log:
        temp(
            "results/logs/correlation/correlation_correlate_feature.{correlation}.{featureA}.{featureB}.log"
        ),
    wrapper:
        getWrapper("evaluate/correlation")


# combine correlated features
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
            "Feature_A=%s" % f[0]
            for f in correlation_getCorrelateFeatures(wc.correlation)
        ],
    log:
        temp(
            "results/logs/correlation/correlation_combine_correlate_feature.{correlation}.log"
        ),
    wrapper:
        getWrapper("file_manipulation/concat")


# combine correlated features
rule correlation_plot:
    input:
        lambda wc: expand(
            "results/annotation/{variant_set}/{variant_set}.{feature_set}.{missing_value}.sorted.tsv.gz",
            variant_set=config["correlation"][wc.correlation]["A"]["variants"],
            feature_set=config["correlation"][wc.correlation]["A"]["feature_set"],
            missing_value=config["correlation"][wc.correlation]["B"]["missing_value"],
        ),
        lambda wc: expand(
            "results/annotation/{variant_set}/{variant_set}.{feature_set}.{missing_value}.sorted.tsv.gz",
            variant_set=config["correlation"][wc.correlation]["B"]["variants"],
            feature_set=config["correlation"][wc.correlation]["B"]["feature_set"],
            missing_value=config["correlation"][wc.correlation]["B"]["missing_value"],
        ),
    output:
        "results/correlation/{correlation}/feature.correlate.{method}.png",
    params:
        columns=lambda wc: [
            f[0] for f in correlation_getCorrelateFeatures(wc.correlation)
        ],
        arrange="ID",
        order=False,
        method=lambda wc: wc.method,
        xlab=lambda wc: getCorrPlotXLab(wc.correlation),
        ylab=lambda wc: getCorrPlotYLab(wc.correlation),
    log:
        temp("results/logs/correlation/correlation_plot.{correlation}.{method}.log"),
    wrapper:
        getWrapper("plots/ggcorrplot")
