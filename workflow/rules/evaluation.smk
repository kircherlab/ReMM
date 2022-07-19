###################################
#### evaluation sub workflow  ####
###################################


"""
Results will be saved in results/evaluation/<training_run>

Evaluation of training runs like AUPRC, AUROC or curves. 
Also mean AURPRC and AUROC for repetetive (100 times) CV training.
"""


include: "evaluation_defs.smk"


### metrics ###


# create AUC metrics for a run
rule evaluation_auc:
    input:
        "results/training/{training_run}/predictions/predictions.with_labels.tsv.gz",
    output:
        "results/evaluation/{training_run}/metrics/auc.tsv.gz",
    params:
        label_column="LABEL",
        prediction_column="SCORE",
        positive_label=1,
    log:
        temp("results/logs/evaluation/auc.{training_run}.log"),
    wrapper:
        getWrapper("evaluate/auc")


# create AUC metrics for a repetitive run
rule evaluation_aucRepetitive:
    input:
        "results/training/{training_run}/predictions/repetitive/predictions.{seed}.with_labels.tsv.gz",
    output:
        temp(
            "results/evaluation/{training_run}/metrics/repetitive/auc_tmp.{seed}.tsv.gz"
        ),
    params:
        label_column="LABEL",
        prediction_column="SCORE",
        positive_label=1,
    log:
        temp("results/logs/evaluation/aucRepetitive.{training_run}.{seed}.log"),
    wrapper:
        getWrapper("evaluate/auc")


# rename header for repetetive metrics (before merging)
rule evaluation_aucRepetitiveRename:
    input:
        "results/evaluation/{training_run}/metrics/repetitive/auc_tmp.{seed}.tsv.gz",
    output:
        temp("results/evaluation/{training_run}/metrics/repetitive/auc.{seed}.tsv.gz"),
    params:
        columns=lambda wc: {"value": "seed_%d" % int(wc.seed)},
    log:
        temp("results/logs/evaluation/aucRepetitiveRename.{training_run}.{seed}.log"),
    wrapper:
        getWrapper("file_manipulation/rename")


# combine repetitive metrics
rule evaluation_aucRepetitiveCombine:
    input:
        lambda wc: expand(
            "results/evaluation/{{training_run}}/metrics/repetitive/auc.{seed}.tsv.gz",
            seed=getSeedsForTraining(wc.training_run, 100),
        ),
    output:
        "results/evaluation/{training_run}/metrics/repetitive/auc_all.tsv.gz",
    params:
        index="metric",
    log:
        temp("results/logs/evaluation/aucRepetitiveCombine.{training_run}.log"),
    wrapper:
        getWrapper("file_manipulation/concat")


# create mean max min for all repetitive runs
rule evaluation_aucRepetitiveMean:
    input:
        "results/evaluation/{training_run}/metrics/repetitive/auc_all.tsv.gz",
    output:
        "results/evaluation/{training_run}/metrics/repetitive/auc_all.mean_std_min_max.tsv.gz",
    params:
        columns=lambda wc: [
            "seed_%d" % seed for seed in getSeedsForTraining(wc.training_run, 100)
        ],
        new_columns=["mean", "std", "min", "max"],
        operations=["mean", "std", "min", "max"],
    log:
        temp("results/logs/evaluation/aucRepetitiveMean.{training_run}.log"),
    wrapper:
        getWrapper("file_manipulation/summarize_columns")


# create metrcis table per theshold for single run


rule evaluation_metrics_per_threshold:
    input:
        "results/training/{training_run}/predictions/predictions.with_labels.tsv.gz",
    output:
        "results/evaluation/{training_run}/metrics/metrics_per_threshold.tsv.gz",
    params:
        label_column="LABEL",
        prediction_column="SCORE",
        positive_label=1,
        decimals=3,
    log:
        temp("results/logs/evaluation/metrics_per_threshold.{training_run}.log"),
    wrapper:
        getWrapper("evaluate/metrics_per_threshold")


### plots ###


# plot prc and roc curve for training
rule evaluation_plot_metrics:
    input:
        "results/training/{training_run}/predictions/predictions.with_labels.tsv.gz",
    output:
        "results/evaluation/{training_run}/metrics/curves_prc_roc.png",
    params:
        label_column="LABEL",
        score_column="SCORE",
        positive_label=1,
    log:
        temp("results/logs/evaluation/plot_metrics.{training_run}.log"),
    wrapper:
        getWrapper("plots/metric_curves")


rule evaluation_plot_pre_re_f1_f2:
    input:
        "results/evaluation/{training_run}/metrics/metrics_per_threshold.tsv.gz",
    output:
        "results/evaluation/{training_run}/metrics/pre_re_f1_f2.png",
    params:
        xname="'ReMM score'",
    log:
        temp("results/logs/evaluation/plot_pre_re_f1_f2.{training_run}.log"),
    wrapper:
        getWrapper("plots/pre_re_f1_f2")


rule evaluation_plot_combined:
    input:
        lambda wc: expand(
            "results/evaluation/{model}/metrics/metrics_per_threshold.tsv.gz",
            model=config["evaluation_combined"][wc.evaluation]["models"],
        ),
    output:
        "results/evaluation_combined/{evaluation}/{type}.png",
    params:
        xname=lambda wc: "'Precision'" if wc.type == "PR" else "'False positive rate'",
        yname=lambda wc: "'Recall'" if wc.type == "PR" else "'True positive rate'",
        names=lambda wc: config["evaluation_combined"][wc.evaluation]["names"][wc.type],
        type=lambda wc: wc.type,
    log:
        temp("results/logs/evaluation/evaluation_plot_combined.{evaluation}.{type}.log"),
    wrapper:
        getWrapper("plots/PR_ROC_curves_metric")
