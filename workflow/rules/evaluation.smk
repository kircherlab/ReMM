rule evaluation_auc:
    input:
        "results/training/{training_run}/predictions/predictions_with_labels.tsv.gz",
    output:
        "results/evaluation/{training_run}/metrics/auc.tsv.gz",
    params:
        label_column="LABEL",
        prediction_column="SCORE",
        positive_label=1,
    wrapper:
        getWrapper("evaluate/auc")


rule evaluation_plot_metrics:
    input:
        "results/training/{training_run}/predictions/predictions_with_labels.tsv.gz",
    output:
        "results/evaluation/{training_run}/metrics/curves_prc_roc.png",
    params:
        label_column="LABEL",
        score_column="SCORE",
        positive_label=1,
    wrapper:
        getWrapper("plots/metric_curves")
