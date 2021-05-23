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
