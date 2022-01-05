

rule featureImportance_header:
    input:
        "results/training/{training}/predictions/models/{model_number}.out.importance",
    output:
        "results/training/{training}/feature_importance/single/{model_number}.importance.tsv.gz",
    params:
        model_number=lambda wc: wc.model_number,
    shell:
        """
        (
            echo -e "Feature\\tImportance_{params.model_number}";
            cat {input} | sed 's/\:\s/\\t/g';
        ) | gzip -c > {output}
        """


rule featureImportance_concat:
    input:
        expand(
            "results/training/{{training}}/feature_importance/single/{model_number}.importance.tsv.gz",
            model_number=list(range(0, 100)),
        ),
    output:
        "results/training/{training}/feature_importance/feature_importance.tsv.gz",
    params:
        index="Feature",
    wrapper:
        getWrapper("file_manipulation/concat")


rule featureImportance_replace:
    input:
        "results/training/{training}/feature_importance/feature_importance.tsv.gz",
    output:
        "results/training/{training}/feature_importance/feature_importance.names.tsv.gz",
    params:
        columns=lambda wc: [
            "Feature" for i in range(len(getFeaturesOfTraining(wc.training)))
        ],
        pat=lambda wc: [
            "%d" % i for i in range(len(getFeaturesOfTraining(wc.training)))
        ],
        replace=lambda wc: getFeaturesOfTraining(wc.training),
    wrapper:
        getWrapper("file_manipulation/replace")


rule featureImportance_aucRepetitiveMean:
    input:
        "results/training/{training}/feature_importance/feature_importance.names.tsv.gz",
    output:
        "results/training/{training}/feature_importance/feature_importance.names.mean_max_min.tsv.gz",
    params:
        columns=lambda wc: ["Importance_%d" % num for num in range(100)],
        new_columns=["mean", "max", "min"],
        operations=["mean", "max", "min"],
    wrapper:
        getWrapper("file_manipulation/summarize_columns")
