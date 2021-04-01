rule training_parSMURF_combineInputData:
    input:
        positives=lambda wc: expand(
            "results/annotation/{variant_set_positive}/{variant_set_positive}.{feature_set}.tsv.gz",
            variant_set_positive=config["training"][wc.training_run]["positives"],
            feature_set=config["training"][wc.training_run]["feature_set"],
        ),
        negatives=lambda wc: expand(
            "results/annotation/{variant_set_negative}/{variant_set_negative}.{feature_set}.tsv.gz",
            variant_set_negative=config["training"][wc.training_run]["negatives"],
            feature_set=config["training"][wc.training_run]["feature_set"],
        ),
    output:
        "results/training/{training_run}/input/parsmurf.combined.txt.gz",
    shell:
        """
        set +o pipefail;
        (
            zcat {input.positives} | tail -n +2 | awk '{{$0="1\t"$0}}1';
            zcat {input.negatives} | tail -n +2 |awk '{{$0="0\t"$0}}1';
        ) | bgzip -c > {output};
        """


rule training_parSMURF_createParsmurfInput:
    input:
        cb=lambda wc: getTrainingRunFolds(wc.training_run),
        f="results/training/{training_run}/input/parsmurf.combined.txt.gz",
    output:
        d="results/training/{training_run}/input/parsmurf.data.txt",
        l="results/training/{training_run}/input/parsmurf.labels.txt",
        f="results/training/{training_run}/input/parsmurf.folds.txt",
    script:
        "../../scripts/createParsmurfInput.py"


rule training_parSMURF_generateConfig:
    input:
        data="results/training/{training_run}/input/parsmurf.data.txt",
        labels="results/training/{training_run}/input/parsmurf.labels.txt",
        folds="results/training/{training_run}/input/parsmurf.folds.txt",
        scaffold="resources/scaffold.json",
    output:
        config="results/training/{training_run}/input/parsmurf.config.json",
    params:
        predictions="results/training/{training_run}/predictions/predictions.txt",
        name="{training_run}",
        seed={"1"},
        mode=lambda wc: config["training"][wc.training_run]["config"]["mode"],
    script:
        "../../scripts/generateParsmurfConfig.py"


rule training_parSMURF_run:
    input:
        data="results/training/{training_run}/input/parsmurf.data.txt",
        labels="results/training/{training_run}/input/parsmurf.labels.txt",
        folds="results/training/{training_run}/input/parsmurf.folds.txt",
        config="results/training/{training_run}/input/parsmurf.config.json",
    output:
        "results/training/{training_run}/predictions/predictions.txt",
    shell:
        """
        workflow/bin/parSMURF1 --cfg {input.config}
        """
