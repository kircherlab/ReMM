
rule training_parSMURF_createParsmurfInput:
    input:
        cb=lambda wc: getTrainingRunFolds(wc.training_run),
        positives=lambda wc: expand(
            "results/annotation/{variant_set_positive}/{variant_set_positive}.{feature_set}.sorted.tsv.gz",
            variant_set_positive=config["training"][wc.training_run]["positives"],
            feature_set=config["training"][wc.training_run]["feature_set"],
        ),
        negatives=lambda wc: expand(
            "results/annotation/{variant_set_negative}/{variant_set_negative}.{feature_set}.sorted.tsv.gz",
            variant_set_negative=config["training"][wc.training_run]["negatives"],
            feature_set=config["training"][wc.training_run]["feature_set"],
        ),
    output:
        data="results/training/{training_run}/input/parsmurf.data.txt",
        labels="results/training/{training_run}/input/parsmurf.labels.txt",
        folds="results/training/{training_run}/input/parsmurf.folds.txt",
    shell:
        """
        workflow/scripts/createParsmurfInput.py \
        --folds {input.cb} --positives {input.positives} --negatives {input.negatives} \
        --output-folds {output.folds} --output-data {output.data} --output-labels {output.labels}
        """


rule training_parSMURF_generateConfig:
    input:
        data="results/training/{training_run}/input/parsmurf.data.txt",
        labels="results/training/{training_run}/input/parsmurf.labels.txt",
        folds="results/training/{training_run}/input/parsmurf.folds.txt",
        scaffold="resources/scaffold.json",
    output:
        config="results/training/{training_run}/input/parsmurf.config.{mode}.json",
    params:
        predictions="results/training/{training_run}/predictions/predictions.txt",
        name="{training_run}",
        models="results/training/{training_run}/predictions",
        seed={"1"},
        mode=lambda wc: wc.mode,
    script:
        "../../scripts/generateParsmurfConfig.py"


rule training_parSMURF_cv:
    input:
        data="results/training/{training_run}/input/parsmurf.data.txt",
        labels="results/training/{training_run}/input/parsmurf.labels.txt",
        folds="results/training/{training_run}/input/parsmurf.folds.txt",
        config="results/training/{training_run}/input/parsmurf.config.cv.json",
    output:
        "results/training/{training_run}/predictions/predictions.txt",
    shell:
        """
        workflow/bin/parSMURF1 --cfg {input.config}
        """


rule training_parSMURF_train:
    input:
        data="results/training/{training_run}/input/parsmurf.data.txt",
        labels="results/training/{training_run}/input/parsmurf.labels.txt",
        folds="results/training/{training_run}/input/parsmurf.folds.txt",
        config="results/training/{training_run}/input/parsmurf.config.train.json",
    output:
        "results/training/{training_run}/predictions/models/0.out.forest",
    params:
        models="results/training/{training_run}/predictions/models/",
    shell:
        """
        mkdir -p {params.models};
        workflow/bin/parSMURF1 --cfg {input.config}
        """
