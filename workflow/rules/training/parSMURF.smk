
rule training_parSMURF_createParsmurfInput:
    input:
        cb=lambda wc: getTrainingRunFolds(wc.training_run),
        positives=lambda wc: getTrainingPositives(wc.training_run),
        negatives=lambda wc: getTrainingNegatives(wc.training_run),
    output:
        data="results/training/{training_run}/input/parsmurf.data.txt",
        labels="results/training/{training_run}/input/parsmurf.labels.txt",
        folds="results/training/{training_run}/input/parsmurf.folds.txt",
    shell:
        """
        python workflow/scripts/createParsmurfInput.py \
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
        predictions=(
            "results/training/{training_run}/predictions/parsmurf/predictions.txt"
        ),
        name="{training_run}",
        models="results/training/{training_run}/predictions/models",
        seed=lambda wc: config["training"][wc.training_run]["config"]["seed"],
        mode=lambda wc: wc.mode,
    script:
        "../../scripts/generateParsmurfConfig.py"


rule training_parSMURF_generateRepetitiveConfig:
    input:
        data="results/training/{training_run}/input/parsmurf.data.txt",
        labels="results/training/{training_run}/input/parsmurf.labels.txt",
        folds="results/training/{training_run}/input/parsmurf.folds.txt",
        scaffold="resources/scaffold.json",
    output:
        config="results/training/{training_run}/input/repetitive/parsmurf.config.{seed}.{mode}.json",
    params:
        predictions="results/training/{training_run}/predictions/parsmurf/repetitive/predictions.{seed}.txt",
        name=lambda wc: wc.training_run,
        models="results/training/{training_run}/predictions/models",
        seed=lambda wc: wc.seed,
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
        "results/training/{training_run}/predictions/parsmurf/predictions.txt",
    shell:
        """
        workflow/bin/parSMURF1 --cfg {input.config}
        """


rule training_parSMURF_repetitiveCV:
    input:
        data="results/training/{training_run}/input/parsmurf.data.txt",
        labels="results/training/{training_run}/input/parsmurf.labels.txt",
        folds="results/training/{training_run}/input/parsmurf.folds.txt",
        config="results/training/{training_run}/input/repetitive/parsmurf.config.{seed}.cv.json",
    output:
        "results/training/{training_run}/predictions/parsmurf/repetitive/predictions.{seed}.txt",
    shell:
        """
        workflow/bin/parSMURF1 --cfg {input.config}
        """


rule training_parSMURF_combine:
    input:
        predictions=(
            "results/training/{training_run}/predictions/parsmurf/{predictions_alone_or_repetitive}.txt"
        ),
        positives=lambda wc: getTrainingPositives(wc.training_run),
        negatives=lambda wc: getTrainingNegatives(wc.training_run),
    output:
        "results/training/{training_run}/predictions/{predictions_alone_or_repetitive}.tsv.gz",
    wildcard_constraints:
        predictions_alone_or_repetitive="(repetitive/predictions.\d+)|(predictions)",
    shell:
        """
        paste \
        <(zcat {input.positives} {input.negatives} | egrep -v "^CHR\sPOSITION\sID" | cut -f 1,2) \
        <(cat {input.predictions} | cut -f 1 ) | \
        sort -k 1,1 -k2,2n | \
        bgzip -c > {output}
        """


rule training_parSMURF_combine_labels:
    input:
        predictions=(
            "results/training/{training_run}/predictions/parsmurf/{predictions_alone_or_repetitive}.txt"
        ),
        positives=lambda wc: getTrainingPositives(wc.training_run),
        negatives=lambda wc: getTrainingNegatives(wc.training_run),
    output:
        "results/training/{training_run}/predictions/{predictions_alone_or_repetitive}.with_labels.tsv.gz",
    wildcard_constraints:
        predictions_alone_or_repetitive="(repetitive/predictions.\d+)|(predictions)",
    shell:
        """
        (
            echo -e "CHROM\\tPOS\\tLABEL\\tSCORE";
            paste \
            <(
                zcat {input.positives} | egrep -v "^CHR\sPOSITION\sID" | awk -v 'OFS=\\t' '{{print $1,$2,1}}';
                zcat {input.negatives} | egrep -v "^CHR\sPOSITION\sID" | awk -v 'OFS=\\t' '{{print $1,$2,0}}';
            ) \
            <(cat {input.predictions} | cut -f 1 ) | \
            sort -k 1,1 -k2,2n;
        ) | \
        bgzip -c > {output}
        """


rule training_parSMURF_train:
    input:
        data="results/training/{training_run}/input/parsmurf.data.txt",
        labels="results/training/{training_run}/input/parsmurf.labels.txt",
        folds="results/training/{training_run}/input/parsmurf.folds.txt",
        config="results/training/{training_run}/input/parsmurf.config.train.json",
    output:
        expand(
            "results/training/{{training_run}}/predictions/models/{model_number}.out.forest",
            model_number=list(range(0, 100)),
        ),
        expand(
            "results/training/{{training_run}}/predictions/models/{model_number}.out.importance",
            model_number=list(range(0, 100)),
        ),
    params:
        models="results/training/{training_run}/predictions/models/",
    shell:
        """
        mkdir -p {params.models};
        workflow/bin/parSMURF1 --cfg {input.config}
        """
