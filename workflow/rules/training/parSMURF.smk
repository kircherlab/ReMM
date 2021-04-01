rule training_parSMURF_combineInputData:
    input:
        positives="results/annotation/{variant_set_positive}/{variant_set_positive}.{feature_set}.tsv.gz",
        negatives="results/annotation/{variant_set_negative}/{variant_set_negative}.{feature_set}.tsv.gz",
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
        cb="results/folds/{genomeBuild}/folds.{genomeBuild}.txt.gz",
        f="results/training/{training_run}/input/parsmurf.combined.txt.gz",
    output:
        d="results/training/{training_run}/input/parsmurf.data.txt",
        l="results/training/{training_run}/input/parsmurf.labels.txt",
        f="results/training/{training_run}/input/parsmurf.folds.txt",
    script:
        "../../scripts/createParsmurfInput.py"


rule training_parSMURF_generateConfig:
    input:
        data="results/features/annotated/{genomeBuild}/SNVs.{genomeBuild}.data.txt",
        labels=(
            "results/features/annotated/{genomeBuild}/SNVs.{genomeBuild}.labels.txt"
        ),
        folds="results/features/annotated/{genomeBuild}/SNVs.{genomeBuild}.folds.txt",
        scaffold="resources/scaffold.json",
    output:
        config="config/SNVs.{genomeBuild}.{mode}.json",
    params:
        predictions="results/predictions/{genomeBuild}/SNVs.{genomeBuild}.{mode}.predictions.txt",
        name="SNVs.{genomeBuild}",
        #seed = lambda wc: {output.split('_')[-1]} if output.split('_')[-1].isdigit() else {'1'},
        seed={"1"},
        mode="{mode}",
    script:
        "../../scripts/generateParsmurfConfig.py"


rule runParSMURF:
    input:
        "results/features/annotated/{genomeBuild}/SNVs.{genomeBuild}.data.txt",
        "results/features/annotated/{genomeBuild}/SNVs.{genomeBuild}.labels.txt",
        "results/features/annotated/{genomeBuild}/SNVs.{genomeBuild}.folds.txt",
        "config/SNVs.{genomeBuild}.{mode}.json",
    output:
        "results/predictions/{genomeBuild}/SNVs.{genomeBuild}.{mode}.predictions.txt",
    shell:
        """
        workflow/bin/parSMURF1 --cfg {input[3]}
