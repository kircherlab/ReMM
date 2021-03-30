# TODO missing script createHyperSmurfInput.py


rule combineHyperSmurfInputData:
    input:
        f="results/variants/{genomeBuild}/SNVs.{genomeBuild}.{type}.annotated.tsv.gz",
        cb="results/folds/{genomeBuild}/folds.{genomeBuild}.txt.gz",
    output:
        "results/features/annotated/hyperSMURF/SNVs.{genomeBuild}.hyperSMURF.{type}.txt.gz",
    script:
        "../../scripts/createHyperSmurfInput.py"


rule hyperSMURF:
    input:
        hyperSMURF="workflow/bin/remm-cli-0.0.4-SNAPSHOT-jar-with-dependencies.jar",
        config_positives="config/features.positives.properties",
        config_negatives="config/features.negatives.properties",
        config_hypersmurf="config/classifier.hyperSMURF.properties",
    output:
        model=(
            "results/predictions/{genomeBuild}/trained.hyperSMURF.{genomeBuild}.model"
        ),
        arff="results/predictions/{genomeBuild}/trained.hyperSMURF.{genomeBuild}.{genomeBuild}.arff.gz",
        probabilities=(
            "results/predictions/{genomeBuild}/trained.{genomeBuild}.hyperSMURF.tsv.gz"
        ),
        log="results/predictions/{genomeBuild}/trained.hyperSMURF.{genomeBuild}.log",
    params:
        threads=30,
        mem="100g",
        seed=1234,
        folds=10,
    shell:
        """
        java -Xmx{params.mem} -jar {input.hyperSMURF} build --cores {params.threads} --seed {params.seed} \
        --classifier {input.config_hypersmurf} \
        --evaluate-folds {params.folds} \
        --data {input.config_positives} --data {input.config_negatives} \
        --evaluate-output  {output.arff} --output {output.probabilities} --serialize-classifier {output.model} \
        > {output.log}
        """
