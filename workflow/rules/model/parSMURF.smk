rule generateParsmurfConfig:
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
        """


# output/predictions/hg38/SNVs.hg38.predictions.txt
