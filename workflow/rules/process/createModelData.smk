rule process_createModelData_combineInputData:
    conda:
        "../../envs/default.yml"
    input:
        p="results/variants/{genomeBuild}/SNVs.{genomeBuild}.positive.annotated.tsv.gz",
        n="results/variants/{genomeBuild}/SNVs.{genomeBuild}.negative.annotated.tsv.gz",
    output:
        o="results/features/annotated/{genomeBuild}/SNVs.{genomeBuild}.combined.txt.gz",
        p=temp("results/features/annotated/{genomeBuild}/p.txt"),
        n=temp("results/features/annotated/{genomeBuild}/n.txt"),
    log:
        temp("results/logs/process/createModelData/combineInputData.{genomeBuild}.log"),
    shell:
        """
        zcat {input.p} | tail -n +2 | awk '{{$0="1\t"$0}}1'  > {output.p};
        zcat {input.n} | tail -n +2 |awk '{{$0="0\t"$0}}1'  > {output.n};
        cat {output.p} {output.n} | bgzip -c > {output.o};        
        """


rule process_createModelData_createParsmurfInput:
    conda:
        "../../envs/default.yml"
    input:
        cb="results/folds/{genomeBuild}/folds.{genomeBuild}.txt.gz",
        f="results/features/annotated/{genomeBuild}/SNVs.{genomeBuild}.combined.txt.gz",
    output:
        d="results/features/annotated/{genomeBuild}/SNVs.{genomeBuild}.data.txt",
        l="results/features/annotated/{genomeBuild}/SNVs.{genomeBuild}.labels.txt",
        f="results/features/annotated/{genomeBuild}/SNVs.{genomeBuild}.folds.txt",
    log:
        temp(
            "results/logs/process/createModelData/createParsmurfInput.{genomeBuild}.log"
        ),
    script:
        "../../scripts/createParsmurfInput.py"
