rule scores_parSMURF_data:
    conda:
        "../../envs/default.yml"
    input:
        "results/scores/{score_name}/input/annotation/{split}.sorted.tsv.gz",
    output:
        temp("results/scores/{score_name}/input/parsmurf/parsmurf.data.{split}.txt"),
    log:
        temp("results/logs/scores/parSMURF_data.{score_name}.{split}.log"),
    shell:
        """
        zcat {input} | egrep -v "^CHR\sPOSITION\sID" | cut -f 4- > {output} 2> {log}
        """


rule scores_parSMURF_labels:
    conda:
        "../../envs/default.yml"
    input:
        "results/scores/{score_name}/input/annotation/{split}.sorted.tsv.gz",
    output:
        temp("results/scores/{score_name}/input/parsmurf/parsmurf.labels.{split}.txt"),
    log:
        temp("results/logs/scores/parSMURF_labels.{score_name}.{split}.log"),
    shell:
        """
        zcat {input} | egrep -v "^CHR\sPOSITION\sID" | awk '{{print 1}}' > {output} 2> {log}
        """


rule scores_parSMURF_conf:
    conda:
        "../../envs/default.yml"
    input:
        data="results/scores/{score_name}/input/parsmurf/parsmurf.data.{split}.txt",
        models=(
            lambda wc: "results/training/%s/predictions/models/0.out.forest"
            % config["scores"][wc.score_name]["training"]
        ),
        labels=(
            "results/scores/{score_name}/input/parsmurf/parsmurf.labels.{split}.txt"
        ),
        scaffold="resources/scaffold.json",
    output:
        config=temp(
            "results/scores/{score_name}/input/model/parsmurf.config.{split}.json"
        ),
    params:
        predictions=(
            "results/scores/{score_name}/predictions/parsmurf/predictions_{split}.txt"
        ),
        name="{score_name}",
        models=(
            lambda wc: "results/training/%s/predictions/models"
            % config["scores"][wc.score_name]["training"]
        ),
        seed="1",
        mode="predict",
        ensThrd="30",
    log:
        temp("results/logs/scores/parSMURF_conf.{score_name}.{split}.log"),
    script:
        "../../scripts/generateParsmurfConfig.py"


rule scores_parSMURF_test:
    conda:
        "../../envs/default.yml"
    input:
        data="results/scores/{score_name}/input/parsmurf/parsmurf.data.{split}.txt",
        config="results/scores/{score_name}/input/model/parsmurf.config.{split}.json",
        labels=(
            "results/scores/{score_name}/input/parsmurf/parsmurf.labels.{split}.txt"
        ),
        models=(
            lambda wc: "results/training/%s/predictions/models/0.out.forest"
            % config["scores"][wc.score_name]["training"]
        ),
    output:
        temp("results/scores/{score_name}/predictions/parsmurf/predictions_{split}.txt"),
    log:
        temp("results/logs/scores/parSMURF_test.{score_name}.{split}.log"),
    shell:
        """
        workflow/bin/parSMURF1 --cfg {input.config} > {log}
        """


rule scores_parSMURF_combine:
    conda:
        "../../envs/default.yml"
    input:
        predictions=(
            "results/scores/{score_name}/predictions/parsmurf/predictions_{split}.txt"
        ),
        positions="results/scores/{score_name}/input/annotation/{split}.sorted.tsv.gz",
    output:
        temp("results/scores/{score_name}/predictions/split/predictions_{split}.tsv.gz"),
    log:
        temp("results/logs/scores/parSMURF_combine.{score_name}.{split}.log"),
    shell:
        """
        paste \
        <(zcat {input.positions} | egrep -v "^CHR\sPOSITION\sID" | cut -f 1,2) \
        <(cat {input.predictions} | cut -f 1 ) | \
        bgzip -c > {output} 2> {log}
        """
