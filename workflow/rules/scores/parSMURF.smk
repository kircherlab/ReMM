rule scores_parSMURF_data:
    input:
        "results/scores/{score_name}/input/annotation/{split}.sorted.tsv.gz",
    output:
        temp("results/scores/{score_name}/input/parsmurf/parsmurf.data.{split}.txt"),
    shell:
        """
        zcat {input} | egrep -v "^CHR\sPOSITION\sID" | cut -f 4- > {output}
        """


rule scores_parSMURF_conf:
    input:
        data="results/scores/{score_name}/input/parsmurf/parsmurf.data.{split}.txt",
        models=(
            lambda wc: "results/training/%s/predictions/models/0.out.forest"
            % config["scores"][wc.score_name]["training"]
        ),
        scaffold="resources/scaffold.json",
    output:
        config="results/scores/{score_name}/input/model/parsmurf.config.{split}.json",
    params:
        predictions=(
            "results/scores/{score_name}/predictions/split/predictions_{split}.txt"
        ),
        name="{score_name}",
        models=(
            lambda wc: "results/training/%s/predictions/models"
            % config["scores"][wc.score_name]["training"]
        ),
        seed={"1"},
        mode="test",
    script:
        "../../scripts/generateParsmurfConfig.py"


rule scores_parSMURF_test:
    input:
        data="results/scores/{score_name}/input/parsmurf/parsmurf.data.{split}.txt",
        config="results/scores/{score_name}/input/model/parsmurf.config.{split}.json",
    output:
        "results/scores/{score_name}/predictions/parsmurf/predictions_{split}.txt",
    shell:
        """
        workflow/bin/parSMURF1 --cfg {input.config}
        """


rule scores_parSMURF_combine:
    input:
        predictions=(
            "results/scores/{score_name}/predictions/parsmurf/predictions_{split}.txt"
        ),
        positions="results/scores/{score_name}/input/annotation/{split}.sorted.tsv.gz",
    output:
        "results/scores/{score_name}/predictions/split/predictions_{split}.tsv.gz",
    shell:
        """
        paste \
        <(zcat {input.positions | egrep -v "^CHR\sPOSITION\sID" | cut -f 1,2}) \
        <(cat {input.predictions} | cut -f 1) | \
        bgzip -c > {output}
        """
