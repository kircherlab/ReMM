###############################################
# Create a final whole genome wide score
###############################################

scores_NumberOfSplits = 3000


rule scores_getNotNRegions:
    conda:
        "../envs/default.yml"
    input:
        lambda wc: config["global_files"]["genome_builds"][wc.genomeBuild]["reference"],
    output:
        "resources/{genomeBuild}/notNRegions.bed",
    log:
        temp("results/logs/scores/getNotNRegions.{genomeBuild}.log"),
    shell:
        """
        python workflow/scripts/notNRegionsOfGenome.py --input {input} --output {output} &> {log}
        """


rule scores_split_notN_regions:
    conda:
        "../envs/default.yml"
    input:
        lambda wc: expand(
            "resources/{genomeBuild}/notNRegions.bed",
            genomeBuild=getTrainingRunGenomeBuild(
                config["scores"][wc.score_name]["training"]
            ),
        ),
    output:
        regions=temp(
            expand(
                "results/scores/{{score_name}}/input/bed/split.{split}.bed",
                split=["%.5d" % i for i in range(1, scores_NumberOfSplits + 1)],
            )
        ),
    params:
        prefix=lambda wc: "results/scores/%s/input/bed/split" % wc.score_name,
        dir_name=lambda wc: "results/scores/%s/input/bed" % wc.score_name,
        files=scores_NumberOfSplits,
    log:
        temp("results/logs/scores/split_notN_regions.{score_name}.log"),
    shell:
        """
        mkdir -p {params.dir_name};
        bedtools makewindows -b {input} -w 500000 | \
        bedtools split -i - -n {params.files} -p {params.prefix} &> {log}
        """


rule scores_extract_features:
    conda:
        "../envs/default.yml"
    input:
        feature_set=lambda wc: expand(
            "results/features/feature_sets/{feature_set}.vcf.gz",
            feature_set=config["training"][
                config["scores"][wc.score_name]["training"]
            ]["feature_set"],
        ),
        feature_set_idx=lambda wc: expand(
            "results/features/feature_sets/{feature_set}.vcf.gz.tbi",
            feature_set=config["training"][
                config["scores"][wc.score_name]["training"]
            ]["feature_set"],
        ),
        regions="results/scores/{score_name}/input/bed/split.{split}.bed",
    output:
        temp("results/scores/{score_name}/input/annotation/{split}.unsorted.tsv.gz"),
    log:
        temp("results/logs/scores/extract_features.{score_name}.{split}.log"),
    shell:
        """
        python workflow/scripts/getAnnotationsByInterval.py --input {input.feature_set} --output {output} \
        --regions {input.regions} &> {log}
        """


rule scores_sort_features:
    conda:
        "../envs/default.yml"
    input:
        "results/scores/{score_name}/input/annotation/{split}.unsorted.tsv.gz",
    output:
        temp("results/scores/{score_name}/input/annotation/{split}.sorted.tsv.gz"),
    params:
        features=lambda wc: " ".join(
            [
                "--feature %s %f"
                % (
                    feature,
                    getFeatureMissingValue(
                        feature,
                        getTrainingRunGenomeBuild(
                            config["scores"][wc.score_name]["training"]
                        ),
                        config["training"][
                            config["scores"][wc.score_name]["training"]
                        ]["missing_value"],
                    ),
                )
                for feature in getFeaturesOfScore(wc.score_name)
            ]
        ),
    log:
        temp("results/logs/scores/sort_features.{score_name}.{split}.log"),
    shell:
        """
        python workflow/scripts/sortAnnotationFile.py --input {input} --output {output} {params.features} &> {log}
        """


include: "scores/parSMURF.smk"


rule scores_combineScores:
    conda:
        "../envs/default.yml"
    input:
        scores=expand(
            "results/scores/{{score_name}}/predictions/split/predictions_{split}.tsv.gz",
            split=["%.5d" % i for i in range(1, scores_NumberOfSplits + 1)],
        ),
    output:
        prediction="results/scores/{score_name}/release/{score_name}.biased.tsv.gz",
    log:
        temp("results/logs/scores/combineScores.{score_name}.log"),
    shell:
        """
        export LC_ALL=C;
        zcat {input.scores} | sort -k1,1 -k2,2n | bgzip -c > {output.prediction} 2> {log};
        """


rule scores_replaceScores:
    conda:
        "../envs/default.yml"
    input:
        biased_score="results/scores/{score_name}/release/{score_name}.biased.tsv.gz",
        predictions=lambda wc: expand(
            "results/training/{training_run}/predictions/predictions.tsv.gz",
            training_run=config["scores"][wc.score_name]["training"],
        ),
    output:
        score="results/scores/{score_name}/release/{score_name}.unbiased.tsv.gz",
    params:
        comments=lambda wc: " ".join(
            ["--comment '%s'" % i for i in config["scores"][wc.score_name]["comments"]]
        ),
    log:
        temp("results/logs/scores/replaceScores.{score_name}.log"),
    shell:
        """
        zcat {input.biased_score} | \
        python workflow/scripts/replaceScores.py \
        --replace-score-file {input.predictions} \
        {params.comments} \
        | uniq | bgzip -c > {output.score} 2> {log};
        """


rule scores_indexScore:
    conda:
        "../envs/default.yml"
    input:
        "results/scores/{score_name}/release/{score_name}.unbiased.tsv.gz",
    output:
        "results/scores/{score_name}/release/{score_name}.unbiased.tsv.gz.tbi",
    log:
        temp("results/logs/scores/indexScore.{score_name}.log"),
    shell:
        """
        tabix -s 1 -b 2 -e 2 -c "#" {input} &> {log}
        """
