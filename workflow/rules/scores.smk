###############################################
# Create a final whole genome wide score
###############################################

scores_NumberOfSplits = 3000


def getFeaturesOfScore(score_name):
    training = config["scores"][score_name]["training"]
    feature_set = config["training"][training]["feature_set"]
    return getFeaturesOfFeatureSet(feature_set)


rule scores_getNotNRegions:
    input:
        lambda wc: config["global_files"]["genome_builds"][wc.genomeBuild]["reference"],
    output:
        "resources/{genomeBuild}/notNRegions.bed",
    shell:
        """
        python workflow/scripts/notNRegionsOfGenome.py --input {input} --output {output} 
        """


rule scores_split_notN_regions:
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
        prefix="results/scores/{score_name}/input/bed/split",
        dir_name="results/scores/{score_name}/input/bed",
        files=scores_NumberOfSplits,
    shell:
        """
        mkdir -p {params.dir_name};
        bedtools makewindows -b {input} -w 500000 | \
        bedtools split -i - -n {params.files} -p {params.prefix}
        """


rule scores_extract_features:
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
    shell:
        """
        python workflow/scripts/getAnnotationsByInterval.py --input {input.feature_set} --output {output} \
        --regions {input.regions}
        """


rule scores_sort_features:
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
    shell:
        """
        python workflow/scripts/sortAnnotationFile.py --input {input} --output {output} {params.features}
        """


include: "scores/parSMURF.smk"


rule scores_combineScores:
    input:
        scores=expand(
            "results/scores/{{score_name}}/predictions/split/predictions_{split}.tsv.gz",
            split=["%.5d" % i for i in range(1, scores_NumberOfSplits + 1)],
        ),
    output:
        prediction="results/scores/{score_name}/release/{score_name}.biased.tsv.gz",
    shell:
        """
        export LC_ALL=C;
        zcat {input.scores} | sort -k1,1 -k2,2n | bgzip -c > {output.prediction};
        """


rule scores_replaceScores:
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
    shell:
        """
        zcat {input.biased_score} | \
        python workflow/scripts/replaceScores.py \
        --replace-score-file {input.predictions} \
        {params.comments} \
        | bgzip -c > {output.score};
        """


rule scores_indexScore:
    input:
        "results/scores/{score_name}/release/{score_name}.unbiased.tsv.gz",
    output:
        "results/scores/{score_name}/release/{score_name}.unbiased.tsv.gz.tbi",
    shell:
        """
        tabix -s 1 -b 2 -e 2 -c "#" {input}
        """
