###############################################
# Create a final whole genome wide score
###############################################


def getFeaturesOfScore(score_name):
    training = config["scores"][score_name]["training"]
    feature_set = config["training"][training]["feature_set"]
    return getFeaturesOfFeatuteSet(feature_set)


rule scores_getNotNRegions:
    input:
        lambda wc: config["global_files"]["genome_builds"][wc.genomeBuild]["reference"],
    output:
        "resources/{genomeBuild}/notNRegions.bed",
    shell:
        """
        python workflow/scripts/notNRegionsOfGenome.py --input {input} --output {output} 
        """


checkpoint scores_split_notN_regions:
    input:
        lambda wc: expand(
            "resources/{genomeBuild}/notNRegions.bed",
            genomeBuild=getTrainingRunGenomeBuild(
                config["scores"][wc.score_name]["training"]
            ),
        ),
    output:
        regions=directory("results/scores/{score_name}/input/bed"),
    params:
        prefix="results/scores/{score_name}/input/bed/split",
        dir_name="results/scores/{score_name}/input/bed",
    shell:
        """
        mkdir {params.dir_name};
        bedtools makewindows -b {input} -w 500000 | \
        bedtools split -i - -n 3000 -p {params.prefix}
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


# TODO add default value for feature
rule scores_sort_features:
    input:
        "results/scores/{score_name}/input/annotation/{split}.unsorted.tsv.gz",
    output:
        "results/scores/{score_name}/input/annotation/{split}.sorted.tsv.gz",
    params:
        features=lambda wc: " ".join(
            [
                "--feature %s %f"
                % (
                    feature,
                    getFeatureMissingValue(
                        feature,
                        config["scores"][wc.score_name]["genome_build"],
                        config["training"][config["scores"][wc.score]["training"]][
                            "missing_value"
                        ],
                    ),
                )
                for feature in getFeaturesOfScore(wc.feature_set)
            ]
        ),
    shell:
        """
        python workflow/scripts/sortAnnotationFile.py --input {input} --output {output} {params.features}
        """


include: "scores/parSMURF.smk"


def aggregate_PredictedScoresPerInterval(wc):
    checkpoint_output = checkpoints.scores_split_notN_regions.get(**wc).output[0]
    return expand(
        "results/scores/{score_name}/predictions/split/predictions_{split}.tsv.gz",
        score_name=wc.score_name,
        split=glob_wildcards(os.path.join(checkpoint_output, "split.{i}.bed")).i,
    )


rule scores_combineScores:
    input:
        aggregate_PredictedScoresPerInterval,
    output:
        prediction="results/scores/{score_name}/release/{score_name}.biased.tsv.gz",
    shell:
        """
        export LC_ALL=C;
        zcat {input} | sort -k1,1 -k2,2n | bgzip -c > {output.prediction};
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
