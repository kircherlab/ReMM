################################
#### training sub workflow  ####
################################


"""
Results will be saved in results/training/<training_run>
"""

def getFeaturesOfTraining(training_run):
    feature_set = config["training"][training_run]["feature_set"]
    return getFeaturesOfFeatureSet(feature_set)

def getTrainingRunGenomeBuild(training_run):
    """
    Get the genome buld of training (not specified by training but variant set and feature set config)
    """
    conf = config["training"][training_run]

    genomeBuild_positives = getVariantSetGenomeBuild(conf["positives"])
    genomeBuild_negatives = getVariantSetGenomeBuild(conf["negatives"])
    genomeBuild_feature_set = config["feature_sets"][conf["feature_set"]][
        "genome_build"
    ]

    if (
        genomeBuild_positives != genomeBuild_negatives
        or genomeBuild_negatives != genomeBuild_feature_set
    ):
        raise Exception(
            "genome builds are different for training %s. positives: %s negatives: %s feature set: %s"
            % (
                training_run,
                genomeBuild_positives,
                genomeBuild_negatives,
                genomeBuild_feature_set,
            )
        )
    else:
        return genomeBuild_positives


# check if positie/negatve and feature set have the same genome build
for training_run in config["training"].keys():
    getTrainingRunGenomeBuild(training_run)


def getTrainingData(training_run, label):
    return expand(
        "results/annotation/{variant_set_positive}/{variant_set_positive}.{feature_set}.{missing_value}.sorted.tsv.gz",
        variant_set_positive=config["training"][training_run][label],
        feature_set=config["training"][training_run]["feature_set"],
        missing_value=config["training"][training_run]["missing_value"],
    )


def getTrainingPositives(training_run):
    return getTrainingData(training_run, "positives")


def getTrainingNegatives(training_run):
    return getTrainingData(training_run, "negatives")


rule getFolds:
    input:
        m="resources/ReMM.20171122.partition.mapping.tsv.gz",
        c="resources/{genomeBuild}/cytoBand.txt.gz",
    output:
        "resources/{genomeBuild}/folds.txt.gz",
    script:
        "../scripts/createFolds.py"


def getTrainingRunFolds(training_run):
    return "resources/%s/folds.txt.gz" % getTrainingRunGenomeBuild(training_run)


include: "training/parSMURF.smk"
