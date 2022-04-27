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


# check if positive/negative and feature set have the same genome build
for training_run in config["training"].keys():
    getTrainingRunGenomeBuild(training_run)


def getTrainingRunFolds(training_run):
    return "resources/%s/folds.txt.gz" % getTrainingRunGenomeBuild(training_run)
