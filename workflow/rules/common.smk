"""
Common snakemake rule file.

This includes all function.
"""

########################
### Global functions ###
########################


def getWrapper(wrapper):
    """
    Get directory for snakemake wrappers.
    """
    return "file:%s/%s/wrapper.py" % (
        config["global_files"]["wrapper_directory"],
        wrapper,
    )


###########################
### ALL rules functions ###
###########################


def getAllVariantsInput():
    """
    Getter for all variants.
    """
    output = []
    for variant_set in config["variants"]:
        output += getVariantsInput(variant_set, "annotate")
    return output


def getAllAverageFilesForFeatures(genomeBuild):
    """
    Get a list of features for a genome build that are defined in the feature config file as well as in the config file.
    Use them to to return the average feature path.
    """
    featureList = getDefinedFeatures(genomeBuild)
    return expand(
        "results/features/single_vcf/{feature}/{genomeBuild}/single/{feature}.avg.tsv.gz",
        feature=featureList,
        genomeBuild=genomeBuild,
    )


def getAllAnnotations():
    """
    Get all positive and negatuve training annotations
    TODO: Update for other annotations like benchmark sets
    """
    output = []
    for training in config["training"].values():
        output += [
            "results/annotation/%s/%s.%s.sorted.tsv.gz"
            % (training["positives"], training["positives"], training["feature_set"],)
        ]
        output += [
            "results/annotation/%s/%s.%s.sorted.tsv.gz"
            % (training["negatives"], training["negatives"], training["feature_set"],)
        ]
    return output


def getCorrelationPlots():
    output = []
    for correlation, values in config["correlation"].items():
        if "plot" in values["correlate"]:
            output += expand(
                "results/correlation/{correlation}/{input_type}.correlate.{method}.png",
                correlation=correlation,
                input_type=values["correlate"]["plot"],
                method=["spearman", "pearson"],
            )
    return output


#######################
### FEATURE getters ###
#######################


def getDefinedFeatures(genomeBuild):
    """
    Get a list of features for a genome build that are defined in the feature config file as well as in the config file.
    """
    featureList = []
    for feature, build in features.items():
        if genomeBuild in build:
            featureList += [feature]
    return featureList


def getFeaturesOfFeatureSet(feature_set):
    """
    Return the features of a feature set
    """
    return config["feature_sets"][feature_set]["features"]


def getFeatureMissingValue(feature, genome_build, missing_value_config):
    """
    return the default value of a feature (and its genome build)
    """
    return features[feature][genome_build]["missing_value"][missing_value_config]


def getFeaturesOfTraining(training_run):
    feature_set = config["training"][training_run]["feature_set"]
    return getFeaturesOfFeatureSet(feature_set)


def getFeaturesOfScore(score_name):
    training = config["scores"][score_name]["training"]
    feature_set = config["training"][training]["feature_set"]
    return getFeaturesOfFeatureSet(feature_set)


########################
### Training getters ###
########################


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
