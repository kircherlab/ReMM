def correlation_getCorrelateFeatures(correlation):
    """
    Function to collect features that should be correlated from both feature sets-
    Raises exception when feature not present in feature set.
    """
    features = [
        i.split("=")
        for i in config["correlation"][correlation]["correlate"]["features"]
    ]
    featureSetA = config["feature_sets"][
        config["correlation"][correlation]["A"]["feature_set"]
    ]["features"]
    featureSetB = config["feature_sets"][
        config["correlation"][correlation]["B"]["feature_set"]
    ]["features"]
    for feature in features:
        if feature[0] not in featureSetA:
            raise Exception(
                "Workflow correlation: cannot find feature in %s in feature set for set A!"
                % feature[0]
            )
        if feature[1] not in featureSetB:
            raise Exception(
                "Workflow correlation: cannot find feature in %s in feature set for set A!"
                % feature[0]
            )
    return features


def getCorrPlotLabs(correlation, lab):
    if "labs" in config["correlation"][correlation]["correlate"]:
        return config["correlation"][correlation]["correlate"]["labs"][lab]
    else:
        return ""


def getCorrPlotXLab(correlation):
    return getCorrPlotLabs(correlation, "x")


def getCorrPlotYLab(correlation):
    return getCorrPlotLabs(correlation, "y")
