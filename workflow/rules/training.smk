# training workflow

## check if positie/negatve and feature set have the same genome build


def getVariantSetGenomeBuild(variant_set):
    variantSet_conf = config["variants"][variant_set]
    genomeBuild = variantSet_conf["genome_build"]
    switcher = {"hg38": "hg19", "hg19": "hg38"}
    if "liftover" in variantSet_conf:
        return switcher[genomeBuild]
    else:
        return genomeBuild


for training_run, conf in config["training"].items():
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


rule getFolds:
    input:
        m="resources/ReMM.20171122.partition.mapping.tsv.gz",
        c="resources/{genomeBuild}/cytoBand.txt.gz",
    output:
        "resources/{genomeBuild}/folds.txt.gz",
    script:
        "../scripts/createFolds.py"


def getTrainingRunFolds(training_run):
    genomeBuild = getVariantSetGenomeBuild(
        config["training"][training_run]["positives"]
    )
    return "resources/%s/cytoBand.txt.gz" % genomeBuild


include: "training/parSMURF.smk"
