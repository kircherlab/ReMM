################################
#### training sub workflow  ####
################################


"""
Results will be saved in results/training/<training_run>
"""

include: "training_defs.smk"


rule getFolds:
    input:
        m="resources/ReMM.20171122.partition.mapping.tsv.gz",
        c="resources/{genomeBuild}/cytoBand.txt.gz",
    output:
        "resources/{genomeBuild}/folds.txt.gz",
    script:
        "../scripts/createFolds.py"


include: "training/parSMURF.smk"
