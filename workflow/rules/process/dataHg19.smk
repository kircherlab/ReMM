# TODO Script getDataHg19.py is missing
# TODO Hard coded paths


rule getDataHg19:
    input:
        p="/fast/work/groups/ag_kircher/ReMM/ReMM/data/variants/RegulatoryMendelianMutations/GRCh37/SNVs.all.20160109.vcf.gz",
        na="/fast/work/groups/ag_kircher/ReMM/ReMM/data/training/annotated/HumanDerived/GRCh37/SNVs.noncoding.20180118.ReMM.20180119.tsv.gz",
        pa="/fast/work/groups/ag_kircher/ReMM/ReMM/data/training/annotated/RegulatoryMendelianMutations/GRCh37/SNVs_noncoding_20160101.ReMM.20171122.tsv.gz",
    output:
        p="results/variants/hg19/SNVs.hg19.positive.annotated.tsv.gz",
        n="results/variants/hg19/SNVs.hg19.negative.annotated.tsv.gz",
    script:
        "../../scripts/getDataHg19.py"


# TODO Script getRemmScore.py is missing


rule getReMM:
    output:
        "results/predictions/hg19/remm_with_ps_predictions.csv",
    script:
        "../../scripts/getRemmScore.py"
