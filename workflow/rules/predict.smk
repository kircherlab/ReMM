###############################################
# Create a final whole genome wide score
###############################################


rule predict_getNotNRegions:
    input:
        lambda wc: config["global_files"]["genome_builds"][wc.genomeBuild]["reference"],
    output:
        "resources/{genomeBuild}/notNRegions.bed",
    shell:
        """
        python workflow/scripts/notNRegionsOfGenome.py --input {input} --output {output} 
        """


rule predict_exttact_features:
    input:
        "results/prediction/training_run}/input/genome/{split}.unsorted.tsv.gz",
    output:
        "results/prediction/training_run}/input/genome/{split}.sorted.tsv.gz",
    shell:
        """
        python workflow/getAnnotationByinterval.py --input {input.vcf} --output {output} \
        --regions {input.regions}
        """


rule predict_sort_features:
    input:
        "results/prediction/training_run}/input/genome/{split}.unsorted.tsv.gz",
    output:
        "results/prediction/training_run}/input/genome/{split}.sorted.tsv.gz",
    shell:
        """
        python workflow/sortAnnotationFile.py --input {input} --output {output}
        """


# rule predict_combineScores:
#     input:
#         predictions=getPredictedScoresPerInterval,
#     output:
#         prediction="results/scores/{build}/{build_version}/{release}/score/{build}.{build_version}.{seed}.biased.tsv.gz",
#     shell:
#         """
#         export LC_ALL=C;
#         zcat {input.predictions} | sort -k1,1 -k2,2n | bgzip -c > {output.prediction};
#         """
