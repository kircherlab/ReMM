###############################################
# Create a final whole genome wide score
###############################################


rule predict_combineScores:
    input:
        predictions=getPredictedScoresPerInterval,
    output:
        prediction="results/scores/{build}/{build_version}/{release}/score/{build}.{build_version}.{seed}.biased.tsv.gz",
    shell:
        """
        export LC_ALL=C;
        zcat {input.predictions} | sort -k1,1 -k2,2n | bgzip -c > {output.prediction};
