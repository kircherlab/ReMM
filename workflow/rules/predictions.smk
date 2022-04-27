
include: "predictions_defs.smk"

checkpoint predictions_createInputData:
    input:
        lambda wc: expand(
            "results/annotation/{{variant_set}}/{{variant_set}}.{feature_set}.{missing_value}.sorted.tsv.gz",
            feature_set=config["training"][wc.training]["feature_set"],
            missing_value=config["training"][wc.training]["missing_value"],
        ),
    output:
        temp(
            "results/predictions/{training}/{variant_set}/input/parsmurf/parsmurf.data.{split}.txt"
        ),
    params:
        lines=100000,
        prefix=(
            "results/predictions/{training}/{variant_set}/input/parsmurf/parsmurf.data."
        ),
        suffix=".txt",
    shell:
        """
        split --additional-suffix={params.suffix} -a 4 -l {params.lines} <(
            zcat {input} | egrep -v "^CHR\sPOSITION\sID" | cut -f 4-
        ) {params.prefix}
        """


# rule predictions_createInputData:
#     input:
#         lambda wc: expand(
#             "results/annotation/{{variant_set}}/{{variant_set}}.{feature_set}.{missing_value}.sorted.tsv.gz",
#             feature_set=config["training"][wc.training]["feature_set"],
#             missing_value=config["training"][wc.training]["missing_value"],
#         ),
#     output:
#         temp(
#             "results/predictions/{training}/{variant_set}/input/parsmurf/parsmurf.data.txt"
#         ),
#     shell:
#         """
#         zcat {input} | egrep -v "^CHR\sPOSITION\sID" | cut -f 4- > {output}
#         """


rule predictions_createInputLabels:
    input:
        "results/predictions/{training}/{variant_set}/input/parsmurf/parsmurf.data.{split}.txt",
    output:
        temp(
            "results/predictions/{training}/{variant_set}/input/parsmurf/parsmurf.labels.{split}.txt"
        ),
    shell:
        """
        cat {input} | awk '{{print 1}}' > {output}
        """


rule predictions_parSMURF_conf:
    input:
        data="results/predictions/{training}/{variant_set}/input/parsmurf/parsmurf.data.{split}.txt",
        models=expand(
            "results/training/{{training}}/predictions/models/{model_number}.out.forest",
            model_number=list(range(0, 100)),
        ),
        labels="results/predictions/{training}/{variant_set}/input/parsmurf/parsmurf.labels.{split}.txt",
        scaffold="resources/scaffold.json",
    output:
        config=temp(
            "results/predictions/{training}/{variant_set}/input/model/parsmurf.config.{split}.json"
        ),
    params:
        predictions="results/predictions/{training}/{variant_set}/predictions/parsmurf/predictions.{split}.txt",
        name="{training}_{variant_set}",
        models="results/training/{training}/predictions/models",
        seed=lambda wc: config["training"][wc.training]["config"]["seed"],
        mode="predict",
        ensThrd="30",
    script:
        "../scripts/generateParsmurfConfig.py"


rule predictions_parSMURF_run:
    input:
        data="results/predictions/{training}/{variant_set}/input/parsmurf/parsmurf.data.{split}.txt",
        config="results/predictions/{training}/{variant_set}/input/model/parsmurf.config.{split}.json",
        labels="results/predictions/{training}/{variant_set}/input/parsmurf/parsmurf.labels.{split}.txt",
        models=expand(
            "results/training/{{training}}/predictions/models/{model_number}.out.forest",
            model_number=list(range(0, 100)),
        ),
    output:
        temp(
            "results/predictions/{training}/{variant_set}/predictions/parsmurf/predictions.{split}.txt"
        ),
    shell:
        """
        workflow/bin/parSMURF1 --cfg {input.config}
        """


rule predictions_parSMURF_aggregate:
    input:
        predictions=aggregate_PredictedScores,
        data=aggregate_Data,
    output:
        "results/predictions/{training}/{variant_set}/predictions/parsmurf/predictions.txt",
    shell:
        """
        cat {input.predictions} > {output};
        """


rule predictionsparSMURF_combine_labels:
    input:
        predictions=(
            "results/predictions/{training}/{variant_set}/predictions/parsmurf/predictions.txt",
        ),
        data=lambda wc: expand(
            "results/annotation/{{variant_set}}/{{variant_set}}.{feature_set}.{missing_value}.sorted.tsv.gz",
            feature_set=config["training"][wc.training]["feature_set"],
            missing_value=config["training"][wc.training]["missing_value"],
        ),
    output:
        "results/predictions/{training}/{variant_set}/predictions/predictions.with_IDs.tsv.gz",
    shell:
        """
        (
            echo -e "CHROM\\tPOS\\tID\\tSCORE";
            paste \
            <(
                zcat {input.data} | egrep -v "^CHR\sPOSITION\sID" | awk -v 'OFS=\\t' '{{print $1,$2,$3}}';
            ) \
            <(cat {input.predictions} | cut -f 1 ) | \
            sort -k 1,1 -k2,2n;
        ) | \
        bgzip -c > {output}
        """
