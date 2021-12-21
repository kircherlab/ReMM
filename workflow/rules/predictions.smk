rule predictions_createInputData:
    input:
        lambda wc: expand(
            "results/annotation/{{variant_set}}/{{variant_set}}.{feature_set}.{missing_value}.sorted.tsv.gz",
            feature_set=config["training"][wc.training]["feature_set"],
            missing_value=config["training"][wc.training]["missing_value"],
        ),
    output:
        temp(
            "results/predictions/{training}/{variant_set}/input/parsmurf/parsmurf.data.txt"
        ),
    shell:
        """
        zcat {input} | egrep -v "^CHR\sPOSITION\sID" | cut -f 4- > {output}
        """


rule predictions_createInputLabels:
    input:
        lambda wc: expand(
            "results/annotation/{{variant_set}}/{{variant_set}}.{feature_set}.{missing_value}.sorted.tsv.gz",
            feature_set=config["training"][wc.training]["feature_set"],
            missing_value=config["training"][wc.training]["missing_value"],
        ),
    output:
        temp(
            "results/predictions/{training}/{variant_set}/input/parsmurf/parsmurf.labels.txt"
        ),
    shell:
        """
        zcat {input} | egrep -v "^CHR\sPOSITION\sID" | awk '{{print 1}}' > {output}
        """


rule predictions_parSMURF_conf:
    input:
        data="results/predictions/{training}/{variant_set}/input/parsmurf/parsmurf.data.txt",
        models=("results/training/{training}/predictions/models/0.out.forest"),
        labels=(
            "results/predictions/{training}/{variant_set}/input/parsmurf/parsmurf.labels.txt"
        ),
        scaffold="resources/scaffold.json",
    output:
        config="results/predictions/{training}/{variant_set}/input/model/parsmurf.config.json",
    params:
        predictions=(
            "results/predictions/{training}/{variant_set}/predictions/parsmurf/predictions.txt"
        ),
        name="{training}_{variant_set}",
        models=("results/training/{training}/predictions/models"),
        seed=lambda wc: config["training"][wc.training]["config"]["seed"],
        mode="predict",
        ensThrd="30",
    script:
        "../scripts/generateParsmurfConfig.py"


rule predictions_parSMURF_run:
    input:
        data="results/predictions/{training}/{variant_set}/input/parsmurf/parsmurf.data.txt",
        config="results/predictions/{training}/{variant_set}/input/model/parsmurf.config.json",
        labels=(
            "results/predictions/{training}/{variant_set}/input/parsmurf/parsmurf.labels.txt"
        ),
        models=("results/training/{training}/predictions/models/0.out.forest"),
    output:
        "results/predictions/{training}/{variant_set}/predictions/parsmurf/predictions.txt",
    shell:
        """
        workflow/bin/parSMURF1 --cfg {input.config}
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
