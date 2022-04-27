def aggregate_PredictedScores(wc):
    checkpoint_output = checkpoints.predictions_createInputData.get(
        **wc, split="aaaa"
    ).output[0]
    variables = glob_wildcards(
        os.path.join(os.path.dirname(checkpoint_output), "parsmurf.data.{i}.txt")
    ).i
    variables.sort()
    return expand(
        "results/predictions/{training}/{variant_set}/predictions/parsmurf/predictions.{split}.txt",
        training=wc.training,
        variant_set=wc.variant_set,
        split=variables,
    )


def aggregate_Data(wc):
    checkpoint_output = checkpoints.predictions_createInputData.get(
        **wc, split="aaaa"
    ).output[0]
    return expand(
        "results/predictions/{training}/{variant_set}/input/parsmurf/parsmurf.data.{split}.txt",
        training=wc.training,
        variant_set=wc.variant_set,
        split=glob_wildcards(
            os.path.join(os.path.dirname(checkpoint_output), "parsmurf.data.{i}.txt")
        ).i,
    )
