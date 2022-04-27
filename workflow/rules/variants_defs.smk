def getVariantSetGenomeBuild(variant_set):
    """
    Getter of genomeBuild taking liftover into account.
    """
    variantSet_conf = config["variants"][variant_set]
    genomeBuild = variantSet_conf["genome_build"]
    switcher = {"hg38": "hg19", "hg19": "hg38"}
    if "processing" in variantSet_conf and "liftover" in variantSet_conf["processing"]:
        return switcher[genomeBuild]
    else:
        return genomeBuild


def getVariantsInput(variant_set, step, idx=False):
    """
    input variants are different for each process step becaus ethey are flexible.
    this trys to collect the correct input. The following steps are possible only in that order:
    liftover > jannovar > filters (bcftools > downsample) > annotate
    """

    add = ".tbi" if idx else ""
    variant_set_config = config["variants"][variant_set]
    if variant_set_config["type"] == "file":
        output = expand(
            "{file}{add}", file=variant_set_config["properties"]["file"], add=add
        )
    elif variant_set_config["type"] == "generation":
        output = expand(
            "results/variant_generation/{name}/{genomeBuild}/{name}.vcf.gz{add}",
            name=variant_set_config["properties"]["name"],
            genomeBuild=variant_set_config["genome_build"],
            add=add,
        )
    else:
        raise Exception("Unknown variant type %s" % variant_set_config["type"])
    if step == "liftover":
        return output
    # liftover
    if (
        "processing" in variant_set_config
        and "liftover" in variant_set_config["processing"]
    ):
        output = expand(
            "results/variants/{variant_set}/liftover/{variant_set}.vcf.gz{add}",
            variant_set=variant_set,
            add=add,
        )
    # jannovar
    if step == "jannovar":
        return output
    if (
        "processing" in variant_set_config
        and "jannovar" in variant_set_config["processing"]
    ):
        output = expand(
            "results/variants/{variant_set}/jannovar/{variant_set}.vcf.gz{add}",
            variant_set=variant_set,
            add=add,
        )
    # bcftools
    if step == "bcftools":
        return output
    if (
        "processing" in variant_set_config
        and "filters" in variant_set_config["processing"]
    ):
        if "bcftools" in variant_set_config["processing"]["filters"]:
            output = expand(
                "results/variants/{variant_set}/bcftools/{variant_set}.vcf.gz{add}",
                variant_set=variant_set,
                add=add,
            )
    # downsample
    if step == "downsample":
        return output
    if (
        "processing" in variant_set_config
        and "filters" in variant_set_config["processing"]
    ):
        if "downsample" in variant_set_config["processing"]["filters"]:
            output = expand(
                "results/variants/{variant_set}/downsample/{variant_set}.vcf.gz{add}",
                variant_set=variant_set,
                add=add,
            )
    if step == "annotate":
        return output

    raise Exception("Unknown variant processing step %s" % step)
