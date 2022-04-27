def getRegionFileForVariantGeneration(name):
    """
    Helper to get the correct file for regions or variants
    """
    if config["variant_generation"][name]["type"] == "random":
        return "results/variant_generation/{name}/hg38/random_full.{regions_or_variants}.bed.gz"
    elif config["variant_generation"][name]["type"] == "regions":
        return "results/variant_generation/{name}/hg38/regions_full.{regions_or_variants}.bed.gz"
    elif config["variant_generation"][name]["type"] == "variants":
        return "results/variant_generation/{name}/hg38/variants_full.{regions_or_variants}.bed.gz"
    raise Exception("Unknown type of variant generation name: %s" % name)
