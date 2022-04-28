import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()


rule features_fantom_get_Fantom5Perm_hg19:
    conda:
        "../../envs/default.yml"
    input:
        HTTP.remote(features["Fantom5Perm"]["hg19"]["url"], keep_local=False),
    output:
        "results/features/download/Fantom5Perm/hg19/Fantom5Perm.all.bed.gz",
    log:
        temp("results/logs/features/fantom/get_Fantom5Perm_hg19.log"),
    shell:
        """
        cat {input} | tail -n +2 | sort -k1,1 -k2,2n -k3,3n | bgzip -c > {output}
        """


rule features_fantom_get_Fantom5Perm_hg38:
    conda:
        "../../envs/default.yml"
    output:
        "results/features/download/Fantom5Perm/hg38/Fantom5Perm.all.bed.gz",
    params:
        url=features["Fantom5Perm"]["hg38"]["url"],
    log:
        temp("results/logs/features/fantom/get_Fantom5Perm_hg38.log"),
    shell:
        """
        curl {params.url} > {output}
        """


rule features_fantom_get_Fantom5Robust_hg19:
    conda:
        "../../envs/default.yml"
    input:
        HTTP.remote(features["Fantom5Robust"]["hg19"]["url"], keep_local=False),
    output:
        "results/features/download/Fantom5Robust/hg19/Fantom5Robust.all.bed.gz",
    log:
        temp("results/logs/features/fantom/Fantom5Robust_hg19.log"),
    shell:
        """
        cat {input} | tail -n +2 | sort -k1,1 -k2,2n -k3,3n | bgzip -c > {output}
        """


rule features_fantom_get_Fantom5Robust_hg38:
    conda:
        "../../envs/default.yml"
    output:
        "results/features/download/Fantom5Robust/hg38/Fantom5Robust.all.bed.gz",
    params:
        url=features["Fantom5Robust"]["hg38"]["url"],
    log:
        temp("results/logs/features/fantom/Fantom5Robust_hg38.log"),
    shell:
        """
        curl {params.url} | zcat | sort -k1,1 -k2,2n -k3,3n | bgzip -c  > {output}
        """


rule features_fantom_download_fantom_getFWD:
    conda:
        "../../envs/default.yml"
    input:
        "results/features/download/{feature_fantom}/hg38/{feature_fantom}.all.bed.gz",
    output:
        "results/features/download/{feature_fantom}/hg38/{feature_fantom}.fwd.bed.gz",
    wildcard_constraints:
        feature_fantom="(Fantom5Robust)|(Fantom5Perm)",
    log:
        temp("results/logs/features/fantom/download_fantom_getFWD.{feature_fantom}.log"),
    shell:
        """
        zcat {input} | awk '{{if ($6 == "+") {{print $0}}}}' | \
        bgzip -c > {output}
        """


rule features_fantom_download_fantom_getREV:
    conda:
        "../../envs/default.yml"
    input:
        "results/features/download/{feature_fantom}/hg38/{feature_fantom}.all.bed.gz",
    output:
        "results/features/download/{feature_fantom}/hg38/{feature_fantom}.rev.bed.gz",
    wildcard_constraints:
        feature_fantom="(Fantom5Robust)|(Fantom5Perm)",
    log:
        temp("results/logs/features/fantom/download_fantom_getREV.{feature_fantom}.log"),
    shell:
        """
        zcat {input} | awk '{{if ($6 == "-") {{print $0}}}}' | \
        bgzip -c > {output}
        """
