import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

rule get_Fantom5Perm_hg19:
    input:
        HTTP.remote(features['Fantom5Perm']['hg19']["url"], keep_local=False)
    output:
        "results/features/download/Fantom5Perm/hg19/Fantom5Perm.all.bed.gz"
    shell:
        """
        cat {input} | tail -n +2 | sort -k1,1 -k2,2n -k3,3n | bgzip -c > {output}
        """

rule get_Fantom5Perm_hg38:
    output:
        "results/features/download/Fantom5Perm/hg38/Fantom5Perm.all.bed.gz"
    params:
        url=features['Fantom5Perm']['hg38']["url"]
    shell:
        """
        curl {params.url} > {output}
        """

rule get_Fantom5Robust_hg19:
    input:
        HTTP.remote(features['Fantom5Robust']['hg19']["url"], keep_local=False)
    output:
        "results/features/download/Fantom5Robust/hg19/Fantom5Robust.all.bed.gz"
    shell:
        """
        cat {input} | tail -n +2 | sort -k1,1 -k2,2n -k3,3n | bgzip -c > {output}
        """

rule get_Fantom5Robust_hg38:
    output:
        "results/features/download/Fantom5Robust/hg38/Fantom5Robust.all.bed.gz"
    params:
        url=features['Fantom5Robust']['hg38']["url"]
    shell:
        """
        curl {params.url} | zcat | sort -k1,1 -k2,2n -k3,3n | bgzip -c  > {output}
        """

rule features_download_fantom_getFWD:
    input:
        "results/features/download/{feature_fantom}/hg38/{feature_fantom}.all.bed.gz"
    output:
        "results/features/download/{feature_fantom}/hg38/{feature_fantom}.fwd.bed.gz"
    wildcard_constraints:
        feature_fantom="(Fantom5Robust)|(Fantom5Perm)"
    shell:
        """
        zcat {input} | awk '{{if ($6 == "+") {{print $0}}}}' | \
        bgzip -c > {output}
        """

rule features_download_fantom_getREV:
    input:
        "results/features/download/{feature_fantom}/hg38/{feature_fantom}.all.bed.gz"
    output:
        "results/features/download/{feature_fantom}/hg38/{feature_fantom}.rev.bed.gz"
    wildcard_constraints:
        feature_fantom="(Fantom5Robust)|(Fantom5Perm)"
    shell:
        """
        zcat {input} | awk '{{if ($6 == "-") {{print $0}}}}' | \
        bgzip -c > {output}
        """
