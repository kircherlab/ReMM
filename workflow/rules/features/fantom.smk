rule get_Fantom5Perm:
    output:
        "results/features/download/Fantom5Perm/hg38/Fantom5Perm.all.bed.gz"
    params:
        url=features['Fantom5Perm']['hg38']["url"]
    shell:
        """
        curl {params.url} > {output}
        """

rule get_Fantom5Robust:
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