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
