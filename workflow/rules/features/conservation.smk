# feature subworkflow for conservation scores


rule features_conservation_convertToBigWig_hg19:
    conda:
        "../../envs/default.yaml"
    input:
        "results/features/download/{feature}/hg19/{chr}.{extension}.wigFix.gz",
    output:
        "results/features/download/{feature}/hg19/{chr}.{extension}.bw.gz",
    log:
        temp(
            "results/logs/features/conservation/convertToBigWig_hg19.{feature}.{chr}.{extension}.log"
        ),
    shell:
        """
        bigWigToBedGraph {input} utils/{wildcards.genomeBuild}.chrom.sizes {output}
        """


rule features_conservation_PriPhyloP_hg19_download_process:
    conda:
        "../../envs/default.yaml"
    output:
        "results/features/download/priPhyloP/hg19/priPhyloP.{chr}.wig.gz",
    params:
        url=lambda wildcards: "%s%s.phyloP46way.primate.wigFix.gz"
        % (
            features["priPhyloP"]["hg19"]["url"],
            wildcards.chr,
        ),
    log:
        temp(
            "results/logs/features/conservation/PriPhyloP_hg19_download_process.{chr}.log"
        ),
    shell:
        """
        curl {params.url} > {output}
        """


rule features_conservation_PriPhyloP_hg38_download_process:
    conda:
        "../../envs/default.yaml"
    output:
        "results/features/download/priPhyloP/hg38/priPhyloP.all.wig.gz",
    params:
        url=features["priPhyloP"]["hg38"]["url"],
    log:
        temp("results/logs/features/conservation/PriPhyloP_hg38_download_process.log"),
    shell:
        """
        curl {params.url} | zcat | \
        awk -v 'put="F"' '{{
          if ($0 ~ /chrom/) {{
            if ($0 ~ /(\schrom=(chr([0-9]+|X|Y|M))\sstart=[0-9]+\sstep=[0-9]+)/) {{
              print $0;
              put="T";
            }} else {{
              put="F";
            }}
          }} else {{
            if (put=="T") {{
              print $0;
            }}
          }} }}' | bgzip -c > {output}
        """


rule features_conservation_PriPhastCons_hg19_download_process:
    conda:
        "../../envs/default.yaml"
    output:
        "results/features/download/priPhastCons/hg19/priPhastCons.{chr}.wig.gz",
    params:
        url=lambda wildcards: "%s%s.phastCons46way.primates.wigFix.gz"
        % (
            features["priPhastCons"]["hg19"]["url"],
            wildcards.chr,
        ),
    log:
        temp(
            "results/logs/features/conservation/PriPhastCons_hg19_download_process.{chr}.log"
        ),
    shell:
        """
        curl {params.url} > {output}
        """


rule features_conservation_PriPhastCons_hg38_download_process:
    conda:
        "../../envs/default.yaml"
    output:
        "results/features/download/priPhastCons/hg38/priPhastCons.all.wig.gz",
    params:
        url=features["priPhastCons"]["hg38"]["url"],
    log:
        temp(
            "results/logs/features/conservation/PriPhastCons_hg38_download_process.log"
        ),
    shell:
        """
        curl {params.url} | zcat | \
        awk -v 'put="F"' '{{
          if ($0 ~ /chrom/) {{
            if ($0 ~ /(\schrom=(chr([0-9]+|X|Y|M))\sstart=[0-9]+\sstep=[0-9]+)/) {{
              print $0;
              put="T";
            }} else {{
              put="F";
            }}
          }} else {{
            if (put=="T") {{
              print $0;
            }}
          }} }}' | bgzip -c > {output}
        """


rule features_conservation_VerPhyloP_hg38_download_process:
    conda:
        "../../envs/default.yaml"
    output:
        "results/features/download/verPhyloP/hg38/verPhyloP.{chr}.wig.gz",
    params:
        url=lambda wildcards: "%s%s.phyloP100way.wigFix.gz"
        % (
            features["verPhyloP"]["hg38"]["url"],
            wildcards.chr,
        ),
    log:
        temp(
            "results/logs/features/conservation/VerPhyloP_hg38_download_process.{chr}.log"
        ),
    shell:
        """
        curl {params.url} > {output}
        """


rule features_conservation_VerPhyloP_hg19_download_process:
    conda:
        "../../envs/default.yaml"
    output:
        "results/features/download/verPhyloP/hg19/verPhyloP.{chr}.wig.gz",
    params:
        url=lambda wildcards: "%s%s.phyloP46way.wigFix.gz"
        % (
            features["verPhyloP"]["hg19"]["url"],
            wildcards.chr,
        ),
    log:
        temp(
            "results/logs/features/conservation/VerPhyloP_hg19_download_process.{chr}.log"
        ),
    shell:
        """
        curl {params.url} > {output}
        """


rule features_conservation_VerPhastCons_hg38_download_process:
    conda:
        "../../envs/default.yaml"
    output:
        "results/features/download/verPhastCons/hg38/verPhastCons.{chr}.wig.gz",
    params:
        url=lambda wildcards: "%s%s.phastCons100way.wigFix.gz"
        % (
            features["verPhastCons"]["hg38"]["url"],
            wildcards.chr,
        ),
    log:
        temp(
            "results/logs/features/conservation/VerPhastCons_hg38_download_process.{chr}.log"
        ),
    shell:
        """
        curl {params.url} > {output}
        """


rule features_conservation_VerPhastCons_hg19_download_process:
    conda:
        "../../envs/default.yaml"
    output:
        "results/features/download/verPhastCons/hg19/verPhastCons.{chr}.wig.gz",
    params:
        url=lambda wildcards: "%s%s.phastCons46way.wigFix.gz"
        % (
            features["verPhastCons"]["hg19"]["url"],
            wildcards.chr,
        ),
    log:
        temp(
            "results/logs/features/conservation/VerPhastCons_hg19_download_process.{chr}.log"
        ),
    shell:
        """
        curl {params.url} > {output}
        """


rule features_conservation_MamPhastCons_hg38_download_process:
    conda:
        "../../envs/default.yaml"
    output:
        "results/features/download/mamPhastCons/hg38/mamPhastCons.{chr}.wig.gz",
    params:
        url=lambda wildcards: "%s%s.phastCons30way.wigFix.gz"
        % (
            features["mamPhastCons"]["hg38"]["url"],
            wildcards.chr,
        ),
    log:
        temp(
            "results/logs/features/conservation/MamPhastCons_hg38_download_process.{chr}.log"
        ),
    shell:
        """
        curl {params.url} > {output}
        """


rule features_conservation_MamPhastCons_hg19_download_process:
    conda:
        "../../envs/default.yaml"
    output:
        "results/features/download/mamPhastCons/hg19/mamPhastCons.{chr}.wig.gz",
    params:
        url=lambda wildcards: "%s%s.phastCons46way.placental.wigFix.gz"
        % (
            features["mamPhastCons"]["hg19"]["url"],
            wildcards.chr,
        ),
    log:
        temp(
            "results/logs/features/conservation/MamPhastCons_hg19_download_process.{chr}.log"
        ),
    shell:
        """
        curl {params.url} > {output}
        """


rule features_conservation_MamPhyloP_hg38_download_process:
    conda:
        "../../envs/default.yaml"
    output:
        "results/features/download/mamPhyloP/hg38/mamPhyloP.{chr}.wig.gz",
    params:
        url=lambda wildcards: "%s%s.phyloP30way.wigFix.gz"
        % (
            features["mamPhyloP"]["hg38"]["url"],
            wildcards.chr,
        ),
    log:
        temp(
            "results/logs/features/conservation/MamPhyloP_hg38_download_process.{chr}.log"
        ),
    shell:
        """
        curl {params.url} > {output}
        """


rule features_conservation_MamPhyloP_hg19_download_process:
    conda:
        "../../envs/default.yaml"
    output:
        "results/features/download/mamPhyloP/hg19/mamPhyloP.{chr}.wig.gz",
    params:
        url=lambda wildcards: "%s%s.phyloP46way.placental.wigFix.gz"
        % (
            features["mamPhyloP"]["hg19"]["url"],
            wildcards.chr,
        ),
    log:
        temp(
            "results/logs/features/conservation/MamPhyloP_hg19_download_process.{chr}.log"
        ),
    shell:
        """
        curl {params.url} > {output}
        """


# TODO hard coded path!


rule features_conservation_GERP_hg38_process:
    conda:
        "../../envs/default.yaml"
    input:
        features["GerpRS"]["hg38"]["url"]
    output:
        "results/features/download/gerpElement/hg38/gerpElement.all.bed.gz",
    log:
        temp("results/logs/features/conservation/GERP_hg38_process.log"),
    shell:
        """
        zcat {input} | \
        awk -v 'OFS=\\t' '{{print "chr"$0}}' | bgzip -c > {output}
        """


# cut -f 1,2,76-77 | \
# uniq | grep -v NA | \
# awk -v 'OFS=\\t' 'BEGIN {{chr=$1;start=$2-1;end=$2;rs=$3;p=$4}} {{ if (chr==$1 && end == ($2-1) && rs==$3 && p==$4) {{end=$2;}} else {{print "chr"chr,start,end,rs,p; chr=$1;start=($2-1);end=$2;rs=$3;p=$4}} }}' | \
# tail -n +2 | bgzp -c >
