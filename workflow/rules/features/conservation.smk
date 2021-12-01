rule convertToBigWig_hg19:
    input:
        "results/features/download/{feature}/hg19/{chr}.{extension}.wigFix.gz",
    output:
        "results/features/download/{feature}/hg19/{chr}.{extension}.bw.gz",
    shell:
        """
        bigWigToBedGraph {input} utils/{wildcards.genomeBuild}.chrom.sizes {output}
        """


rule features_PriPhyloP_hg19_download_process:
    output:
        "results/features/download/priPhyloP/hg19/priPhyloP.{chr}.wig.gz",
    params:
        url=lambda wildcards: "%s%s.phyloP46way.primate.wigFix.gz"
        % (
            features["priPhyloP"]["hg19"]["url"],
            wildcards.chr,
        ),
    shell:
        """
        curl {params.url} > {output}
        """


rule features_PriPhyloP_hg38_download_process:
    output:
        "results/features/download/priPhyloP/hg38/priPhyloP.all.wig.gz",
    params:
        url=features["priPhyloP"]["hg38"]["url"],
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


rule features_PriPhastCons_hg19_download_process:
    output:
        "results/features/download/priPhastCons/hg19/priPhastCons.{chr}.wig.gz",
    params:
        url=lambda wildcards: "%s%s.phastCons46way.primate.wigFix.gz"
        % (
            features["priPhastCons"]["hg19"]["url"],
            wildcards.chr,
        ),
    shell:
        """
        curl {params.url} > {output}
        """


rule features_PriPhastCons_hg38_download_process:
    output:
        "results/features/download/priPhastCons/hg38/priPhastCons.all.wig.gz",
    params:
        url=features["priPhastCons"]["hg38"]["url"],
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


rule features_VerPhyloP_hg38_download_process:
    output:
        "results/features/download/verPhyloP/hg38/verPhyloP.{chr}.wig.gz",
    params:
        url=lambda wildcards: "%s%s.phyloP100way.wigFix.gz"
        % (
            features["verPhyloP"]["hg38"]["url"],
            wildcards.chr,
        ),
    shell:
        """
        curl {params.url} > {output}
        """


rule features_VerPhyloP_hg19_download_process:
    output:
        "results/features/download/verPhyloP/hg19/verPhyloP.{chr}.wig.gz",
    params:
        url=lambda wildcards: "%s%s.phyloP46way.wigFix.gz"
        % (
            features["verPhyloP"]["hg19"]["url"],
            wildcards.chr,
        ),
    shell:
        """
        curl {params.url} > {output}
        """


rule features_VerPhastCons_hg38_download_process:
    output:
        "results/features/download/verPhastCons/hg38/verPhastCons.{chr}.wig.gz",
    params:
        url=lambda wildcards: "%s%s.phastCons100way.wigFix.gz"
        % (
            features["verPhastCons"]["hg38"]["url"],
            wildcards.chr,
        ),
    shell:
        """
        curl {params.url} > {output}
        """


rule features_VerPhastCons_hg19_download_process:
    output:
        "results/features/download/verPhastCons/hg19/verPhastCons.{chr}.wig.gz",
    params:
        url=lambda wildcards: "%s%s.phastCons46way.wigFix.gz"
        % (
            features["verPhastCons"]["hg19"]["url"],
            wildcards.chr,
        ),
    shell:
        """
        curl {params.url} > {output}
        """


rule features_MamPhastCons_hg38_download_process:
    output:
        "results/features/download/mamPhastCons/hg38/mamPhastCons.{chr}.wig.gz",
    params:
        url=lambda wildcards: "%s%s.phastCons30way.wigFix.gz"
        % (
            features["mamPhastCons"]["hg38"]["url"],
            wildcards.chr,
        ),
    shell:
        """
        curl {params.url} > {output}
        """


rule features_MamPhastCons_hg19_download_process:
    output:
        "results/features/download/mamPhastCons/hg19/mamPhastCons.{chr}.wig.gz",
    params:
        url=lambda wildcards: "%s%s.phastCons46way.placental.wigFix.gz"
        % (
            features["mamPhastCons"]["hg19"]["url"],
            wildcards.chr,
        ),
    shell:
        """
        curl {params.url} > {output}
        """


rule features_MamPhyloP_hg38_download_process:
    output:
        "results/features/download/mamPhyloP/hg38/mamPhyloP.{chr}.wig.gz",
    params:
        url=lambda wildcards: "%s%s.phyloP30way.wigFix.gz"
        % (
            features["mamPhyloP"]["hg38"]["url"],
            wildcards.chr,
        ),
    shell:
        """
        curl {params.url} > {output}
        """

rule features_MamPhyloP_hg19_download_process:
    output:
        "results/features/download/mamPhyloP/hg19/mamPhyloP.{chr}.wig.gz",
    params:
        url=lambda wildcards: "%s%s.phastCons46way.placental.wigFix.gz"
        % (
            features["mamPhyloP"]["hg19"]["url"],
            wildcards.chr,
        ),
    shell:
        """
        curl {params.url} > {output}
        """


# TODO hard coded path!


rule features_GERP_hg38_process:
    input:
        "/fast/groups/ag_kircher/CADD/dependencies/annotations/gerp/gerp2_elements_hg38_MAM.bg.gz",
    output:
        "results/features/download/gerpElement/hg38/gerpElement.all.bed.gz",
    shell:
        """
        zcat {input} | \
        awk -v 'OFS=\\t' '{{print "chr"$0}}' | bgzip -c > {output}
        """


# cut -f 1,2,76-77 | \
# uniq | grep -v NA | \
# awk -v 'OFS=\\t' 'BEGIN {{chr=$1;start=$2-1;end=$2;rs=$3;p=$4}} {{ if (chr==$1 && end == ($2-1) && rs==$3 && p==$4) {{end=$2;}} else {{print "chr"chr,start,end,rs,p; chr=$1;start=($2-1);end=$2;rs=$3;p=$4}} }}' | \
# tail -n +2 | bgzp -c >
