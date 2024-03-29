rule get_feature_1KG_19:
    output:
        bed=temp("results/features/hg19/1KG/1KG.{chr}.bed"),
        bed_gz="results/features/hg19/1KG/1KG.{chr}.bed.gz",
    params:
        chr="{chr}",
    conda:
        "../../envs/GCContent.yaml"
    shell:
        """
        curl http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.{params.chr}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz  | \
        zcat  | cut -f 1-9 | \
        ruby scripts/rareVariantFractionInWindow.rb 500 0.005 {output.bed};
        cat {output.bed}| bgzip -c > {output.bed_gz}
        """
