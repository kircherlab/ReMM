rule get1KG:
    output:
        bed=temp("results/features/hg38/1KG/1KG.{chr}.bed"),
        bed_gz="results/features/hg38/1KG/1KG.{chr}.bed.gz",
    params:
        chr="{chr}",
    conda:
        "../../envs/GCContent.yaml"
    shell:
        """
        set +o pipefail;
        curl -s --connect-timeout 540 --retry 20 --retry-all-errors http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20200515_EBI_Freebayescalls/ALL.{params.chr}.freebayes.20200518.snps_indels.NYhc.GRCh38.vcf.gz | \
        zcat  | cut -f 1-9 | \
        ruby scripts/rareVariantFractionInWindow.rb 500 0.005 {output.bed};
        cat {output.bed}| bgzip -c > {output.bed_gz}
        """


# http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
# http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL{params.chr}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz  | \
#        zcat  | cut -f 1-9 | \
#        ruby scripts/rareVariantFractionInWindow.rb 500 0.005 {output.bed};
#        cat {output.bed}| bgzip -c > {output.bed_gz}

# rule convert_feature_1kG:
#    input:
#       "input/features/{genomeBuild}/1KG/1KG.{chr}.vcf.gz"
#    output:
#        "input/features/{genomeBuild}/1KG/1KG.{chr}.bed.gz"
#    conda:
#        "../../envs/GCContent.yml"
#    shell:
#        """
#        ruby scripts/rareVariantFractionInWindow.rb 500 0.005 {input} | bgzip -c > {output}
#        """
