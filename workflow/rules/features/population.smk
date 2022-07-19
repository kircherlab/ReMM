import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()


rule features_1KG_hg38_download:
    conda:
        "../../envs/default.yml"
    output:
        vcf=temp("results/features/download/1KG/hg38/1KG.{chr}.vcf.gz"),
    params:
        url="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20200515_EBI_Freebayescalls/ALL.{chr}.freebayes.20200518.snps_indels.NYhc.GRCh38.vcf.gz",
    log:
        "results/logs/features/population/1KG_hg38_download.{chr}.log",
    shell:
        """
        wget --timeout=10000 --tries=0 --continue -O {output.vcf} {params.url} &> {log};
        #curl -s --connect-timeout 540 --retry 20 --retry-all-errors; 
        """


rule features_1KG_hg38_process:
    input:
        vcf="results/features/download/1KG/hg38/1KG.{chr}.vcf.gz",
    output:
        bed=temp("results/features/download/1KG/hg38/1KG.{chr}.bed"),
        bed_gz="results/features/download/1KG/hg38/1KG.{chr}.bed.gz",
    conda:
        "../../envs/ruby.yaml"
    log:
        "results/logs/features/population/1KG_hg38_process.{chr}.log",
    shell:
        """
        zcat {input.vcf} | cut -f 1-9 | ruby workflow/scripts/rareVariantFractionInWindow.rb 500 0.005 {output.bed} &> {log};
        cat {output.bed} | bgzip -c > {output.bed_gz} 2>> {log};
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
#        "../../envs/GCContent.yaml"
#    shell:
#        """
#        ruby scripts/rareVariantFractionInWindow.rb 500 0.005 {input} | bgzip -c > {output}
#        """
