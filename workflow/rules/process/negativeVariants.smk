# TODO Hard coded path


rule getNegativeVariants:
    input:
        "/fast/groups/ag_kircher/CADD/cadd_v1.3/training_data/GRCh38/humanDerived/annotated/SNVs.vcf.gz",
    output:
        "results/variants/hg38/SNVs.hg38.negative.vcf.gz",
    shell:
        """
        (
            echo -e "##fileformat=VCFv4.3\\n#CHROM\\tPOS\\tID\tREF\\tALT\\tQUAL\\tFILTER\\tINFO";
            zcat {input}  | awk -v 'OFS=\\t' '{{print $1="chr" $1,$2,$3,$4,$5,".","PASS","."}}';
        ) | bgzip -c > {output}
        """
