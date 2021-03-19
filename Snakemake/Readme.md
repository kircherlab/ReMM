## Project structure

The folder Snakemake contains several subfolders:

- config: configuration files for download and processing of annotations as well as configuration files for running parSMURF
- rules: snakemake rules separated into subfolders according to the workflow steps
- env: a yaml file for creating the main working environment 'ReMM' as well as files for some additional environments required by certain rules
- scripts: (external) scripts used by snakemake rules
- utils: diverse reference files and scaffolds used during the workflow. 

Some directories contain more detailed information on files contained.






The Snakemake workflow consists of XX different rules, many of which are excecuted multiple times. 
In the following, we briefly discuss what most of the rules thus and what changes you need to do, to adapt the wokflow to your data. 


Each feature is downloaded and preprocessed by a separate Snakemake rule (rules/features). This modularization is needed due to different data sources and due to the non-identical preprocessing steps. The download  is managed by a config file (config/featuresConfig38.json) that contains the download links and the nessasary meta information for the processing steps.



Absolute paths:
/fast/groups/ag_kircher/CADD/dependencies/annotations/gerp/gerp2_elements_hg38_MAM.bg.gz
/fast/groups/ag_kircher/ReMM/MA_Lusi/Snakemake/scripts/createIntervals.py
/fast/groups/ag_kircher/ReMM/MA_Lusi/Snakemake/scripts/numTFBSConserved.py
/fast/groups/ag_kircher/ReMM/MA_Lusi/Snakemake/scripts/createIntervals.py
/fast/work/groups/ag_kircher/ReMM/MA_Lusi/Snakemake/scripts/generateParsmurfConfig.py
/fast/work/groups/ag_kircher/ReMM/MA_Lusi/Snakemake/scripts/createParsmurfInput.py
/fast/groups/ag_kircher/ReMM/MA_Lusi/Snakemake/scripts/createFolds.py
/fast/groups/ag_kircher/CADD/cadd_v1.3/training_data/GRCh38/humanDerived/annotated/SNVs.vcf.gz
/fast/work/groups/ag_kircher/ReMM/ReMM/data/variants/RegulatoryMendelianMutations/GRCh37/SNVs.all.20160109.vcf.gz
"/fast/groups/ag_kircher/ReMM/MA_Lusi/Snakemake/scripts/liftOver.py

in
Snakemake/scripts/createHyperSmurfInput.py
Snakemake/scripts/createParsmurfInput.py
Snakemake/scripts/liftOver.py
Snakemake/scripts/numTFBSConserved.py
