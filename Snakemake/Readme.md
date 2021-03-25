## Project workflow
### Folder structure

The folder Snakemake contains several subfolders:

- config: configuration files for downloading and processing of annotations as well as configuration files for running parSMURF
- rules: snakemake rules separated into subfolders according to the workflow steps
- env: a yaml file for creating the main working environment 'ReMM' as well as files for some additional environments required by certain rules
- scripts: (external) scripts used by snakemake rules
- utils: diverse reference files and scaffolds used during the workflow
- input: files that serve as input for cetain rules (raw fetaure files, variants)
- output: computed files - output of rules (VCF feature files, training sets, predictions)

Some directories contain more detailed information on files contained.

### Adding or removing features

The first step of the wokflow is to download the raw feature data. If you do not want to change any features and stick with the 26 fetaures (list can be found in XX), you do need to take any actions. If you want to add or remove features, you should first take look at the file *config/featuresConfig38.json*. Here, each feature has a entry with following specifications:

```
 "feature name": {
            "url": link to the file that should be downloaded,
            "file" : name to give to the downloaded file,
            "files" : if a feature consists of many files that have to be downloaded, than you have to write the names of the files in a list; if only one file will be downloaded, write "all",
            "type": type of the file  before it is converted into VCF by AttibutedDB (see Section X),
            "method": methods for processing the files by AttibutedDB (see Section X)
            "description": short description of the feature
        }
```
After you have prepared this, add the name of the feature in the list at the top of the JSON file. Next, you have to add a snakemake rule for downloading the fetaure. For that, create a new text file *featureName.snakerule* in the folder rules/features and add a rule that uses the new entry in the *config/featuresConfig38.json* for donwloading the feature and naming the file(s). You can use one of the existing rules as a scaffold. Here, you can preprocess the data if needed or use a separate rule for that (use the same snakerule file for that). At the end, you need to output a comma-separated file, that will be used in the next step.

To remove a feature from computation of ReMM, you need only to remove its name from the top list in *featuresConfig38.json*. 

### Converting into VCF
After you have downloaded the raw feature and 


The Snakemake workflow consists of XX different rules, many of which are excecuted multiple times. 
In the following, we briefly discuss what most of the rules thus and what changes you need to do, to adapt the wokflow to your data. 


Each feature is downloaded and preprocessed by a separate Snakemake rule (*rules/features*). This modularization is needed due to different data sources and due to the non-identical preprocessing steps. The download  is managed by a config file (config/featuresConfig38.json) that contains the download links and the nessasary meta information for the processing steps.



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
