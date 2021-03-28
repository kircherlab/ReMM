## Project workflow
### Folder structure

The folder Snakemake contains several subfolders:

- config: configuration files for downloading and processing of annotations as well as configuration files for running parSMURF
- rules: snakemake rules separated into subfolders according to the workflow steps
- env: a yaml file for creating the main working environment 'ReMM' as well as files for some additional environments required by certain rules
- scripts: (external) scripts used by snakemake rules
- utils: diverse reference files and scaffolds used during the workflow
- input: files that serve as input for  rules (raw fetaure files, variants)
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
            "method": methods for processing the files by AttibutedDB (see Section X),
            "description": short description of the feature
        }
```
After you have prepared this, add the name of the feature in the list *feature_set["hg38"]* at the top of the JSON file. Next, you have to add a snakemake rule for downloading the fetaure. For that, create a new text file *featureName.snakerule* in the folder *rules/features* and add a rule that uses the new entry in the *config/featuresConfig38.json* for donwloading the feature and naming the file(s). You can use one of the existing rules as a scaffold. Here, you can preprocess the data if needed or use a separate snakemake rule for that (use the same snakerule file). At the end, you need to output a comma-separated file in the input/featureName folder, that will be used in the next step.

To remove a feature from computation of ReMM, you need only to remove its name from the top list *feature_set["hg38"]* in *featuresConfig38.json*. Features not defined in the list, will be not further processed even if the feature is defined as shown above. 

### Converting into VCF
After raw feature files are downloaded and processed, they have to be converted into VCF format. This is done by the progamm AttributeDB that needs a proporty file for processing features. The file contains following information:

```
name = name of the feature 
file = link to the file(s)
type = type of the file (the programm accepts bed and wig files) 
method = upload or max. If a feature is defined in more than ones for one and the same position across multipple files (e.g. histone modification in different cell types), the programm will retrieve the maximum value of the feature from the files at that position when 'max' is specified. Otherwiese use 'upload' for simple upload of values from one file.
description = short description of the feature
column = column number where the feature values are defined, usually 4 but it can diverge.

```
The property file is created by the rule *createPropertyFile* in the Snakefile. It processes the content of the feature configuration file and creates a property file  *input/hg38/PropertyFiles/featureName.properties*. The VCf files are created by the rule *createSingleFeatureVCF* and indexed by *indexSingleFeatreVCF* from Snakefile.

### Creating feature set
The individual feature files are merged into one VCF after all features have been downloaded, processed and converted into VCF. The rule *mergeSingleFeatureVCF* in Snakemake file uses bcftools for that. It ouputs a file of 53.9Gb (for existing 26 fetaures) and creates and index file. Both are used in the next step for annotating variants.

### Variants
#### Pathogenic variants
The set of 456 positive variants is available in the folder *utils* under *SNVs.all.20160109.vcf.gz*. These are Indels and SNP curated for ReMM GRCh37 study, so in the first step the positition have to be lifted over to GRCh38 coordniates, which is done by the rule *liftOverPositive* (in *rules/process*) using UCSC liftOver tool. After that, the SNPS are filtered and saved into a VCF file by the Script *filter.py*  applied by the rule *getPositiveSNPs*.

#### Benign variants
A set of ca. 14 million variants can be downloaded here (LINK??). The file contains coding and non-coding variants, so we first need to annotate the positions and then filter for non-coding variants. Annotation is done by Jannovar (in rule *annotateJannovarNegative*) that needs a refseq library saved in *utils/data/hg38_refseq.ser*. The annotated file in processed by the rule  *jannovarFilter* to filter out the coding variants. The resulting file *input/variants/hg38/SNVs.hg38.negative.refseq.filtered.vcf.gz* contains around 13.8 million non-coding variants.

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
