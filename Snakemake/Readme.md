To DoS:
- Complete the ReMM environment file, create it during the workflow
- Add a link for downloading the negative set (in 4.2 Benign variants)

## Project workflow
### Folder structure

The folder Snakemake contains several subfolders:

- config: configuration files for downloading and processing of annotations as well as configuration files for running parSMURF
- rules: Snakemake rules separated into subfolders according to the workflow steps
- env: a YAML file for creating the main working environment 'ReMM' as well as files for some additional environments required by certain rules
- scripts: (external) scripts used by Snakemake rules
- utils: diverse reference files and scaffolds used during the workflow
- input: files that serve as input for  rules (raw feature files, variants)
- output: computed files - output of rules (VCF feature files, training sets, predictions)

Some directories contain README files with more detailed information on the content of folders. A DAG (directed acyclical graph) of the workflow (schematic representation) starting from the step *4. Variants* can be found in WorkflowDAG.jpg. The DAG for steps 1-3 is not shown for the sake of clarity because they contain more than one hundred rule executions.

### 1. Collecting annotations

Each annotation is downloaded and preprocessed by a separate Snakemake rule (*rules/features*). This modularization is needed due to different data sources and due to the non-identical preprocessing steps. The download is managed by a **configuration file** (*config/featuresConfig38.json*) that contains the download links and the necessary meta information for the processing steps. Each feature has the following entry:


```
 "feature name": {
            "url": link to the file that should be downloaded,
            "file" : name to give to the downloaded file,
            "files" : if feature data is spread across files that have to be downloaded, than you have to write the names of the files in a list; if only one file will be downloaded, write "all",
            "type": type of the file  before it is converted into VCF by AttibuteDB (see below),
            "method": methods for processing the files by AttibuteDB (see below),
            "description": short description of the feature
        }
```
If you want to add a feature for the calculation of the ReMM score, you have to define it as descibed above and add the name of the feature in the list *feature_set["hg38"]* at the top of the JSON configuration file. Next, you have to add a snakemake rule for downloading the feature. For that, create a new text file *featureName.snakerule* in the folder *rules/features* and add a rule that uses the new entry in the *config/featuresConfig38.json* for downloading the feature and naming the file(s). You can use one of the existing rules as a scaffold. In the rule, you can preprocess the data if needed or use a separate snakemake rule for that (use the same snakerule file). At the end, you need to output a comma-separated file in the input/featureName folder, that will be used in the next step. All following steps do not need to be modified in order to add a new feature.

To remove a feature from computation of ReMM, you need only to remove its name from the top list *feature_set["hg38"]* in *featuresConfig38.json*. Features not defined in the list will be not further processed even if the feature has an entry in the configuration file. 

### 2. Converting into VCF
After raw feature files are downloaded and processed, they have to be converted into VCF format. This is done by the program **AttributeDB** that needs a **property file** for processing features. The file contains following information:

```
name = name of the feature 
file = link to the file(s)
type = type of the file (the program accepts bed and wig files) 
method = upload or max. If a feature has multiple values for one position defined in multiple files (e.g. histone modification in different cell types), the program will retrieve the maximum value of the feature from the files at that position when 'max' is specified. Otherwise, use 'upload' for simple upload of values from one file.
column = column number where the feature values are defined, usually 4, but it can diverge.
description = short description of the feature

```
The property file is created by the rule *createPropertyFile* in the Snakefile. It processes the content of the feature configuration file and creates a property file  *input/hg38/PropertyFiles/featureName.properties*. The VCF files are created by the rule *createSingleFeatureVCF* and indexed by *indexSingleFeatreVCF* from Snakefile.

### 3. Creating feature set
The individual feature files are merged into one VCF after all features have been downloaded, processed and converted into VCF. The rule *mergeSingleFeatureVCF* in Snakemake file uses bcftools for that. It outputs a file containing all 26 features and creates an index file for it. Both are used in the next step for annotating variants.

### 4. Variants
#### 4.1 Pathogenic variants
The set of 456 positive variants is available in the folder *utils* under *SNVs.all.20160109.vcf.gz*. These are Indels and SNP curated for ReMM GRCh37 study, so in the first step the positions have to be lifted over to GRCh38 coordinates, which is done by the rule *liftOverPositive* (in *rules/process*) using the UCSC liftOver tool. After that, the SNPs are filtered and saved into a VCF file by the script *filter.py*  applied by the rule *getPositiveSNPs*.

#### 4.2 Benign variants
A set of 14 million variants can be downloaded here (LINK??). The file contains coding and non-coding variants, so we, first, need to annotate the positions and then filter for non-coding variants. Annotation is done by Jannovar (in rule *annotateJannovarNegative*) that needs a refseq library saved in *utils/data/hg38_refseq.ser*. The annotated file in processed by the rule  *jannovarFilter* to filter out the coding variants. The resulting file *input/variants/hg38/SNVs.hg38.negative.refseq.filtered.vcf.gz* contains around 13.9 million non-coding variants. 
If you want to use a different set, you need to keep the structure of the new file identical to original one as it is important for further steps.

### 5. Variant annotation
Next step of the workflow is to annotate the positive and negative variants with features. This is done in the rule *annotateFeatures* (Snakefile) that hands over the VCF files of positive and negative variants together with the feature VCF file to the **AttibuteDB** that annotates the files. This is very time consuming and needs many hours for completion. The annotated variants are saved in *input/variants/*.

### 6. Training set 
The positive and negative sets are first combined into one file in *combineInputData* and then processed by the rule *createParsmurfInput* to create the input for parSMURF. The input consists of three parts: features, labels and folds. Folds are needed for a special cross-validation technique that handles the locally correlated structure of variants by cross-validating on folds that contain no correlated data. The folds are created according to cytogenic bands in the script *createParsmurfInput.py*.

### 7. Training and Cross-validation
Training  and validation are performed by **parSMURF**. It applies the hyperSMURF method for training with unbalanced data (see paper *Imbalance-Aware Machine Learning for Predicting Rare and Common Disease-Associated Non-Coding Variants*). The executable of parSMURF is saved in the *bin* folder. ParSMURF runs in three different modi: training, cross-validation and prediction (train, cv, predcit in the cofig file). Which modus is to be used as well as other details are defined in the **parSMURF configuration file**. A scaffold of it can be found in *utils*. The rule *generateParsmurfConfig* creates a configuration file basing on the scaffold and the name of the output file that is handed over to Snakemake. The name contains the modus and a seed; if no seed is defined, the default see *1* is used. For example, calling snakemake with the file name *output/predictions/hg38/SNVs.hg38.cv.predictions.txt* will make parSMURF to run in cross-validation modus with the default seed. If you want to change the seed, you can define the file name as  *output/predictions/hg38/SNVs.hg38.cv.predictions.txt_seed*

ParSMURF outputs a tab-delimited file with containing two columns. The first column is the probability *p* of a variant to be pathogenic at the position, the second is the probability *1 - p*, so that each row sums up to *1*. 

### 8. Genome-wide ReMM score
Not done yet

---------------------


Absolute paths:
/fast/groups/ag_kircher/CADD/dependencies/annotations/gerp/gerp2_elements_hg38_MAM.bg.gz
/fast/groups/ag_kircher/CADD/cadd_v1.3/training_data/GRCh38/humanDerived/annotated/SNVs.vcf.gz
/fast/work/groups/ag_kircher/ReMM/ReMM/data/variants/RegulatoryMendelianMutations/GRCh37/SNVs.all.20160109.vcf.gz
