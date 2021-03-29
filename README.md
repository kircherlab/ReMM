## ReMM for GRCh38

The Regulatory Mendelian Mutation (ReMM) score was created for relevance prediction of non-coding variations (SNVs and small InDels) in the human genome (GRCh37) in terms of Mendelian diseases. This project updates the ReMM score for the genome build GRCh38.

## Pre-requirements

### Conda
We use Conda as software and dependency management tool. Conda installation guidelines can be found here:

https://conda.io/projects/conda/en/latest/user-guide/install/index.html

### Additional programs
These programs are used during the workflow. They usually need to be compiled, however, the repository already contains the executables or generated files.

- AttributeDB (https://github.com/visze/attributedb)
- Jannovar (https://doc-openbio.readthedocs.io/projects/jannovar/en/master/)
- parSMURF (https://github.com/AnacletoLAB/parSMURF)

### Snakemake

The workflow is managed by Snakemake - a workflow management system used to create reproducible and scalable data analyses. To install Snakemake as well as all other required packages, you need to create a working environment according to the description in the file env/ReMM.yaml. For that, first

Clone the repository
```
git clone https://github.com/kircherlab/ReMM-GRCh38
cd ReMM-GRCh38

```

Create a working environment and activate it

```
conda env create -n ReMM --file env/ReMM.yaml
conda activate ReMM
```

All paths are relative to the Snakemake file so you do not need to change any path variables. Additionally, Snakemake creates all missing directories, so no need to create any aditional folders either.

## Workflow

The workflow consists of four main parts:

- Download of data
- Data processing and cleansing
- Model training and validation
- Calculation of ReMM for the whole genome 

Snakemake folder contains a graph of the workflow and more detailed information on the most important steps.

To launch a snakemake workflow, you need to tell snakemake which file you want to generate. For example, to compute the cross-validation ReMM scores for GRCh38, you need to run:

```
snakemake output/predictions/hg38/SNVs.hg38.cv.predictions.txt
```

To execute any step separately (see Readme in the Snakemake folder for details on workflow steps), you need to look up the name of the desired output file in the scripts and call Snakemake with the exact name. Using a flag -n, you can initiate a 'dry run': Snakemake will check the consistency of all rules and files and show the number of steps. However, a clean dry run does not necessarily mean that no errors will occur during a normal run. ReMM score is not allele-specific so that you get only one score independent of the variant itself. The workflow from the download of data up to computing the scores may take several days or weeks depending on the computing power. On a XX machine with XX cores, the running time is XX hours.

