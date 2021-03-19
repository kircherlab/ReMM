## MA-ReMM

The Regulatory Mendelian Mutation (ReMM) score was created for relevance prediction of non-coding variations (SNVs and small InDels) in the human genome (GRCh37) in terms of Mendelian diseases. This project updates the ReMM score for the genome build GRCh38. 

## Pre-requirements
### Conda
We use Conda as software/dependency managment tool. Conda installation guidlines can be found here:

https://conda.io/projects/conda/en/latest/user-guide/install/index.html

### Snakemake
The workflow is managed by Snakemake - a workflow management system used to create reproducible and scalable data analyses. To install Snakemake as well as all other required packages, you need to create a working environment  according to the description in the file env/ReMM.yaml. For that, first

Clone the repository 
```
git clone https://github.com/nazaretl/ReMM-GRCh38
cd ReMM-GRCh38

```

Create  working environment

```
conda env create -n ReMM --file env/ReMM.yaml

```




