The Regulatory Mendelian Mutation (ReMM) score was created for relevance prediction of non-coding variations (SNVs and small InDels) in the human genome (GRCh37) in terms of Mendelian diseases. This project updates the ReMM score for the genome build GRCh38. 

## Add the snakemake env file 

The entire workflow is managed by Snakemake - a workow management system used to create reproducible and scalable data analyses. The workflow consist of the following parts:
    - Download of data
    - Data processing and cleansing 
    - Model training and validation 
    
These steps are performed by number of scripts that are managed in the main workflow file called Snakemake. 
    
Download of data.
Each feature is downloaded and preprocessed by a separate Snakemake rule (rules/features). This modularization is needed due to different data sources and due to the non-identical preprocessing steps. The download  is managed by a config file (config/featuresConfig38.json) that contains the download links and the nessasary meta information for the processing steps.
