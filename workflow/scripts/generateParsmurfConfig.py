import json

with open(snakemake.input.scaffold) as json_file:
    file = json.load(json_file)
    
    file['name'] = snakemake.params.name
    file['data']['dataFile'] = snakemake.input.data
    if "labels" in snakemake.input.keys():
        file['data']['labelFile'] = snakemake.input.labels
    if "folds" in snakemake.input.keys():
        file['data']['foldFile'] = snakemake.input.folds
    file['data']['outFile'] = snakemake.params.predictions
    file['data']['forestDir'] = snakemake.params.models
    if "ensThrd" in snakemake.params.keys():
        file['exec']['ensThrd'] = snakemake.params.ensThrd
    file['exec']['seed'] = int(snakemake.params.seed)
    file['exec']['mode'] =snakemake.params.mode
    
    
        
with open(snakemake.output.config, 'w') as outfile:
    json.dump(file, outfile)

