import json

with open(snakemake.input.scaffold) as json_file:
    file = json.load(json_file)
    
    file['name'] = snakemake.params.name
    file['data']['dataFile'] = snakemake.input.data
    file['data']['labelFile'] = snakemake.input.labels
    file['data']['foldFile'] = snakemake.input.folds
    file['data']['outFile'] = snakemake.params.predictions
    file['exec']['seed'] = int(list(snakemake.params.seed)[0])
    file['exec']['mode'] =snakemake.params.mode
    file['data']['forestDir'] = "models/"+snakemake.params.genomeBuild
    
        
with open(snakemake.output.config, 'w') as outfile:
    json.dump(file, outfile)

