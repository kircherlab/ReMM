
import json

with open(snakemake.input.scaffold) as json_file:
    file = json.load(json_file)
    
    file['name'] = snakemake.params.name
    file['data']['dataFile'] = snakemake.input.features
    file['data']['labelFile'] = snakemake.input.labels
    file['data']['foldFile'] = snakemake.input.folds
    file['data']['outFile'] = 'output/' + snakemake.params.name + '/' + snakemake.params.name + '_predictions.txt'

with open(snakemake.output.config, 'w') as outfile:
    json.dump(file, outfile)

