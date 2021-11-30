# TODO Lusi


### Rule all zum testen in meinem Ordner hg38 = test_hg38 (project=genomeBuild, bei dir anders)
### habe als Beispiel für {N} drin gelassen
#rule all:
 #   input:
  #      expand("results/validation/identical.region{N}.{genomeBuild}.csv.gz", genomeBuild = ['hg38'],       N =[0,1,2] )
       #  "results/validation/random.hg38.csv.gz"     

    
rule getRandomPositions:
     # ich würde eventuell diesen pfad auch aus der config entnehmen, weil wenn der Workflow schief läuft, wird Snakemake versuchen diese Datei auch zu löschen
    input:
        "results/features/feature_sets/{project}.vcf.gz",
        "results/features/feature_sets/{project}.vcf.gz.tbi"
    output:
        "results/validation/random.{project}.csv.gz" # werden 3 temp Dateien erstellt, die im script direkt gelöscht werden          
    params:
        genome = lambda wc: config['global_files']['genome_builds'][wc.project]['genome'],
        N = 120000  
    shell:
        """scripts/getRandomPositions.sh {input} {params.N} {params.genome} {output}"""
        
        
rule getIdenticalRegions:
    input:
        "results/features/feature_sets/{project}.vcf.gz",
        "results/features/feature_sets/{project}.vcf.gz.tbi"
    output:
        "results/validation/identical.region{N}.{project}.csv.gz",
        temp( "results/validation/identical.region{N}.{project}.csv.gz.temp.gz"),
        temp( "results/validation/identical.region{N}.{project}.csv.gz.temp.gz.tbi")# N ist Listeneintrag aus der Config (0,1,2)
    params:
        region = lambda wc: config['feature_sets'][wc.project]['validation'][int(wc.N)]
    shell:
         """scripts/getIdenticalRegions.sh {params.region} {input} {output}"""
