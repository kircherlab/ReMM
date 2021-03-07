#!/bin/bash
# Set a name for the job (-J or --job-name).
#SBATCH --job-name=snakemake0

# Set the file to write the stdout and stderr to (if -e is not set; -o or --output).
#SBATCH --output=/fast/groups/ag_kircher/ReMM/MA_Lusi/Snakemake/logs/%x-%j.log

# Set the number of cores (-n or --ntasks).
#SBATCH --ntasks=1

# Force allocation of the two cores on ONE node.
#SBATCH --nodes=1

# Set the memory per CPU. Units can be given in T|G|M|K.
#SBATCH --mem-per-cpu=50G

# Set the partition to be used (-p or --partition).
#SBATCH --partition=long

# Set the expected running time of your job (-t or --time).
# Formats are MM:SS, HH:MM:SS, Days-HH, Days-HH:MM, Days-HH:MM:SS
#SBATCH --time=200:30:00


awk 'BEGIN{ OFS="\t"} NR==FNR{ a[$1 FS $2];next} ($1 FS $2) in a {print "chr"$1, $2, $3} ' <(zcat output/predictions/lifted/new/hg38.predictions.lifted.sorted.txt.gz) <(zcat input/variants/hg19/ReMM.v0.3.1.tsv.gz) |gzip -c  > output/predictions/lifted/new/SNVs.lifted.remm.txt.gz 

