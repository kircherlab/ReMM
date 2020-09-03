#!/bin/bash

# Set a name for the job (-J or --job-name).
#SBATCH --job-name=convertToBed

# Set the file to write the stdout and stderr to (if -e is not set; -o or --output).
#SBATCH --output=/fast/groups/ag_kircher/ReMM/MA_Lusi/Snakemake/logs/%x-%j.log

# Set the number of cores (-n or --ntasks).
#SBATCH --ntasks=1

# Force allocation of the two cores on ONE node.
#SBATCH --nodes=1

# Set the memory per CPU. Units can be given in T|G|M|K.
#SBATCH --mem-per-cpu=20G

# Set the partition to be used (-p or --partition).
#SBATCH --partition=medium

# Set the expected running time of your job (-t or --time).
# Formats are MM:SS, HH:MM:SS, Days-HH, Days-HH:MM, Days-HH:MM:SS
#SBATCH --time=10:30:00
conda activate snakemake-tutorial;

export LC_ALL=C;
java -Xmx5g -jar bin/attributedb-cli-0.0.1-jar-with-dependencies.jar \
vcf -p input/features/hg38/PropertyFiles/mamPhyloXway.properties --output output/features/single/hg38/mamPhyloXway.vcf.gz;
