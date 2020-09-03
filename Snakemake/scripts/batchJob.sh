#!/bin/bash

# Set a name for the job (-J or --job-name).
#SBATCH --job-name=convertToBed

# Set the file to write the stdout and stderr to (if -e is not set; -o or --output).
#SBATCH --output=/fast/groups/ag_kircher/ReMM/MA_Lusi/Snakemake/logs/%x-%j.log

# Set the number of cores (-n or --ntasks).
#SBATCH --ntasks=2

# Force allocation of the two cores on ONE node.
#SBATCH --nodes=1

# Set the memory per CPU. Units can be given in T|G|M|K.
#SBATCH --mem-per-cpu=10G

# Set the partition to be used (-p or --partition).
#SBATCH --partition=medium

# Set the expected running time of your job (-t or --time).
# Formats are MM:SS, HH:MM:SS, Days-HH, Days-HH:MM, Days-HH:MM:SS
#SBATCH --time=5:30:00


conda activate snakemake-tutorial;
conda activate snakemake-tutorial;

zcat  /fast/groups/ag_kircher/ReMM/MA_Lusi/Snakemake/input/features/hg38/priPhastConsXway/hg38.phastCons17way.wig.gz | wig2bed | bgzip  >  /fast/groups/ag_kircher/ReMM/MA_Lusi/Snakemake/input/features/hg38/priPhastConsXway/_hg38.phastCons17way.bed.gz;

zcat /fast/groups/ag_kircher/ReMM/MA_Lusi/Snakemake/input/features/hg38/priPhastConsXway/_hg38.phastCons17way.bed.gz; | cut -f 1-3,5 | bgzip -c  > /fast/groups/ag_kircher/ReMM/MA_Lusi/Snakemake/input/features/hg38/priPhastConsXway/__hg38.phastCons17way.bed.gz;