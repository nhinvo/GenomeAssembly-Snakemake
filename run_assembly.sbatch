#!/usr/bin/env bash
#SBATCH --job-name=SMK-assembly
#SBATCH --time 5-0                        
#SBATCH --partition sched_mit_chisholm              
#SBATCH --cpus-per-task 1                                
#SBATCH -o logs/main.%j.out
#SBATCH -e logs/main.%j.err

eval "$(conda shell.bash hook)"
conda activate snakemake

snakemake --profile profile 
