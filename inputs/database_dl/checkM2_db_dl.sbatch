#!/usr/bin/env bash
#SBATCH --job-name=checkM2_db-dl
#SBATCH --time 12:00:00                  
#SBATCH --partition [YOUR_PARTITION_NAME]              
#SBATCH --cpus-per-task 1          
#SBATCH --mem-per-cpu=50000     
#SBATCH -o checkM2_db-dl.%j.out
#SBATCH -e checkM2_db-dl.%j.err

date

# install checkm2 (comment out installation if already done)
mamba create --name genome-assembly-checkm2 --channel bioconda checkm2 numpy=1.23.4

# activate checkm2 environment 
mamba activate checkm2

# download database for checkm2
# NOTE: edit '.' to full path of preference 
checkm2 database --download --path .

echo Done!
date