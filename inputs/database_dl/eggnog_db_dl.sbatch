#!/usr/bin/env bash
#SBATCH --job-name=eggnog_db-dl
#SBATCH --time 12:00:00                  
#SBATCH --partition [YOUR_PARTITION_NAME]              
#SBATCH --cpus-per-task 1          
#SBATCH --mem-per-cpu=50000     
#SBATCH -o eggnog_db-dl.%j.out
#SBATCH -e eggnog_db-dl.%j.err

date

# install eggnog-mapper (comment out installation if already done)
mamba create --name eggnog --channel bioconda eggnog-mapper

# activate eggnog environment 
mamba activate eggnog

# download database for eggnog
download_eggnog_data.py 

echo Done!
date