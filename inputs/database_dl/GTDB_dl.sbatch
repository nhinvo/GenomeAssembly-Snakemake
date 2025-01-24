#!/usr/bin/env bash
#SBATCH --job-name=GTDB-dl
#SBATCH --time 12:00:00                  
#SBATCH --partition [YOUR_PARTITION_NAME]              
#SBATCH --cpus-per-task 1          
#SBATCH --mem-per-cpu=50000     
#SBATCH -o GTDB-dl.%j.out
#SBATCH -e GTDB-dl.%j.err

date

# NOTE: if needed, change directory into desired folder, or run script in that folder
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
tar xvzf gtdbtk_data.tar.gz

echo done
date