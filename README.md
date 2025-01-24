# GenomeAssembly-Snakemake
Genome Assembly Pipeline for Illumina short read and Oxford Nanopore long read. 

## Setup
### 1. Install Snakemake and Conda/Mamba  
Install Snakemake and Conda/Mamba following the instructions at this [link](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#:~:text=for%20installing%20Snakemake.-,Installation%20via%20Conda/Mamba,-This%20is%20the). 

### 2. Download Databases  
#### 1. CheckM2 Database  
Install the checkM2 Database by running script in `inputs/database_dl/` directory:   
    ```
    sbatch checkM2_db_dl.sbatch  
    ```  
Note: this script creates a new conda environment name "genome-assembly-checkm2" and downloads database into the same directory. Change path in script or run script in another folder if preferred. 

#### 2. GTDB-tk Database  
Install the GTDB Database by running script in `inputs/database_dl/` directory:   
    ```
    sbatch GTDB_dl.sbatch  
    ```  
Note: this script downloads database into the same directory. Change path in script or run script in another folder if preferred.  

### 3.Set up Snakemake Pipeline
#### 1. Experimental Configurations
Edit **config.yaml** file in the `inputs/` Directory:
  - Edit relative path to sample table in line 4
  - Edit "seq data type" in line 5. Options: "illumina short read"  or "nanopore long read"
  - Edit paths to reference databases in lines 9 & 10
  - Edit paths to scratch (intermediate) directory and results directory in lines 12 & 13. 

Create **samples.tsv** file in the `inputs/` Directory: 
  - Create samples.tsv file for your samples with the following required columns (and any other columns for your samples): 
  - Required columns for Illumina samples: ["forward read", "reverse read", "sample"]
    - "forward read" and "reverse read": paths to forward and reverse read .fastq file
    - "sample": unique sample identifier
  - Required columns for Nanopore samples: ["read path", "sample"]
    - "read path": path to .fastq file
    - "sample": unique sample identifier

#### 2. Resource Specifications 

## Running Pipeline 

## Troubleshooting 

## Interpreting Results

