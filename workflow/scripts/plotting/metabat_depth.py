"""
Purpose: to plot bin completeness vs. contamination. 

Note: not useful data - not integrated into pipeline but will keep script. 
"""
import pandas as pd 
import numpy as np
from pathlib import Path 
import seaborn as sns 
import matplotlib.pyplot as plt 

def process_bins(bin_dirs):
    """
    """
    for dir_path in bin_dirs:
        for fpath in Path(dir_path).glob('*fa'):
            bin_name = fpath.stem
            print(bin_name)
            # with open(fpath, )
            print(fpath)
            with open(fpath, 'r') as file:
                for line in file:
                    if line.startswith('>'):
                        print(line)

            break
        break

def process_depth(depth_fpaths):
    """
    """
    df = pd.concat([
        pd.read_table(fpath, usecols=['contigLen', 'totalAvgDepth']) for fpath in depth_fpaths
    ], ignore_index=True)

    df = df[df['totalAvgDepth'] > 1]

    sample_size = 10000
    if len(df) > sample_size:
        df = df.sample(sample_size, random_state=42)

    df['contigLen'] = np.log(df['contigLen'] + 1e-6)
    df['totalAvgDepth'] = np.log(df['totalAvgDepth'] + 1e-6)
        
    return df


def plot(df):
    """
    Plot Completeness vs. Contamination for all bins. 
    """
    scatter_fig = sns.jointplot(
        data=df, 
        x="contigLen", 
        y="totalAvgDepth", 
        s=20, alpha=0.5, edgecolor="black", linewidth=0.2,
        marginal_ticks=True,  
    )

    scatter_fig.fig.suptitle("Metabat2 Contig Length vs. Depth\nFor All Assembled Contigs", fontsize=13)
    scatter_fig.fig.subplots_adjust(top=0.90)  

    plt.savefig("trash.png")
    plt.close()


def main():
    # import contig names in bins
    bin_dirs = [dir_path for dir_path in Path("../../../scratch/metabat_binning/bins").glob('*')]  # TEMP
    bin_df = process_bins(bin_dirs)

    return False

    # checkm_fpaths = snakemake.input[0]
    depth_fpaths = [fpath for fpath in Path("../../../scratch/metabat_binning/depth").glob('*txt')]  # TEMP

    df = process_depth(depth_fpaths)

    # plot completeness vs. contamination 
    plot(df)

main()