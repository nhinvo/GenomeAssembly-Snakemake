"""
Purpose: to plot bin completeness vs. contamination. 
"""
import pandas as pd 
import numpy as np
from pathlib import Path 
import seaborn as sns 
import matplotlib.pyplot as plt 

def process_checkm(checkm_fpaths):
    """
    """
    df = pd.concat([
        pd.read_table(fpath) for fpath in checkm_fpaths
    ])

    # determine whether bin passes quality threshold 
    completeness_threshold = snakemake.config['CheckM2']['completeness threshold']
    contamination_threshold = snakemake.config['CheckM2']['contamination threshold']

    df['Status'] = np.where(
        (df['Completeness'] > completeness_threshold) & 
        (df['Contamination'] < contamination_threshold), 
        "Passed CheckM Thresholds", "Failed CheckM Thresholds"
    )

    return df


def plot(df):
    """
    Plot Completeness vs. Contamination for all bins. 
    """
    # dynamic Contamination y-lim 
    y_range = df['Contamination'].max() - df['Contamination'].min()
    ylim = (
        df['Contamination'].min()-(y_range*0.1),
        df['Contamination'].max()+(y_range*0.1),
    )

    # custom colors
    palette = {
        'Passed CheckM Thresholds': '#70cc73', 
        'Failed CheckM Thresholds': '#d46161', 
    }

    scatter_fig = sns.jointplot(
        data=df, 
        x="Completeness", 
        y="Contamination", 
        hue="Status", 
        s=50, edgecolor="black", linewidth=0.2,
        xlim=(-3,103), ylim=ylim, marginal_ticks=True, 
        palette=palette, 
    )

    scatter_fig.fig.suptitle("CheckM2 Completeness vs. Contamination\nFor All Metabat2 Bins", fontsize=13)
    scatter_fig.fig.subplots_adjust(top=0.90)  

    plt.savefig(snakemake.output['plot_out'])
    plt.close()


def main():
    # import all quality output tables & process 
    checkm_fpaths = snakemake.input['quals']
    df = process_checkm(checkm_fpaths)

    # plot completeness vs. contamination 
    plot(df)

main()