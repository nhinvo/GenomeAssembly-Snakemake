import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def process_final_tsv(fpath):
    """
    Process final unfiltered assembly df for plotting. 
    """
    df = pd.read_table(fpath)[['sample', 'genus', 'Completeness', 'Contamination']]

    # sort by completeness but keep sample order
    df = df.sort_values(by=['sample', 'Completeness'], ascending=[True, False])  
    df = df.drop_duplicates(subset=['sample', 'genus'], keep='first') 

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
    Plot classification of each sample. 
    """
    # custom color pallete 
    palette = {
        'Passed CheckM Thresholds': '#70cc73', 
        'Failed CheckM Thresholds': '#d46161', 
    }

    # order the genus by value counts 
    x_axis_value_order = df['genus'].value_counts().index

    # obtain plot size
    # height = number of samples + space for x-axis labels
    height = df['sample'].nunique() * 0.2 + 2
    # width = number of unique classification + space for legend
    width = df['genus'].nunique() * 0.5 + 2  

    # plotting
    plt.figure(figsize=(width, height))

    # bubble / categorical scatter of classification 
    sns.stripplot(
        data=df, 
        x='genus', y='sample', hue='Status', 
        size=8, palette=palette, edgecolor='black', linewidth=1,
        order=x_axis_value_order, jitter=False, 
    )

    # grid for viewing ease 
    plt.grid(True, linestyle='--', linewidth=0.5, color='grey', alpha=0.5)

    # ticks labels and legends
    plt.xticks(rotation=90)
    plt.xlabel("Taxonomic Classification (Genus)")
    plt.ylabel("Sample Name")
    plt.title("Classification of Metabat2 Bins")
    plt.legend(title="CheckM2 Threshold", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    # save plot 
    plt.savefig(snakemake.output['plot_out'])
    plt.close()



def main():
    # import final table with unfiltered classification and quality 
    final_tsv = snakemake.input['table']
    df = process_final_tsv(final_tsv)

    plot(df)


main()
