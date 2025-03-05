"""
Purpose: to parse and combine: binning stats, bin quality, bin classification. 
"""
import pandas as pd 
from pathlib import Path 

### import paths to data tables ### 
# inputs
bindir_list = snakemake.input['bin_dirs']  # list of directories containin bins 
bin_depths = snakemake.input['bin_depths']   # list of bin depth paths 
quality_list = snakemake.input['qualities']  # list of checkm quality outputs
classification_list = snakemake.input['classifications']  # list of gtdb-tk classification outputs 
# outputs
output = snakemake.output['final_tsv']  # path to output final .tsv file 
output_filtered = snakemake.output['final_tsv_filtered']  # path to output final .tsv file with checkM filtering

def process_sample_bin(depth_fpath):
    """
    Return a df with cols: [sample, bin_id, contig_name] for 
    all bins in sample. 
    """
    depth_fpath = Path(depth_fpath)
    sample_name = depth_fpath.stem.replace('_depth', '')
    bin_dir = f"{depth_fpath.parent.parent}/bins/{sample_name}"

    # obtain contigs (headers) in each bin 
    df_data = {
        'sample': [], 'bin_id': [], 'contigName': [], 
    }

    for bin_fpath in Path(bin_dir).glob('*.fa'):
        bin_id = bin_fpath.stem.split('.')[-1]

        with open(bin_fpath, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    contig_name = line.strip().replace('>', '')
                    df_data['sample'].append(sample_name)
                    df_data['bin_id'].append(bin_id)
                    df_data['contigName'].append(contig_name)

    df = pd.DataFrame(df_data)

    return df


def process_bins(bin_depths):
    """
    Returns a df of information on all bins from metabat2. 
    """
    dfs = []  # list to store all sample bin contig data 

    # 1. For each bin, obtain depth and contig data 
    for depth_fpath in bin_depths: 
        # import average depth for each contig 
        depth_df = pd.read_table(depth_fpath)[['contigName', 'totalAvgDepth', 'contigLen']]
        print(depth_df.head(5))

        # obtain all contigs in each bin 
        bin_df = process_sample_bin(depth_fpath)
        print(bin_df.head(5))

        # merge to obtain depth and length for each contig (that were binned)
        df = pd.merge(bin_df, depth_df, on=['contigName'], how='left')
        print(df.head(5))

        dfs.append(df)

    df = pd.concat(dfs)

    # 2. obtain average depth for each bin and total length
    groups = df.groupby(['sample', 'bin_id'])

    dfs = []

    for index, gdf in groups:
        # skip average and sum if bin has only 1 contig 
        if len(gdf) == 1:
            gdf = gdf.drop(columns=['contigName'])
            gdf['bin_contig_num'] = 1
            dfs.append(gdf)
            continue 

        gdf = pd.DataFrame({
            'sample': [index[0]], 
            'bin_id': [index[1]], 
            'totalAvgDepth': [gdf['totalAvgDepth'].mean()],
            'contigLen': [gdf['contigLen'].sum()],
            'bin_contig_num': len(gdf), 
        })

        dfs.append(gdf)

    df = pd.concat(dfs)

    df = df.rename(columns={'totalAvgDepth': 'bin_avg_depth', 'contigLen':'bin_length'})

    return df 

def process_classification(classification_list):
    """
    Returns df of aggregate data from GTDB-tk classification. 
    """
    df = pd.concat(pd.read_table(fpath) for fpath in classification_list)

    df['sample'] = df['user_genome'].str.split('.').str[0]
    df['bin_id'] = df['user_genome'].str.split('.').str[-1]

    # split up classification 
    taxonomic_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    for counter, rank in enumerate(taxonomic_ranks):
        df['classification'] = df['classification'].str.replace('Unclassified Bacteria', 'd__Unclassified Bacteria')
        df[rank] = df['classification'].str.split(pat=';', expand=True)[counter].str[3:]

    return df 

def process_quality(quality_list):
    """
    Returns df of aggregated data from checkM2 quality assessment. 
    """
    df = pd.concat(pd.read_table(fpath) for fpath in quality_list)

    df['sample'] = df['Name'].str.split('.').str[0]
    df['bin_id'] = df['Name'].str.split('.').str[1]

    return df


def main():
    # obtain binning, quailty, and classification data for all samples
    bin_df = process_bins(bin_depths)
    quality_df = process_quality(quality_list)
    classification_df = process_classification(classification_list)

    # comebine all data
    df = pd.merge(bin_df, classification_df, on=['sample', 'bin_id'], how='outer')
    df = pd.merge(df, quality_df, on=['sample', 'bin_id'], how='outer')

    # col filtering 
    cols = [
        'sample',
        'bin_id',
        'bin_contig_num',
        'bin_length',
        'bin_avg_depth', 
        'Completeness', 
        'Contamination', 
        'kingdom',
        'phylum',
        'class',
        'order',
        'family',
        'genus',
        'species',
        'closest_genome_reference', 
        'closest_genome_ani',
        'warnings',
    ]

    df = df[cols]

    # sort df 
    df = df.sort_values(by=['sample', 'Completeness'], ascending=False)

    # save unfiltered data 
    df.to_csv(output, sep='\t', index=False)

    # filter df by checkm results
    df = df[
        (df['Completeness'] > snakemake.config['CheckM2']['completeness threshold']) &
        (df['Contamination'] < snakemake.config['CheckM2']['contamination threshold'])
    ]

    # save filtered data 
    df.to_csv(output_filtered, sep='\t', index=False)

main()