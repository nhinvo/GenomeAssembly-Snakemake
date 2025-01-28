"""
Purpose: to parse and combine: binning stats, bin quality, bin classification. 
"""
import pandas as pd 
from pathlib import Path 

### import paths to data tables ### 
bindir_list = snakemake.input['bin_dirs']  # list of directories containin bins 
quality_list = snakemake.input['qualities']  # list of checkm quality outputs
classification_list = snakemake.input['classifications']  # list of gtdb-tk classification outputs 
output = snakemake.output[0]  # path to output final .tsv file 

def process_bins(bindir_list):
    """
    Returns a df of information on all bins from metabat2. 
    """
    data = []

    # cycle through every bins 
    for bindir_fpath in bindir_list:
        bindir_fpath = Path(bindir_fpath) 
        sample = bindir_fpath.name

        for bin_fpath in Path(bindir_fpath).glob('*fa'):
            bin_id = bin_fpath.stem.split('.')[-1]  # e.g. '1' from full/file/path/Paragon14_34.1.fa

            # count number of contigs and length 
            contig_count = 0
            bin_len = 0
            with open(bin_fpath, 'r') as file:
                for line in file:
                    if line.startswith('>'):
                        contig_count += 1
                    else:
                        line_len = len(line)
                        bin_len += line_len

            bin_data = {
                'sample': sample, 
                'bin_id': bin_id, 
                'bin_contig_num': contig_count, 
                'bin_length': bin_len, 
            }

            data.append(bin_data)

    df = pd.DataFrame(data)

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
    bin_df = process_bins(bindir_list)
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

    # save data 
    df.to_csv(output, sep='\t', index=False)

main()