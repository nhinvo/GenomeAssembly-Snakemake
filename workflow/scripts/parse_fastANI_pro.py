"""
Purpose: to obtain fastANI results between Pro bins. 
"""
import pandas as pd 

### import paths to data tables ### 
aggregate_table = snakemake.input['aggregate_table']  # table of: bin, classification, quality 
fastANI_fpath = snakemake.input['fastANI_result']  # table of many-to-many fastANI result 
output = snakemake.output[0]  # path to output final .tsv file 

def import_fastANI(fpath):
    """
    """
    df = pd.read_table(
        fpath, names=['query', 'reference', 'ANI', 'fragment_mapping', 'total_query_fragment'], 
    )

    cols = ['query', 'reference']
    for col in cols:
        df[col] = df[col].str.split('/').str[-1].str.replace('.fa', '')

    return df

def import_aggregate(fpath):
    """
    """
    df = pd.read_table(fpath)

    # create columns for mapping with fastANI 
    df['query'] = df['sample'] + '.' + df['bin_id'].astype(str)
    df['reference'] = df['query']

    df['query_genus'] = df['genus']
    df['reference_genus'] = df['genus']

    return df


def main():
    # import and pre-process fastANI and bin/classification data
    fastANI_df = import_fastANI(fastANI_fpath)
    aggregate_df = import_aggregate(aggregate_table)

    # obtain query data
    df = pd.merge(fastANI_df, aggregate_df[['query', 'query_genus']], on=['query'], how='inner')

    # obtain refernce data
    df = pd.merge(df, aggregate_df[['reference', 'reference_genus']], on=['reference'], how='inner')

    # filter for Pro mappings only 
    df = df[
        (df['query_genus'].str.contains('Prochlorococcus')) & 
        (df['reference_genus'].str.contains('Prochlorococcus'))
    ]

    # sort ANI
    df = df.sort_values(by=['ANI'], ascending=False)

    df.to_csv(output, sep='\t', index=False)

main()