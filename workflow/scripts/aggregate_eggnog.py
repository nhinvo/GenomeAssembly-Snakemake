"""
Purpose: to map bins to prokka .gff outputs. 
"""
import pandas as pd 
from pathlib import Path 

# convert snakemake input into list if string 
bin_annot_paths = snakemake.input[0]
bin_annot_paths = [bin_annot_paths] if isinstance(bin_annot_paths, str) else bin_annot_paths

def main():
    data = []
    for fpath in bin_annot_paths:
        # obtain parent dir
        print(fpath)
        fpath = Path(fpath)
        parent_dir = fpath.parent

        # obtain paths to all gff files 
        for gff_fpath in parent_dir.glob('*.gff'):
            fname = gff_fpath.stem
            file_data = {
                'filename': fname,
                'gff_path': gff_fpath, 
            }

            data.append(file_data)

    df = pd.DataFrame(data)
    print(df)
    df['sample'] = df['filename'].str.split('.').str[0]
    df['bin_id'] = df['filename'].str.split('.').str[-1]

    df.to_csv(snakemake.output[0], sep='\t', index=False)

main()