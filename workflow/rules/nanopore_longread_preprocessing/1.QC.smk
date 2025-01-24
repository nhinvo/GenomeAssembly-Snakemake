rule read_filtering:
    """
    Filter raw reads by length and mean quality.

    min_length: minimum read length 
    min_mean_q: minimum mean quality threshold (note: NOT Phred score)
        - refer to filtlong documentation for more info on how mean is calculated 
    """
    input: lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'read_path'],
    output: scratch_dict["QC"] / "{sample}_filtered.fastq",
    conda: "../../envs/filtlong.yaml"
    shell:
        """
        filtlong \
            --min_length 1000 \
            --min_mean_q 70 \
            {input} > {output}
        """