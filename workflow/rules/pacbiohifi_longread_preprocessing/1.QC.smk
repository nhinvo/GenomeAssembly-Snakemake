rule bam_to_fastq:
    """
    Convert hifi bam files into fastq for downstream filtering/mapipng. 
    """
    input: lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'bam_path'],
    output: temp(scratch_dict["QC"] / "{sample}.fastq"),
    conda: "../../envs/samtools.yaml"
    shell:
        """
        samtools fastq {input} \
            --threads {resources.cpus_per_task} \
            > {output}
        """

rule read_filtering:
    """
    Filter raw reads by length and mean quality.

    min_length: minimum read length 
    min_mean_q: minimum mean quality threshold (note: NOT Phred score)
        - refer to filtlong documentation for more info on how mean is calculated 
    """
    input: scratch_dict["QC"] / "{sample}.fastq",
    output: temp(scratch_dict["QC"] / "{sample}_filtered.fastq"),
    conda: "../../envs/filtlong.yaml"
    shell:
        """
        filtlong \
            --min_length 1000 \
            --min_mean_q 70 \
            {input} > {output}
        """