rule map_reads:
    """
    Index reference assembly and map reads to obtain coverage for binning. 

    Credit: Konnor von Emster.
    """
    input: 
        trimmed_r1 = scratch_dict["QC"] / "{sample}_1_trimmed.fastq.gz",
        trimmed_r2 = scratch_dict["QC"] / "{sample}_2_trimmed.fastq.gz",
        reference_assembly = scratch_dict["genome_assembly"] / "{sample}" / "scaffolds.fasta",
    output: 
        temp(scratch_dict["read_mapping"] / "{sample}.sam"),
    conda: 
        "../../envs/bowtie2.yaml"
    shell: 
        """
        # index reference assembly
        bowtie2-build \
            --threads {resources.cpus_per_task} \
            {input.reference_assembly} \
            {input.reference_assembly}

        # map reads 
        bowtie2 \
            --threads {resources.cpus_per_task} \
            -x {input.reference_assembly} \
            -1 {input.trimmed_r1} \
            -2 {input.trimmed_r2} \
            -S {output} 
        """