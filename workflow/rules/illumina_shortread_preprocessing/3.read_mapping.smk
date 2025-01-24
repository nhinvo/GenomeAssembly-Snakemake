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

rule samtools_view:
    """
    Convert sam to bam.
    """
    input: scratch_dict["read_mapping"] / "{sample}.sam", 
    output: temp(scratch_dict["read_mapping"] / "{sample}.bam"), 
    conda: "../../envs/samtools.yaml"
    shell:
        """
        samtools view \
            --bam \
            --threads {resources.cpus_per_task} \
            --output {output} \
            {input}
        """

rule samtools_sort_index:
    """
    Sort index bam. 
    """
    input: 
        scratch_dict["read_mapping"] / "{sample}.bam", 
    output: 
        sorted_bam = scratch_dict["read_mapping"] / "{sample}_sorted.bam", 
        indexed_bam = scratch_dict["read_mapping"] / "{sample}_sorted.bam.bai", 
    conda: 
        "../../envs/samtools.yaml"
    shell:
        """
        # sort bam 
        samtools sort \
            --threads {resources.cpus_per_task} \
            -o {output.sorted_bam} \
            {input}

        # index the sorted_bam
        samtools index \
            --threads {resources.cpus_per_task} \
            --bai \
            --output {output.indexed_bam} \
            {output.sorted_bam}
        """