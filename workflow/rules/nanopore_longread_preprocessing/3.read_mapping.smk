rule map_reads:
    """
    Map reads to assemblies to obtain coverage.

    -ax map-ont: for Oxford Nanopore reads
    """
    input: 
        filtered_reads = scratch_dict["QC"] / "{sample}_filtered.fastq",
        reference_assembly = scratch_dict["genome_assembly"] / "{sample}" / "assembly.fasta",
    output: temp(scratch_dict["read_mapping"] / "{sample}.sam"),
    conda: "../../envs/minimap2.yaml"
    shell: 
        """
        minimap2 \
            -ax map-ont \
            -t {resources.cpus_per_task} \
            {input.reference_assembly} \
            {input.filtered_reads} \
            -o {output}
        """

rule samtools_view:
    """
    Description: convert sam to bam.
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
    Description: sort index bam. 
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