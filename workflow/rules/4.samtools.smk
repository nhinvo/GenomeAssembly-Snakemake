rule samtools_view:
    """
    Convert sam to bam.
    """
    input: scratch_dict["read_mapping"] / "{sample}.sam", 
    output: temp(scratch_dict["read_mapping"] / "{sample}_unsorted.bam"), 
    conda: "../envs/samtools.yaml"
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
        scratch_dict["read_mapping"] / "{sample}_unsorted.bam", 
    output: 
        sorted_bam = scratch_dict["read_mapping"] / "{sample}_sorted.bam", 
        indexed_bam = scratch_dict["read_mapping"] / "{sample}_sorted.bam.bai", 
    conda: 
        "../envs/samtools.yaml"
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

rule samtools_mapping_stats:
    """
    Obtain mapping statistics.
    """
    input: 
        scratch_dict["read_mapping"] / "{sample}_sorted.bam", 
    output: 
        idxstats = scratch_dict["read_mapping"] / "{sample}_idxstats.tsv", 
        stats = scratch_dict["read_mapping"] / "{sample}_stats.tsv", 
    conda: 
        "../envs/samtools.yaml"
    shell:
        """
        # run idxstats to obtain contig mapping count 
        samtools idxstats \
            --threads {resources.cpus_per_task} \
            {input} > {output.idxstats}

        # run stats to obtain total read count 
        samtools stats \
            --threads {resources.cpus_per_task} \
            {input} | grep ^SN | cut -f 2- \
            > {output.stats}
        """

rule aggregate_mapping_stats:
    input: expand(scratch_dict["read_mapping"] / "{sample}_stats.tsv", sample=SAMPLES)
    output: results_dict['mapping_stats'],
    conda: "../envs/data.yaml"
    shell: "touch {output}"  # temp 
    # script: "../script/mapping_stats.py"