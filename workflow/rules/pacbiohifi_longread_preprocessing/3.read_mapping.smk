rule map_reads:
    """
    Map reads to assemblies to obtain coverage.

    -ax map-hifi: for PacBio Hifi reads
    """
    input: 
        filtered_reads = scratch_dict["QC"] / "{sample}_filtered.fastq",
        reference_assembly = scratch_dict["genome_assembly"] / "{sample}" / "assembly.fasta",
    output: temp(scratch_dict["read_mapping"] / "{sample}.sam"),
    conda: "../../envs/minimap2.yaml"
    shell: 
        """
        minimap2 \
            -ax map-hifi \
            -t {resources.cpus_per_task} \
            {input.reference_assembly} \
            {input.filtered_reads} \
            -o {output}
        """