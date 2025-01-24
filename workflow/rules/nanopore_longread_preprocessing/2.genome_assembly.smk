rule assembly_metaflye:
    """
    Assemble filtered long reads using flye.

    meta: metagenome / uneven coverage mode
    nano-corr: ONT reads that were corrected with other methods (<3% error)

    credit: Konnor von Emster. 
    """
    input: scratch_dict["QC"] / "{sample}_filtered.fastq",
    output: scratch_dict["genome_assembly"] / "{sample}" / "assembly.fasta",
    conda: "../../envs/flye.yaml"
    shell: 
        """
        flye \
            --threads {resources.cpus_per_task} \
            --meta --nano-corr {input} \
            --out-dir $(dirname {output})
        """