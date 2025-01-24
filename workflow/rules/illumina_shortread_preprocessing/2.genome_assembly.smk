rule assembly_spades:
    """
    Short read assembly using SPAdes.
    """
    input:
        r1 = scratch_dict["QC"] / "{sample}_1_trimmed.fastq.gz",
        r2 = scratch_dict["QC"] / "{sample}_2_trimmed.fastq.gz",
    output:
        scratch_dict["genome_assembly"] / "{sample}" / "scaffolds.fasta",
    conda:
        "../../envs/spades.yaml"
    shell:
        "spades.py --threads {resources.cpus_per_task} -1 {input.r1} -2 {input.r2} -o $(dirname {output})"