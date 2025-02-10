rule assembly_metaflye:
    """
    Assemble filtered long reads using flye.

    meta: metagenome / uneven coverage mode
    pacbio-hifi: PacBio HiFi reads (<1% error)

    credit: Konnor von Emster. 
    """
    input: scratch_dict["QC"] / "{sample}_filtered.fastq",
    output: scratch_dict["genome_assembly"] / "{sample}" / "assembly.fasta",
    conda: "../../envs/flye.yaml"
    params: lambda wildcards: str(SAMPLE_TABLE.loc[wildcards.sample, 'metagenomic']),  # is sample metagenomic?
    shell: 
        """
        meta_param=""
        if [ "{params}" == "True" ]; then
            meta_param="--meta "
        fi

        flye \
            --threads {resources.cpus_per_task} $meta_param\
            --pacbio-hifi {input} \
            --out-dir $(dirname {output})
        """