rule metabat2_binning:
    """
    Obtain coverage (depth) from mapped bam file, and bin 
    assembled contigs.
    """
    input: 
        assembly = scratch_dict["genome_assembly"] / "{sample}" / f"{assembly_fname}.fasta",
        sorted_bam = scratch_dict["read_mapping"] / "{sample}_sorted.bam", 
    output: 
        bam_depth = scratch_dict["metabat_binning"] / "depth" / "{sample}_depth.txt", 
        bin_dir = directory(scratch_dict["metabat_binning"] / "bins" / "{sample}"), 
    conda: 
        "../envs/metabat2.yaml"
    shell:
        """
        # Generate depth file from bam file 
        jgi_summarize_bam_contig_depths \
            --outputDepth {output.bam_depth} \
            {input.sorted_bam}

        # Run metabat2
        metabat2 \
            --inFile {input.assembly} \
            --numThreads {resources.cpus_per_task} \
            --abdFile {output.bam_depth} \
            --outFile {output.bin_dir}/{wildcards.sample}
        """