rule bin_similarity_fastANI:
    """
    fastANI on all bins. 
    """
    input: 
        list(fpath for fpath in Path(scratch_dict["metabat_binning"] / "bins").glob('*/*fa'))
    output: 
        bin_paths = scratch_dict['fastANI'] / "bin_paths.txt", 
        fastANI_out = scratch_dict['fastANI'] / "fastANI.tsv", 
    conda: "../envs/fastANI.yaml"
    shell:
        """
        # create file of all bins 
        printf "%s\n" {input} > {output.bin_paths}

        # run fastANI - many to many (all bins)
        fastANI \
            --threads {resources.cpus_per_task} \
            --ql {output.bin_paths} \
            --rl {output.bin_paths} \
            -o {output.fastANI_out}
        """

rule parse_fastANI_Prochlorococcus:
    input: 
        aggregate_table = results_dict['aggregate_table'], 
        fastANI_result = scratch_dict['fastANI'] / "fastANI.tsv", 
    output:
        results_dict['fastANI_Prochlorococcus'], 
    conda: 
        "../envs/data.yaml"
    script:
        "../scripts/parse_fastANI_pro.py"