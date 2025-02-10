rule checkm_bin_quality:
    """
    Assess genome quality using checkM2.
    """
    input: 
        bin_dir = scratch_dict["metabat_binning"] / "bins" / "{sample}", 
        connect = scratch_dict["metabat_binning"] / "depth" / "{sample}_depth.txt",  # input to connect rule to binning 
    output: 
        scratch_dict["checkm_bin_quality"] / "{sample}" / "quality_report.tsv"
    conda: 
        "../envs/checkm2.yaml"
    params:
        checkM_db = config['database']['checkM database'], 
    shell:
        """
        checkm2 predict \
            --threads {resources.cpus_per_task} \
            --database_path {params.checkM_db} \
            --force --remove_intermediates \
            --input {input.bin_dir} --extension fa \
            --output-directory $(dirname {output})
        """