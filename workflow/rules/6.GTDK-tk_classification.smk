rule classify_bins_gtdb_tk:
    """
    Taxonomic classification of assemblies using GTDB-tk
    """
    input: 
        bin_dir = scratch_dict["metabat_binning"] / "bins" / "{sample}", 
        connect = scratch_dict["metabat_binning"] / "depth" / "{sample}_depth.txt",  # input to connect rule to binning 
    output: 
        scratch_dict["gtdb_classification"]  / "{sample}" / "gtdbtk.bac120.summary.tsv", 
    conda: 
        "../envs/GTDB-tk.yaml"
    params:
        GTDB_database = config['database']['GTDB database'], 
    shell:
        """
        # set database path for GTDB-tk 
        export GTDBTK_DATA_PATH={params.GTDB_database}; 

        # run GTDB-tk classifier 
        gtdbtk classify_wf \
            --cpus {resources.cpus_per_task} \
            --genome_dir {input.bin_dir} \
            --extension fa \
            --skip_ani_screen \
            --out_dir $(dirname {output})
        """