rule annotate_gene:
    input: 
        scratch_dict["metabat_binning"] / "bins" / "{sample}",
    output:
        scratch_dict['gene_annotation'] / "{sample}" / "{sample}_done.txt",
    conda: 
        "../envs/prokka.yaml"
    shell:
        """
        # run prokka on each bin in the directory 
        for file in {input}/*; do 
            echo Running prokka on: $file

            prokka \
                --cpus {resources.cpus_per_task} \
                --force \
                --outdir $(dirname {output}) \
                $file
        done

        touch {output}
        """
