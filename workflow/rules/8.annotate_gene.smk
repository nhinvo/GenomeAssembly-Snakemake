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

            filename=$(basename "$file" .fa)

            prokka \
                --cpus {resources.cpus_per_task} \
                --force --prefix $filename \
                --outdir $(dirname {output}) \
                $file

            echo Completed prokka on $filename. 
        done

        touch {output}
        """

rule finalize_annot:
    input: expand(scratch_dict['gene_annotation'] / "{sample}" / "{sample}_done.txt", sample=SAMPLES),
    output: results_dict['gff_path_table'],
    conda: "../envs/data.yaml"
    script: "../scripts/finalize_annot.py"