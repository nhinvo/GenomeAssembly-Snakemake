rule annotate_gene_prokka:
    input: 
        scratch_dict["metabat_binning"] / "bins" / "{sample}",
    output:
        scratch_dict['prokka_gene_annotation'] / "{sample}" / "{sample}_done.txt",
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

rule annotate_gene_eggnog:
    input:
        scratch_dict['prokka_gene_annotation'] / "{sample}" / "{sample}_done.txt",
    output:
        scratch_dict['eggnog_gene_annotation'] / "{sample}" / "{sample}_done.txt",
    conda:
        "../envs/eggnog.yaml"
    params:
        eggnogg_db = config['database']['eggnogg database'], 
        temp_dir = scratch_dict['eggnog_gene_annotation'], 
    shell:
        """
        # obtain parent dir of done.txt file path 
        sample_dir=$(dirname {input})
        echo "Sample directory path: $sample_dir"

        # run eggnog on each bin in the sample directory 
        for file in "$sample_dir"/*.faa; do 
            echo "Running eggnog on .faa: $file"

            filename=$(basename "$file" .faa)

            emapper.py \
                --cpu {resources.cpus_per_task} --override \
                -i $file --itype proteins \
                --decorate_gff yes --report_no_hits \
                --temp_dir {params.temp_dir} \
                --data_dir {params.eggnogg_db} \
                --output $filename \
                --output_dir $(dirname {output})

            echo "Completed eggnog on $filename."
        done

        touch {output}
        """

rule aggregate_prokka:
    """
    """
    input: expand(scratch_dict['prokka_gene_annotation'] / "{sample}" / "{sample}_done.txt", sample=SAMPLES),
    output: results_dict['prokka_final_gff'],
    conda: "../envs/data.yaml"
    script: "../scripts/aggregate_prokka.py"

rule aggregate_eggnog:
    """
    """
    input: expand(scratch_dict['eggnog_gene_annotation'] / "{sample}" / "{sample}_done.txt", sample=SAMPLES),
    output: results_dict['eggnog_final_gff'],
    conda: "../envs/data.yaml"
    script: "../scripts/aggregate_eggnog.py"