rule combine_results:
    """
    Combine binning, quality, and classification results. 
    """
    input:
        bin_dirs = expand(scratch_dict["metabat_binning"] / "bins" / "{sample}", sample=SAMPLES), 
        bin_depths = expand(scratch_dict["metabat_binning"] / "depth" / "{sample}_depth.txt", sample=SAMPLES),
        qualities = expand(scratch_dict["checkm_bin_quality"] / "{sample}" / "quality_report.tsv", sample=SAMPLES), 
        classifications = expand(scratch_dict["gtdb_classification"]  / "{sample}" / "gtdbtk.bac120.summary.tsv", sample=SAMPLES), 
    output:
        final_tsv = results_dict['final_tsv'], 
        final_tsv_filtered = results_dict['final_tsv_filtered'], 
    conda: 
        "../envs/data.yaml"
    script:
        "../scripts/combine_results.py"

rule checkM2_plot:
    """
    Plot distribution of Completeness vs. Contamination of all bins. 
    """
    input: 
        quals = expand(scratch_dict["checkm_bin_quality"] / "{sample}" / "quality_report.tsv", sample=SAMPLES), 
    output: 
        plot_out = results_dict['bin_quality_plot'], 
    conda: "../envs/data.yaml"
    script: "../scripts/plotting/checkM2_plot.py"

rule classification_plot:
    """
    Plot classification of all bins in samples. 
    """
    input: 
        table = results_dict['final_tsv'], 
    output: 
        plot_out = results_dict['classification_plot'], 
    conda: "../envs/data.yaml"
    script: "../scripts/plotting/GTDB-classification_plot.py"