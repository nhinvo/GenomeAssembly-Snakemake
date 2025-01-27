rule aggregate_results:
    """
    Combine binning, quality, and classification results. 
    """
    input:
        bin_dirs = expand(scratch_dict["metabat_binning"] / "bins" / "{sample}", sample=SAMPLES), 
        qualities = expand(scratch_dict["checkm_bin_quality"] / "{sample}" / "quality_report.tsv", sample=SAMPLES), 
        classifications = expand(scratch_dict["gtdb_classification"]  / "{sample}" / "gtdbtk.bac120.summary.tsv", sample=SAMPLES), 
    output:
        results_dict['aggregate_table'], 
    conda: 
        "../envs/data.yaml"
    script:
        "../scripts/aggregate_results.py"