rule run_trim_PE:
    input:
        r1 = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'forward read'],
        r2 = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'reverse read'],
        ref = config["input"]["adapter_file"],
    output:
        o1 = scratch_dict["QC"] / "{sample}_1_trimmed.fastq.gz",
        o2 = scratch_dict["QC"] / "{sample}_2_trimmed.fastq.gz",
    conda:
        "../../envs/bbtools.yaml"
    shell:
        "bbduk.sh threads={resources.cpus_per_task} "
        "in1={input.r1} in2={input.r2} "
        "out1={output.o1} out2={output.o2} "
        "minlen=25 qtrim=rl trimq=10 "
        "ref={input.ref} ktrim=r k=23 mink=11 hdist=1"