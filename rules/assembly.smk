rule spades_coassembly:
    input: 
        fwd = lambda wildcards: expand(join(PROJECT_DIR, "01_processing/03_dehost/{sample}_R1.fq.gz"), sample=batch_samples[wildcards.batch]),
        rev = lambda wildcards: expand(join(PROJECT_DIR, "01_processing/03_dehost/{sample}_R2.fq.gz"), sample=batch_samples[wildcards.batch])
    output: 
        contigs = join(PROJECT_DIR, "02_metaspades/{batch}/contigs.fasta"),
        scaffolds = join(PROJECT_DIR, "02_metaspades/{batch}/scaffolds.fasta"),
        concat_r1 = temp(join(PROJECT_DIR, "02_metaspades/{batch}/concat_R1.fq")),
        concat_r2 = temp(join(PROJECT_DIR, "02_metaspades/{batch}/concat_R2.fq"))
    resources:
        time = lambda wildcards, attempt: 72 * attempt,
        mem = lambda wildcards, attempt: 333 * attempt,
    threads: 16
    singularity: "docker://quay.io/biocontainers/spades:3.15.5--h95f258a_1"
    params:
        outdir = join(PROJECT_DIR, "02_metaspades/{batch}/"),
        tmp_dir = join(PROJECT_DIR, "tmp/spades_{batch}")
    shell: """
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        mkdir -p {params.tmp_dir}
        
        zcat {input.fwd} > {output.concat_r1}
        zcat {input.rev} > {output.concat_r2}

        spades.py --meta \
            -1 {output.concat_r1} \
            -2 {output.concat_r2} \
            -o {params.outdir} \
            -m {resources.mem} \
            -t {threads} \
            --tmp-dir {params.tmp_dir} \
            --only-assembler
        
        rm -rf {params.tmp_dir}
    """

rule quast_spades:
    input:
        join(PROJECT_DIR, "02_metaspades/{batch}/contigs.fasta")
    output:
        join(PROJECT_DIR, "02_metaspades/{batch}/quast/report.tsv")
    singularity: "docker://quay.io/biocontainers/quast:5.0.2--py35pl526ha92aebf_0"
    params:
        outdir = join(PROJECT_DIR, "02_metaspades/{batch}/quast"),
        min_contig = 500
    resources:
        mem = 8,
        time = 2
    threads: 2
    shell: """
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        quast.py -o {params.outdir} {input} \
            --fast \
            --min-contig {params.min_contig} \
            -t {threads}
    """

rule combine_spades_quast_reports:
    input:
        expand(join(PROJECT_DIR, "02_metaspades/{batch}/quast/report.tsv"), batch=batch_list)
    output:
        join(PROJECT_DIR, "02_metaspades/quast_report_merged.tsv")
    params:
        batch_names = batch_list,
        assembly_dir = join(PROJECT_DIR, "02_metaspades/")
    run:
        import pandas as pd
        
        combined_data = []
        for i, report_file in enumerate(input):
            df = pd.read_csv(report_file, sep='\t')
            df['Batch'] = params.batch_names[i]
            combined_data.append(df)
        
        combined_df = pd.concat(combined_data, ignore_index=True)
        combined_df.to_csv(output[0], sep='\t', index=False)
