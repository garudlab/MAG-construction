# Adrian Verster, July 2025

QC_DIR = config.get("scratch_dir") if config.get("scratch_dir") is not None else PROJECT_DIR

def get_processed_reads(sample_or_wildcards):
    sample = sample_or_wildcards.sample if hasattr(sample_or_wildcards, 'sample') else sample_or_wildcards
    if config.get("skip_qc", False):
        return {
            "fwd": join(FASTQ_INDIR, f"{sample}_1.fastq.gz"),
            "rev": join(FASTQ_INDIR, f"{sample}_2.fastq.gz")
        }
    if config.get("remove_host", True):
        return {
            "fwd": join(QC_DIR, f"01_processing/03_dehost/{sample}_R1.fq.gz"),
            "rev": join(QC_DIR, f"01_processing/03_dehost/{sample}_R2.fq.gz")
        }
    else:
        return {
            "fwd": join(QC_DIR, f"01_processing/02_trimmed/{sample}_R1_val_1.fq.gz"),
            "rev": join(QC_DIR, f"01_processing/02_trimmed/{sample}_R2_val_2.fq.gz")
        }

rule spades_coassembly:
    input: 
        fwd = lambda wildcards: [get_processed_reads(s)["fwd"] for s in batch_samples[wildcards.batch]],
        rev = lambda wildcards: [get_processed_reads(s)["rev"] for s in batch_samples[wildcards.batch]]
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


rule spades_single:
    input:
        unpack(get_processed_reads)
    output:
        contigs = join(PROJECT_DIR, "02_metaspades/{sample}/contigs.fasta"),
        scaffolds = join(PROJECT_DIR, "02_metaspades/{sample}/scaffolds.fasta")
    resources:
        time = lambda wildcards, attempt: 24 * attempt,
        mem = lambda wildcards, attempt: 64 * attempt,
    threads: 16
    singularity: "docker://quay.io/biocontainers/spades:3.15.5--h95f258a_1"
    params:
        outdir = join(PROJECT_DIR, "02_metaspades/{sample}/"),
        tmp_dir = join(PROJECT_DIR, "tmp/spades_{sample}")
    shell: """
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        mkdir -p {params.tmp_dir}
        spades.py --meta \
            -1 {input.fwd} \
            -2 {input.rev} \
            -o {params.outdir} \
            -m {resources.mem} \
            -t {threads} \
            --tmp-dir {params.tmp_dir} \
            --only-assembler
        rm -rf {params.tmp_dir}
    """



rule quast_spades:
    input:
        join(PROJECT_DIR, "02_metaspades/{batch_or_sample}/contigs.fasta")
    output:
        join(PROJECT_DIR, "02_metaspades/{batch_or_sample}/quast/report.tsv")
    singularity: "docker://quay.io/biocontainers/quast:5.0.2--py35pl526ha92aebf_0"
    params:
        outdir = join(PROJECT_DIR, "02_metaspades/{batch_or_sample}/quast"),
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
        expand(join(PROJECT_DIR, "02_metaspades/{batch_or_sample}/quast/report.tsv"), batch_or_sample=batch_or_sample_all)
    output:
        join(PROJECT_DIR, "02_metaspades/quast_report_merged.tsv")
    params:
        batch_or_sample_names = batch_or_sample_all,
        assembly_dir = join(PROJECT_DIR, "02_metaspades/")
    run:
        import pandas as pd
        
        combined_data = []
        for i, report_file in enumerate(input):
            df = pd.read_csv(report_file, sep='\t')
            if config['coassembly']:
                df['Batch'] = params.batch_or_sample_names[i]
            else:
                df['Sample'] = params.batch_or_sample_names[i]
            combined_data.append(df)
        
        combined_df = pd.concat(combined_data, ignore_index=True)
        combined_df.to_csv(output[0], sep='\t', index=False)
