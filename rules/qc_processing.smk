# Adrian Verster, July 2025


QC_DIR = config.get("scratch_dir") if config.get("scratch_dir") is not None else PROJECT_DIR

rule deduplicate:
    input:
        fwd = join(FASTQ_INDIR, "{sample}_1.fastq.gz"),
        rev = join(FASTQ_INDIR, "{sample}_2.fastq.gz")
    output:
        fwd = join(QC_DIR, "01_processing/01_dedup/{sample}_R1.fastq.gz"),
        rev = join(QC_DIR, "01_processing/01_dedup/{sample}_R2.fastq.gz"),
    params:
        outdir = join(QC_DIR, "01_processing/01_dedup/")
    threads: 1
    resources:
        mem = lambda wildcards, attempt: attempt * 16, 
        time = 24
    singularity: "docker://dzs74/htstream"
    benchmark: join(QC_DIR, "01_processing/01_dedup/{sample}_time.txt")
    shell: """
        mkdir -p {params.outdir} && cd {params.outdir}
        hts_SuperDeduper -1 {input.fwd} -2 {input.rev} -f {wildcards.sample} -F
    """

rule trim_galore:
    input:
        fwd = rules.deduplicate.output.fwd,
        rev = rules.deduplicate.output.rev,
    output:
        fwd = join(QC_DIR, "01_processing/02_trimmed/{sample}_R1_val_1.fq.gz"),
        rev = join(QC_DIR, "01_processing/02_trimmed/{sample}_R2_val_2.fq.gz"),
    threads: 2
    resources:
        mem=32,
        time=lambda wildcards, attempt: attempt * 24
    params:
        q_min = 30,
        min_len = 60,
        outdir = join(QC_DIR, "01_processing/02_trimmed/"),
    singularity: "docker://quay.io/biocontainers/trim-galore:0.6.7--hdfd78af_0"
    benchmark: join(QC_DIR, "01_processing/02_trimmed/{sample}_time.txt")
    shell: """
        mkdir -p {params.outdir}
        trim_galore --quality {params.q_min} \
            --length {params.min_len} \
            --output_dir {params.outdir} \
            --paired {input.fwd} {input.rev} \
            --cores {threads}
    """

rule host_removal:
    input:
        fwd = rules.trim_galore.output.fwd,
        rev = rules.trim_galore.output.rev,
    output:
        fwd = join(QC_DIR, "01_processing/03_dehost/{sample}_R1.fq.gz"),
        rev = join(QC_DIR, "01_processing/03_dehost/{sample}_R2.fq.gz"),
    params:
        host_ref = config["host_genome"]
    threads: 16
    resources:
        mem=lambda wildcards, attempt: attempt * 16,
        time=lambda wildcards, attempt: attempt * 24
    container: "samtools_bwa.sif"
    shell: """
        bwa mem -t {threads} {params.host_ref} {input.fwd} {input.rev} | samtools fastq -@ {threads} -t -T BX -f 12 -1 {output.fwd} -2 {output.rev}
    """
