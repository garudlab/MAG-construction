#
from os.path import join, abspath, expanduser
import sys
import pandas as pd

configfile: "config.yaml"

################################################################################
# specify project directories
PROJECT_DIR = config["output_directory"]
FASTQ_INDIR = config["input_directory"]

# convert PROJECT_DIR to absolute path
if PROJECT_DIR[0] == '~':
    PROJECT_DIR = expanduser(PROJECT_DIR)
PROJECT_DIR = abspath(PROJECT_DIR)

def get_samples_and_batches(sample_file):
    """Parse metadata CSV to get individual samples and batch groupings"""
    df = pd.read_csv(sample_file)
    
    samples = df['sample_id'].tolist()
    batch_samples = {}
    
    for _, row in df.iterrows():
        sample_id = row['sample_id']
        batch = row['subject_id'].replace(" ","_",-1)
        if batch == "blank":
            continue
        batch_samples.setdefault(batch, []).append(sample_id)
    
    return samples, batch_samples

# get sample and batch info from metadata
samples, batch_samples = get_samples_and_batches(config['sample_table'])
batch_list = list(batch_samples.keys())

print(f"Samples found: {samples}")
print(f"Batches found: {batch_list}")
print(f'Total samples: {len(samples)}')
print(f'Total batches: {len(batch_list)}')

# define output files - include QC outputs
outfiles_qc = []
outfiles_qc.extend(expand(join(PROJECT_DIR, "01_processing/01_dedup/{sample}_R1.fastq.gz"), sample=samples))
outfiles_qc.extend(expand(join(PROJECT_DIR, "01_processing/02_trimmed/{sample}_R1_val_1.fq.gz"), sample=samples))

outfiles_spades_assembly = expand(join(PROJECT_DIR, "02_metaspades/{batch}/contigs.fasta"), batch=batch_list)
outfiles_spades_quast = expand(join(PROJECT_DIR, "02_metaspades/{batch}/quast/report.tsv"), batch=batch_list)

outfiles_all = []
outfiles_all.extend(outfiles_qc)
outfiles_all.extend(outfiles_spades_assembly)
outfiles_all.extend(outfiles_spades_quast)
outfiles_all.append(join(PROJECT_DIR, "02_metaspades/quast_report_merged.tsv"))

def get_spades_reads_command_batch(batch):
    """Build SPAdes command for batch co-assembly using processed reads"""
    batch_sample_list = batch_samples[batch]
    
    cmd_parts = []
    
    # Add paired-end libraries from processed reads
    for i, sample in enumerate(batch_sample_list, 1):
        r1 = join(PROJECT_DIR, f"01_processing/03_dehost/{sample}_R1.fq.gz")
        r2 = join(PROJECT_DIR, f"01_processing/03_dehost/{sample}_R2.fq.gz")
        cmd_parts.append(f"--pe{i}-1 {r1} --pe{i}-2 {r2}")
    
    return " ".join(cmd_parts)

################################################################################

rule all:
    input:
        outfiles_all

################################################################################
######## QC PROCESSING ##########################################################
################################################################################
rule deduplicate:
    input:
        fwd = join(FASTQ_INDIR, "{sample}_1.fastq.gz"),
        rev = join(FASTQ_INDIR, "{sample}_2.fastq.gz")
    output:
        fwd = join(PROJECT_DIR, "01_processing/01_dedup/{sample}_R1.fastq.gz"),
        rev = join(PROJECT_DIR, "01_processing/01_dedup/{sample}_R2.fastq.gz"),
    params:
        #tmp_fwd = '{sample}_1.fastq.gz',
        #tmp_rev = '{sample}_2.fastq.gz',
        outdir = join(PROJECT_DIR, "01_processing/01_dedup/")
    threads: 1
    resources:
        mem = lambda wildcards, attempt: attempt * 16, 
        time = 24
    singularity: "docker://dzs74/htstream"
    benchmark: join(PROJECT_DIR, "01_processing/01_dedup/{sample}_time.txt")
    shell: """
        mkdir -p {params.outdir} && cd {params.outdir}
        hts_SuperDeduper -1 {input.fwd} -2 {input.rev} -f {wildcards.sample} -F
    """

rule trim_galore:
    input:
        fwd = rules.deduplicate.output.fwd,
        rev = rules.deduplicate.output.rev,
    output:
        fwd = join(PROJECT_DIR, "01_processing/02_trimmed/{sample}_R1_val_1.fq.gz"),
        rev = join(PROJECT_DIR, "01_processing/02_trimmed/{sample}_R2_val_2.fq.gz"),
    threads: 2
    resources:
        mem=32,
        time=lambda wildcards, attempt: attempt * 24
    params:
        q_min = 30,  # minimum quality score of 30
        min_len = 60,  # minimum read length of 60
        outdir = join(PROJECT_DIR, "01_processing/02_trimmed/"),
    singularity: "docker://quay.io/biocontainers/trim-galore:0.6.7--hdfd78af_0"
    benchmark: join(PROJECT_DIR, "01_processing/02_trimmed/{sample}_time.txt")
    shell: """
        mkdir -p {params.outdir}
        trim_galore --quality {params.q_min} \
            --length {params.min_len} \
            --output_dir {params.outdir} \
            --paired {input.fwd} {input.rev} \
            --cores {threads}
    """


rule mouse_removal:
    input:
        fwd = rules.trim_galore.output.fwd,
        rev = rules.trim_galore.output.rev,
    output:
        fwd = join(PROJECT_DIR, "01_processing/03_dehost/{sample}_R1.fq.gz"),
        rev = join(PROJECT_DIR, "01_processing/03_dehost/{sample}_R2.fq.gz"),
    params:
        mouse_ref = config["mouse_genome"]
    threads: 16
    resources:
        mem=lambda wildcards, attempt: attempt * 16,
        time=lambda wildcards, attempt: attempt * 24
    container: "/u/home/a/averster/samtools_bwa.sif"
    shell: """
        bwa mem -t {threads} {params.mouse_ref} {input.fwd} {input.rev} | samtools fastq -@ {threads} -t -T BX -f 12 -1 {output.fwd} -2 {output.rev}
    """


################################################################################
######## SPADES CO-ASSEMBLY ####################################################
################################################################################
rule spades_coassembly:
    input: 
        # Depend on all trimmed reads for the batch
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
        #reads_command = lambda wildcards: get_spades_reads_command_batch(wildcards.batch),
        tmp_dir = join(PROJECT_DIR, "tmp/spades_{batch}")
    shell: """
        # Clean up any previous runs
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        mkdir -p {params.tmp_dir}
        
        # Concatenate all R1 and R2 files
        zcat {input.fwd} > {output.concat_r1}
        zcat {input.rev} > {output.concat_r2}

        # Run SPAdes meta co-assembly
        spades.py --meta \
            -1 {output.concat_r1} \
            -2 {output.concat_r2} \
            -o {params.outdir} \
            -m {resources.mem} \
            -t {threads} \
            --tmp-dir {params.tmp_dir} \
            --only-assembler
        
        # Clean up temp directory
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

################################################################################
######## BINNING RULES (using containerized tools) ############################
################################################################################
rule concoct_binning:
    input:
        contigs = join(PROJECT_DIR, "02_metaspades/{batch}/contigs.fasta")
    output:
        join(PROJECT_DIR, "03_binning/concoct/{batch}/clustering_gt1000.csv")
    singularity: "docker://quay.io/biocontainers/concoct:1.1.0--py37h6bb024c_0"
    params:
        outdir = join(PROJECT_DIR, "03_binning/concoct/{batch}/")
    threads: 8
    shell: """
        mkdir -p {params.outdir}
        
        # Cut up contigs
        cut_up_fasta.py {input.contigs} -c 10000 -o 0 --merge_last -b {params.outdir}/contigs_10K.bed > {params.outdir}/contigs_10K.fa
        
        # Map reads back to contigs (you'll need to add mapping step)
        # concoct_coverage_table.py {params.outdir}/contigs_10K.bed {params.outdir}/mapping/*.bam > {params.outdir}/coverage_table.tsv
        
        # Run CONCOCT
        # concoct --composition_file {params.outdir}/contigs_10K.fa --coverage_file {params.outdir}/coverage_table.tsv -b {params.outdir}/
    """

rule maxbin2_binning:
    input:
        contigs = join(PROJECT_DIR, "02_metaspades/{batch}/contigs.fasta")
    output:
        directory(join(PROJECT_DIR, "03_binning/maxbin2/{batch}/"))
    singularity: "docker://quay.io/biocontainers/maxbin2:2.2.7--h9ee0642_2"
    threads: 8
    shell: """
        mkdir -p {output}
        
        # Run MaxBin2 (you'll need abundance files from read mapping)
        # run_MaxBin.pl -contig {input.contigs} -abund abundance_file.txt -out {output}/maxbin
    """
