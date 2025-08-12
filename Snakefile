# Adrian Verster, July 2025
# run it with snakemake --profile cluster -j 20 --use-singularity --singularity-args "--bind /u:/u" -k --keep-going --keep-incomplete --rerun-triggers mtime --use-conda -k
# snakemake -j 20 --use-singularity --singularity-args "--bind /u:/u" -k --keep-going --keep-incomplete --rerun-triggers mtime -k --singularity-prefix containers/
from os.path import join, abspath, expanduser
import os
import glob
import sys
import pandas as pd

configfile: "config.yaml"

PROJECT_DIR = config["output_directory"]
FASTQ_INDIR = config["input_directory"]

if PROJECT_DIR[0] == '~':
    PROJECT_DIR = expanduser(PROJECT_DIR)
PROJECT_DIR = abspath(PROJECT_DIR)

def get_samples_and_batches(sample_file, sample_col='sample_id', batch_col='subject_id'):
    df = pd.read_csv(sample_file)
    
    samples = df[sample_col].tolist()
    batch_samples = {}
    
    for _, row in df.iterrows():
        sample_id = row[sample_col]
        batch = row[batch_col].replace(" ","_",-1)
        if batch == "blank":
            continue
        batch_samples.setdefault(batch, []).append(sample_id)
    
    return samples, batch_samples

samples, batch_samples = get_samples_and_batches(config['sample_table'], config['sample_col'], config['batch_col'])
batch_list = list(batch_samples.keys())

print(f"Samples found: {samples}")
print(f"Batches found: {batch_list}")
print(f'Total samples: {len(samples)}')
print(f'Total batches: {len(batch_list)}')

include: "rules/qc_processing.smk"
include: "rules/assembly.smk"
include: "rules/binning.smk"
include: "rules/annotation.smk"

def all_bin_outputs(batch_list):
    bin_outputs = []
    for batch in batch_list:
        bin_outputs.extend([
            join(PROJECT_DIR, f"03_binning/DAStool/{batch}/bins"),
            join(PROJECT_DIR, f"03_binning/DAStool/{batch}/checkm/checkm.tsv")
        ])
    return bin_outputs

def all_annotation_outputs():
    annotation_outputs = []
    for batch in batch_list:
        annotation_outputs.extend([
            join(PROJECT_DIR, f"03_binning/DAStool/{batch}/bins")
        ])
    return annotation_outputs

outfiles_qc = []
outfiles_qc.extend(expand(join(PROJECT_DIR, "01_processing/01_dedup/{sample}_R1.fastq.gz"), sample=samples))
outfiles_qc.extend(expand(join(PROJECT_DIR, "01_processing/02_trimmed/{sample}_R1_val_1.fq.gz"), sample=samples))
outfiles_qc.extend(expand(join(PROJECT_DIR, "01_processing/03_dehost/{sample}_R1.fq.gz"), sample=samples))

outfiles_assembly = []
outfiles_assembly.extend(expand(join(PROJECT_DIR, "02_metaspades/{batch}/contigs.fasta"), batch=batch_list))
outfiles_assembly.extend(expand(join(PROJECT_DIR, "02_metaspades/{batch}/quast/report.tsv"), batch=batch_list))
outfiles_assembly.append(join(PROJECT_DIR, "02_metaspades/quast_report_merged.tsv"))

outfiles_binning = []
outfiles_binning.extend(expand(join(PROJECT_DIR, "03_binning/DAStool/{batch}/completed.txt"), batch=batch_list))
#outfiles_binning.append(join(PROJECT_DIR, "04_dRep/checkm.tsv"))

print(outfiles_binning)

outfiles_classification = []
outfiles_classification.append(join(PROJECT_DIR, "04_dRep/kraken.krak"))
outfiles_drep = [join(PROJECT_DIR, "04_dRep/dereplicated_genomes/")]

outfiles_all = []
outfiles_all.extend(outfiles_qc)
outfiles_all.extend(outfiles_assembly)
outfiles_all.extend(outfiles_binning)
#outfiles_all.extend(outfiles_classification)
#outfiles_all.extend(outfiles_drep)

rule all:
    input:
        outfiles_all
