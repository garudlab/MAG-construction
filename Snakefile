# Adrian Verster, July 2025
# snakemake --profile cluster -j 20 --use-singularity --singularity-args "--bind /u:/u" -k --keep-going --keep-incomplete --rerun-triggers mtime -k --singularity-prefix containers/

from os.path import join, abspath, expanduser
import pandas as pd

configfile: "config.yaml"

PROJECT_DIR = config["output_directory"]
FASTQ_INDIR = config["input_directory"]

if PROJECT_DIR[0] == '~':
    PROJECT_DIR = expanduser(PROJECT_DIR)
PROJECT_DIR = abspath(PROJECT_DIR)
QC_DIR = config.get("scratch_dir") if config.get("scratch_dir") is not None else PROJECT_DIR

def get_samples_and_batches(sample_file, sample_col='sample_id', batch_col='subject_id'):
    df = pd.read_csv(sample_file)
    
    samples = df[sample_col].tolist()
    if config["samples_ignore"] is not None:
        samples = [x for x in samples if x not in config["samples_ignore"]]

    batch_samples = {}
    
    for _, row in df.iterrows():
        sample_id = row[sample_col]
        batch = row[batch_col].replace(" ","_",-1)
        if batch == "blank":
            continue
        batch_samples.setdefault(batch, []).append(sample_id)
    
    return samples, batch_samples


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


if config["coassembly"]:
    samples, batch_samples = get_samples_and_batches(config['sample_table'], config['sample_col'], config['batch_col'])
    batch_or_sample_all = list(batch_samples.keys())
    print(f"Samples found: {samples}")
    print(f"Batches found: {batch_or_sample_all}")
    print(f'Total samples: {len(samples)}')
    print(f'Total batches: {len(batch_or_sample_all)}')
else:
    all_files = set(os.listdir(FASTQ_INDIR))
    samples = [f.replace('_1.fastq.gz', '') for f in all_files
                       if f.endswith('_1.fastq.gz') and f.replace('_1.fastq.gz', '_2.fastq.gz') in all_files]
    if config["samples_ignore"] is not None:
        samples = [x for x in samples if x not in config["samples_ignore"]]
    if config["samples_sub"] is not None:
        samples = [x for x in samples if x in config["samples_sub"]]

    batch_or_sample_all = samples
    batch_samples = {}
    for sample in batch_or_sample_all:
        batch_samples[sample] = [sample]
    print(f'Total samples: {len(batch_or_sample_all)}')

include: "rules/qc_processing.smk"
include: "rules/assembly.smk"
include: "rules/binning.smk"
include: "rules/annotation.smk"


outfiles_qc = []
if not config.get("skip_qc",False):
    outfiles_qc.extend(expand(join(QC_DIR, "01_processing/01_dedup/{sample}_R1.fastq.gz"), sample=samples))
    outfiles_qc.extend(expand(join(QC_DIR, "01_processing/02_trimmed/{sample}_R1_val_1.fq.gz"), sample=samples))
    if config.get("remove_host",True): 
        outfiles_qc.extend(expand(join(QC_DIR, "01_processing/03_dehost/{sample}_R1.fq.gz"), sample=samples))

outfiles_assembly = []
outfiles_assembly.extend(expand(join(PROJECT_DIR, "02_metaspades/{batch_or_sample}/contigs.fasta"), batch_or_sample=batch_or_sample_all))
outfiles_assembly.extend(expand(join(PROJECT_DIR, "02_metaspades/{batch_or_sample}/quast/report.tsv"), batch_or_sample=batch_or_sample_all))
outfiles_assembly.append(join(PROJECT_DIR, "02_metaspades/quast_report_merged.tsv"))

outfiles_binning = []
outfiles_binning.extend(expand(join(PROJECT_DIR, "03_binning/DAStool/{batch_or_sample}/completed.txt"), batch_or_sample=batch_or_sample_all))
outfiles_binning.append(join(PROJECT_DIR, "04_dRep/checkm.tsv"))

outfiles_classification = []
outfiles_classification.append(join(PROJECT_DIR, "04_dRep/kraken.krak"))
outfiles_drep = [join(PROJECT_DIR, "04_dRep/dereplicated_genomes/")]

outfiles_all = []
outfiles_all.extend(outfiles_qc)
outfiles_all.extend(outfiles_assembly)
outfiles_all.extend(outfiles_binning)
outfiles_all.extend(outfiles_classification)
outfiles_all.extend(outfiles_drep)

rule all:
    input:
        outfiles_all
