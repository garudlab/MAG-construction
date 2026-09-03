# Adrian Verster, July 2025

rule quast_bins:
    input:
        bins_dir = join(PROJECT_DIR, "03_binning/DAStool/{batch_or_sample}/bins")
    output:
        touch(join(PROJECT_DIR, "03_binning/DAStool/{batch_or_sample}/quast_all.done"))
    singularity: "docker://quay.io/biocontainers/quast:5.0.2--py35pl526ha92aebf_0"
    resources:
        mem = 8,
        time = 4
    params:
        quast_root = join(PROJECT_DIR, "03_binning/DAStool/{batch_or_sample}/quast"),
        thresholds = "0,10000,50000,100000,250000,500000,1000000,2000000,3000000"
    shell: """
        for bin_fa in {input.bins_dir}/*.fa; do
            bin=$(basename $bin_fa .fa)
            quast.py -o {params.quast_root}/$bin $bin_fa --contig-thresholds {params.thresholds} --fast
        done
    """

rule aggregate_quast:
    input:
        expand(join(PROJECT_DIR, "03_binning/DAStool/{batch_or_sample}/quast_all.done"), batch_or_sample=batch_or_sample_all)
    output:
        join(PROJECT_DIR, "03_binning/quast_summary.csv")
    params:
        dastool_dir = join(PROJECT_DIR, "03_binning/DAStool")
    shell: """
        python lib/aggregate_quast.py {params.dastool_dir} -o {output}
    """


rule prokka_bins:
    input:
        join(PROJECT_DIR, "03_binning/DAStool/{batch_or_sample}/bins/{bin}.fa")
    output:
        join(PROJECT_DIR, "03_binning/DAStool/{batch_or_sample}/prokka/{bin}/{batch_or_sample}_{bin}.gff")
    singularity: "docker://quay.io/biocontainers/prokka:1.14.5--pl526_0"
    resources:
        mem = 48,
        time = lambda wildcards, attempt: 4 * attempt
    threads: 8
    params:
        prokkafolder = join(PROJECT_DIR, "03_binning/DAStool/{batch_or_sample}/prokka/{bin}"),
        prefix = "{batch_or_sample}_{bin}"
    shell: """
        if [[ {wildcards.bin} =~ "unbinned" ]]; then
            mkdir -p {params.prokkafolder}
            touch {output}
        else
            prokka {input} --outdir {params.prokkafolder} --prefix {params.prefix} \
            --centre X --compliant --force --cpus {threads} --noanno
        fi
    """

rule prokka_bins2:
    input:
        join(PROJECT_DIR, "04_dRep/combined_bins/{bin}.fa")
    output:
        join(PROJECT_DIR, "04_dRep/combined_bins_prokka/{bin}/{bin}.gff")
    singularity: "docker://quay.io/biocontainers/prokka:1.14.6--pl5321hdfd78af_5"
    resources:
	    mem_mb = 16000,
	    runtime = 720,
	    cpus_per_task = 8
    threads: 8
    params:
        prokkafolder = join(PROJECT_DIR, "04_dRep/combined_bins_prokka/{bin}/"),
        prefix = "{bin}"
    shell: """
        if [[ {wildcards.bin} =~ "unbinned" ]]; then
            mkdir -p {params.prokkafolder}
            touch {output}
        else
            prokka {input} --outdir {params.prokkafolder} --prefix {params.prefix} \
            --centre X --compliant --force --cpus {threads} --noanno
        fi
    """

BINS, = glob_wildcards(join(PROJECT_DIR, "04_dRep/combined_bins/{bin}.fa"))
rule run_prokka_bins2:
    input:
        expand(
            join(PROJECT_DIR, "04_dRep/combined_bins_prokka/{bin}/{bin}.gff"),
            bin=BINS
        )


rule aragorn_bins:
    input:
        join(PROJECT_DIR, "03_binning/DAStool/{batch_or_sample}/bins/{bin}.fa")
    output:
        join(PROJECT_DIR, "03_binning/DAStool/{batch_or_sample}/rna/trna/{bin}.txt")
    singularity: "docker://quay.io/biocontainers/prokka:1.14.5--pl526_0"
    resources:
        mem = 8,
        time = 1
    shell: """
        mkdir -p $(dirname {output})
        aragorn -t {input} -o {output}
    """

rule barrnap_bins:
    input:
        join(PROJECT_DIR, "03_binning/DAStool/{batch_or_sample}/bins/{bin}.fa")
    output:
        join(PROJECT_DIR, "03_binning/DAStool/{batch_or_sample}/rna/rrna/{bin}.txt")
    singularity: "docker://quay.io/biocontainers/prokka:1.14.5--pl526_0"
    resources:
        mem = 8,
        time = 1
    shell: """
        mkdir -p $(dirname {output})
        barrnap {input} > {output}
    """

rule kraken2_assembly_strains:
    input:
        join(PROJECT_DIR, "04_dRep/dereplicated_genomes")
    output:
        krak = join(PROJECT_DIR, "04_dRep/kraken.krak"),
        krak_report = join(PROJECT_DIR, "04_dRep/kraken.krak.report")
    singularity: "docker://quay.io/biocontainers/kraken2:2.0.9beta--pl526hc9558a2_0"
    params:
        db = config["kraken2db"]
    resources:
        mem = 256,
        time = 24
    threads: 32
    shell: """
        kraken2 --db {params.db} --threads {threads} \
        --output {output.krak} --report {output.krak_report} {input}/*.fa
    """

rule build_contig_taxonomy_mapping:
    input:
        genomes_dir = join(PROJECT_DIR, "04_dRep/dereplicated_genomes"),
        taxonomy_file = join(config["kraken2db"], "ktaxonomy.tsv")
    output:
        contig_mapping = join(PROJECT_DIR, "04_dRep/contig_to_bin.pkl"),
        taxonomy_mapping = join(PROJECT_DIR, "04_dRep/taxonomy_mapping.pkl")
    resources:
        mem = 16,
        time = 4
    shell: """
        python lib/build_contig_mapping.py \
        --genomes_dir {input.genomes_dir} \
        --taxonomy_file {input.taxonomy_file} \
        --contig-mapping-out {output.contig_mapping} \
        --taxonomy-mapping-out {output.taxonomy_mapping}
    """


rule parse_kraken_annotation:
    input:
        kraken_file = join(PROJECT_DIR, "04_dRep/kraken.krak"),
        contig_mapping = join(PROJECT_DIR, "04_dRep/contig_to_bin.pkl"),
        taxonomy_mapping = join(PROJECT_DIR, "04_dRep/taxonomy_mapping.pkl")
    output:
        join(PROJECT_DIR, "04_dRep/dereplicated_genomes_species_calls.csv")
    resources:
        mem = 16,
        time = 4
    shell: """
        python lib/parse_kraken_annotation.py \
        --kraken_file {input.kraken_file} \
        --contig_mapping {input.contig_mapping} \
        --taxonomy_mapping {input.taxonomy_mapping} \
        --output_file {output}
    """
