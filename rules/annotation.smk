rule quast_bins:
    input:
        join(PROJECT_DIR, "03_binning/DAStool/{batch}/bins/{bin}.fa")
    output:
        join(PROJECT_DIR, "03_binning/DAStool/{batch}/quast/{bin}/report.tsv")
    singularity: "docker://quay.io/biocontainers/quast:5.0.2--py35pl526ha92aebf_0"
    resources:
        mem = 8,
        time = 1
    params:
        quastfolder = join(PROJECT_DIR, "03_binning/DAStool/{batch}/quast/{bin}/"),
        thresholds = "0,10000,50000,100000,250000,500000,1000000,2000000,3000000"
    shell: """
        quast.py -o {params.quastfolder} {input} --contig-thresholds {params.thresholds} --fast
    """

rule prokka_bins:
    input:
        join(PROJECT_DIR, "03_binning/DAStool/{batch}/bins/{bin}.fa")
    output:
        join(PROJECT_DIR, "03_binning/DAStool/{batch}/prokka/{bin}/{batch}_{bin}.gff")
    singularity: "docker://quay.io/biocontainers/prokka:1.14.5--pl526_0"
    resources:
        mem = 48,
        time = lambda wildcards, attempt: 4 * attempt
    threads: 8
    params:
        prokkafolder = join(PROJECT_DIR, "03_binning/DAStool/{batch}/prokka/{bin}"),
        prefix = "{batch}_{bin}"
    shell: """
        if [[ {wildcards.bin} =~ "unbinned" ]]; then
            mkdir -p {params.prokkafolder}
            touch {output}
        else
            prokka {input} --outdir {params.prokkafolder} --prefix {params.prefix} \
            --centre X --compliant --force --cpus {threads} --noanno
        fi
    """

rule aragorn_bins:
    input:
        join(PROJECT_DIR, "03_binning/DAStool/{batch}/bins/{bin}.fa")
    output:
        join(PROJECT_DIR, "03_binning/DAStool/{batch}/rna/trna/{bin}.txt")
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
        join(PROJECT_DIR, "03_binning/DAStool/{batch}/bins/{bin}.fa")
    output:
        join(PROJECT_DIR, "03_binning/DAStool/{batch}/rna/rrna/{bin}.txt")
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
        time = 6
    threads: 16
    shell: """
        kraken2 --db {params.db} --threads {threads} \
        --output {output.krak} --report {output.krak_report} {input}/*.fa
    """

