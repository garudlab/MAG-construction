# Adrian Verster, July 2025

rule bwa_index:
    input:
        join(PROJECT_DIR, "02_metaspades/{batch_or_sample}/contigs.fasta")
    output:
        join(PROJECT_DIR, "03_binning/idx/{batch_or_sample}.fa"),
        join(PROJECT_DIR, "03_binning/idx/{batch_or_sample}.fa.amb"),
        join(PROJECT_DIR, "03_binning/idx/{batch_or_sample}.fa.ann"),
        join(PROJECT_DIR, "03_binning/idx/{batch_or_sample}.fa.bwt"),
        join(PROJECT_DIR, "03_binning/idx/{batch_or_sample}.fa.pac"),
        join(PROJECT_DIR, "03_binning/idx/{batch_or_sample}.fa.sa")
    threads: 1
    resources:
        mem = 8,
        time = 2
    container: "/u/home/a/averster/samtools_bwa.sif"
    shell: """
        cp {input} {output[0]}
        bwa index {output[0]}
    """

rule bwa_align:
    input:
        fwd = join(PROJECT_DIR, "01_processing/02_trimmed/{sample}_R1_val_1.fq.gz"),
        rev = join(PROJECT_DIR, "01_processing/02_trimmed/{sample}_R2_val_2.fq.gz"),
        asm = join(PROJECT_DIR, "03_binning/idx/{batch_or_sample}.fa"),
        amb = join(PROJECT_DIR, "03_binning/idx/{batch_or_sample}.fa.amb"),
        ann = join(PROJECT_DIR, "03_binning/idx/{batch_or_sample}.fa.ann"),
        bwt = join(PROJECT_DIR, "03_binning/idx/{batch_or_sample}.fa.bwt"),
        pac = join(PROJECT_DIR, "03_binning/idx/{batch_or_sample}.fa.pac"),
        sa = join(PROJECT_DIR, "03_binning/idx/{batch_or_sample}.fa.sa")
    output:
        join(PROJECT_DIR, "03_binning/alignments/{batch_or_sample}/{sample}.bam")
    threads: 8
    resources:
        mem = lambda wildcards, attempt: 16 * attempt,
        time = 12
    container: "/u/home/a/averster/samtools_bwa.sif"
    shell: """
        mkdir -p $(dirname {output})
        bwa mem -t {threads} {input.asm} {input.fwd} {input.rev} | \
        samtools sort --threads {threads} -o {output}
        samtools index {output}
    """

rule metabat_depth:
    input:
        lambda wildcards: expand(join(PROJECT_DIR, "03_binning/alignments/{batch_or_sample}/{sample}.bam"), 
            batch_or_sample=wildcards.batch_or_sample, sample=batch_samples[wildcards.batch_or_sample])
    output:
        depth = join(PROJECT_DIR, "03_binning/metabat/{batch_or_sample}/depth.txt"),
        paired = join(PROJECT_DIR, "03_binning/metabat/{batch_or_sample}/paired.txt")
    resources:
        mem = 8,
        time = 6
    singularity: "docker://quay.io/biocontainers/metabat2:2.15--h137b6e9_0"
    shell: """
        mkdir -p $(dirname {output.depth})
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} \
        --pairedContigs {output.paired} --minContigLength 1000 \
        --minContigDepth 1 --percentIdentity 50 \
        {input}
    """

checkpoint metabat:
    input:
        asm = join(PROJECT_DIR, "03_binning/idx/{batch_or_sample}.fa"),
        depth = join(PROJECT_DIR, "03_binning/metabat/{batch_or_sample}/depth.txt")
    output:
        directory(join(PROJECT_DIR, "03_binning/metabat/{batch_or_sample}/bins/"))
    threads: 8
    resources:
        mem = 64,
        time = 24
    singularity: "docker://quay.io/biocontainers/metabat2:2.15--h137b6e9_0"
    params:
        outstring = join(PROJECT_DIR, "03_binning/metabat/{batch_or_sample}/bins/bin")
    shell: """
        metabat2 --seed 1 -t {threads} --unbinned \
        --inFile {input.asm} --outFile {params.outstring} --abdFile {input.depth}

        if [ $(ls {output} | wc -l) == "0" ]; then
            cp {input.asm} {output}/bin.unbinned.fa
        fi

        if [ -f {output}/bin.tooShort.fa ]; then
            if [ $(cat {output}/bin.tooShort.fa | wc -l) == "0" ]; then
                rm {output}/bin.tooShort.fa
            fi
        fi
    """

checkpoint maxbin:
    input:
        contigs = join(PROJECT_DIR, "03_binning/idx/{batch_or_sample}.fa"),
        depth = join(PROJECT_DIR, "03_binning/metabat/{batch_or_sample}/depth.txt")
    output:
        directory(join(PROJECT_DIR, "03_binning/maxbin/{batch_or_sample}/bins"))
    threads: 8
    resources:
        mem = 32,
        time = lambda wildcards, attempt: attempt * 6
    singularity: "docker://quay.io/biocontainers/maxbin2:2.2.7--he1b5a44_1"
    params:
        outfolder = join(PROJECT_DIR, "03_binning/maxbin/{batch_or_sample}/"),
        logfile = join(PROJECT_DIR, "03_binning/maxbin/{batch_or_sample}/maxbin.log"),
        abundance_folder = join(PROJECT_DIR, "03_binning/maxbin/{batch_or_sample}/depth_files"),
        abundance_list = join(PROJECT_DIR, "03_binning/maxbin/{batch_or_sample}/abundance_list.txt"),
        num_samples = lambda wildcards: len(batch_samples[wildcards.batch])
    shell: """
        if [ -d {params.outfolder} ]; then rm -r {params.outfolder}; fi
        mkdir -p {params.outfolder}
        cd {params.outfolder}
        
        mkdir {params.abundance_folder}
        for i in $(seq 1 {params.num_samples}); do
            col=$((2 + i * 2))
            cut -f 1,$col {input.depth} | tail -n +2 > {params.abundance_folder}/$i.txt
        done
        
        ls {params.abundance_folder}/*.txt > {params.abundance_list}
        run_MaxBin.pl -contig {input.contigs} -out maxbin \
        -abund_list {params.abundance_list} -thread {threads} || true

        if $(grep -q "cannot be binned" {params.logfile}); then
            mkdir {params.outfolder}/bins/
            cp {input.contigs} {params.outfolder}/bins/maxbin.unbinned.fasta
        elif ls {params.outfolder}/*.fasta 1> /dev/null 2>&1; then
            mkdir {params.outfolder}/bins/
            mv {params.outfolder}/*.fasta {params.outfolder}/bins/
        else 
            exit 1
        fi
    """

rule concoct_cut:
    input:
        join(PROJECT_DIR, "03_binning/idx/{batch_or_sample}.fa")
    output:
        cut_contigs = join(PROJECT_DIR, "03_binning/concoct/{batch_or_sample}/contigs_10K.fa"),
        bedfile = join(PROJECT_DIR, "03_binning/concoct/{batch_or_sample}/contigs_10K.bed")
    resources:
        mem = 64,
        time = 12
    singularity: "docker://quay.io/biocontainers/concoct:1.1.0--py27h88e4a8a_0"
    shell: """
        mkdir -p $(dirname {output.cut_contigs})
        cut_up_fasta.py {input} -c 10000 -o 0 --merge_last -b {output.bedfile} > {output.cut_contigs}
    """

rule concoct_coverage:
    input:
        bedfile = join(PROJECT_DIR, "03_binning/concoct/{batch_or_sample}/contigs_10K.bed"),
        bams = lambda wildcards: expand(join(PROJECT_DIR, "03_binning/alignments/{batch_or_sample}/{sample}.bam"), 
            batch=wildcards.batch, sample=batch_samples[wildcards.batch])
    output:
        join(PROJECT_DIR, "03_binning/concoct/{batch_or_sample}/coverage_table.tsv")
    resources:
        mem = 64,
        time = 12
    singularity: "docker://quay.io/biocontainers/concoct:1.1.0--py27h88e4a8a_0"
    shell: """
        concoct_coverage_table.py {input.bedfile} {input.bams} > {output}
    """

rule concoct_run:
    input:
        cut_contigs = join(PROJECT_DIR, "03_binning/concoct/{batch_or_sample}/contigs_10K.fa"),
        coverage = join(PROJECT_DIR, "03_binning/concoct/{batch_or_sample}/coverage_table.tsv")
    output:
        join(PROJECT_DIR, "03_binning/concoct/{batch_or_sample}/clustering_gt1000.csv")
    params:
        outdir = join(PROJECT_DIR, "03_binning/concoct/{batch_or_sample}")
    threads: 4
    resources:
        mem = 64, 
        time = 24
    singularity: "docker://quay.io/biocontainers/concoct:1.1.0--py27h88e4a8a_0"
    shell: """
        concoct --composition_file {input.cut_contigs} \
        --coverage_file {input.coverage} -b {params.outdir} --threads {threads}
    """

rule concoct_merge:
    input:
        join(PROJECT_DIR, "03_binning/concoct/{batch_or_sample}/clustering_gt1000.csv")
    output:
        join(PROJECT_DIR, "03_binning/concoct/{batch_or_sample}/clustering_merged.csv")
    resources:
        mem = 64,
        time = 12
    singularity: "docker://quay.io/biocontainers/concoct:1.1.0--py27h88e4a8a_0"
    shell: """
        merge_cutup_clustering.py {input} > {output}
    """

checkpoint concoct_extract_bins:
    input:
        original_contigs = join(PROJECT_DIR, "03_binning/idx/{batch_or_sample}.fa"),
        clustering_merged = join(PROJECT_DIR, "03_binning/concoct/{batch_or_sample}/clustering_merged.csv")
    output:
        directory(join(PROJECT_DIR, "03_binning/concoct/{batch_or_sample}/bins/"))
    resources:
        mem = 64,
        time = 12
    singularity: "docker://quay.io/biocontainers/concoct:1.1.0--py27h88e4a8a_0"
    shell: """
        mkdir -p {output}
        extract_fasta_bins.py {input.original_contigs} {input.clustering_merged} \
        --output_path {output}
    """

def get_metabat_bins(wildcards):
    checkpoint_output = checkpoints.metabat.get(**wildcards).output[0]
    return glob_wildcards(os.path.join(checkpoint_output, "{metabat_bin}.fa")).metabat_bin

def get_maxbin_bins(wildcards):
    checkpoint_output = checkpoints.maxbin.get(**wildcards).output[0]
    return glob_wildcards(os.path.join(checkpoint_output, "{maxbin_bin}.fasta")).maxbin_bin

def get_concoct_bins(wildcards):
    checkpoint_output = checkpoints.concoct_extract_bins.get(**wildcards).output[0]
    return glob_wildcards(os.path.join(checkpoint_output, "{concoct_bin}.fasta")).concoct_bin

rule DAStool:
    input:
        lambda wildcards: expand(join(PROJECT_DIR, "03_binning/metabat/{batch_or_sample}/bins/{metabat_bin}.fa"), 
            metabat_bin=get_metabat_bins(wildcards), batch_or_sample=wildcards.batch_or_sample),
        lambda wildcards: expand(join(PROJECT_DIR, "03_binning/maxbin/{batch_or_sample}/bins/{maxbin_bin}.fasta"), 
            maxbin_bin=get_maxbin_bins(wildcards), batch_or_sample=wildcards.batch_or_sample),
        lambda wildcards: expand(join(PROJECT_DIR, "03_binning/concoct/{batch_or_sample}/bins/{concoct_bin}.fasta"), 
            concoct_bin=get_concoct_bins(wildcards), batch_or_sample=wildcards.batch_or_sample),
        contigs = join(PROJECT_DIR, "03_binning/idx/{batch_or_sample}.fa")
    output: 
        join(PROJECT_DIR, "03_binning/DAStool/{batch_or_sample}/completed.txt")
    #conda: "../envs/das_tool.yaml"
    singularity: "docker://quay.io/biocontainers/das_tool:1.1.3--r41hdfd78af_0"
    params:
        outfolder = join(PROJECT_DIR, "03_binning/DAStool/{batch_or_sample}/"),
        metabat_dir = join(PROJECT_DIR, "03_binning/metabat/{batch_or_sample}/bins/"),
        maxbin_dir = join(PROJECT_DIR, "03_binning/maxbin/{batch_or_sample}/bins/"),
        concoct_dir = join(PROJECT_DIR, "03_binning/concoct/{batch_or_sample}/bins/"),
        metabat_tsv = join(PROJECT_DIR, "03_binning/DAStool/{batch_or_sample}/metabat_scaffold2bin.tsv"),
        maxbin_tsv = join(PROJECT_DIR, "03_binning/DAStool/{batch_or_sample}/maxbin_scaffold2bin.tsv"),
        concoct_tsv = join(PROJECT_DIR, "03_binning/DAStool/{batch_or_sample}/concoct_scaffold2bin.tsv"),
        logfile = join(PROJECT_DIR, "03_binning/DAStool/{batch_or_sample}/_DASTool.log")
    threads: 8
    resources:
        mem = 128,
        time = lambda wildcards, attempt: attempt * 6
    shell: """
        mkdir -p {params.outfolder}
        
        Fasta_to_Scaffolds2Bin.sh -e fa -i {params.metabat_dir} > {params.metabat_tsv}
        Fasta_to_Scaffolds2Bin.sh -e fasta -i {params.maxbin_dir} > {params.maxbin_tsv}
        Fasta_to_Scaffolds2Bin.sh -e fa -i {params.concoct_dir} > {params.concoct_tsv}
        
        DAS_Tool -i {params.metabat_tsv},{params.maxbin_tsv},{params.concoct_tsv} \
        -l metabat,maxbin,concoct -c {input.contigs} -o {params.outfolder} \
        --search_engine diamond --threads {threads} --write_bins 1 --write_unbinned 1 || true

        if $(grep -q "No bins with bin-score >0.5 found" {params.logfile}) || \
           $(grep -q "single copy gene prediction using diamond failed" {params.logfile}); then
            mkdir {params.outfolder}/_DASTool_bins
            cp {input.contigs} {params.outfolder}/_DASTool_bins/unbinned.fa
        fi
        
        touch {output}
    """

checkpoint extract_DAStool:
    input: 
        join(PROJECT_DIR, "03_binning/DAStool/{batch_or_sample}/completed.txt")
    output:
        directory(join(PROJECT_DIR, "03_binning/DAStool/{batch_or_sample}/bins"))
    params:
        old_binfolder = join(PROJECT_DIR, "03_binning/DAStool/{batch_or_sample}/_DASTool_bins"),
        new_binfolder = join(PROJECT_DIR, "03_binning/DAStool/{batch_or_sample}/bins")
    shell: """
        cp -r {params.old_binfolder} {params.new_binfolder}
    """


rule prepare_drep_input:
    input:
        bins_dirs = expand(join(PROJECT_DIR, "03_binning/DAStool/{batch_or_sample}/bins"), batch=batch_or_sample)
    output:
        genome_list = join(PROJECT_DIR, "04_dRep/genome_list.txt"),
        combined_dir = directory(join(PROJECT_DIR, "04_dRep/combined_bins/"))
    params:
        drep_dir = join(PROJECT_DIR, "04_dRep/"),
        combined_dir = join(PROJECT_DIR, "04_dRep/combined_bins/")
    shell: """
        mkdir -p {params.combined_dir}
        rm -f {output.genome_list}
        for bins_dir in {input.bins_dirs}; do
            if [ -d "$bins_dir" ]; then
                batch=$(basename $(dirname $bins_dir))
                find "$bins_dir" -name "*.fa" ! -name "*unbinned*" | while read bin_file; do
                    bin_name=$(basename $bin_file .fa)
                    new_name="${{batch_or_sample}}_${{bin_name}}.fa"
                    cp "$bin_file" "{params.combined_dir}/$new_name"
                    echo "{params.combined_dir}$new_name" >> {output.genome_list}
                done
            fi
        done
    """


rule checkm_DAStool_lineage:
    input:
        combined_dir = join(PROJECT_DIR, "04_dRep/combined_bins/")
    output:
        lineage=join(PROJECT_DIR, "04_dRep/checkm/lineage.ms"),
        dir=join(PROJECT_DIR, "04_dRep/checkm/"),
    singularity: "docker://quay.io/biocontainers/checkm-genome:1.2.4--pyhdfd78af_0"
    resources:
        mem = 32, # was 128, way too much
        time = 12
    threads: 8
    params:
        checkm_data=config["checkm_data"]
    shell: """
        export CHECKM_DATA_PATH={params.checkm_data}
        checkm lineage_wf -t {threads} -x fa --tab_table {input} {output.dir}
    """

rule checkm_DAStool_qa:
    input:
        combined_dir = join(PROJECT_DIR, "04_dRep/combined_bins/"),
        dir=join(PROJECT_DIR, "04_dRep/checkm/"),
        lineage=join(PROJECT_DIR, "04_dRep/checkm/lineage.ms"),
    output:
        file=join(PROJECT_DIR, "04_dRep/checkm.tsv"),
    singularity: "docker://quay.io/biocontainers/checkm-genome:1.2.4--pyhdfd78af_0"
    resources:
        mem = 32, # was 128, way too much
        time = 12
    threads: 8
    params:
        checkm_data=config["checkm_data"]
    shell: """
        export CHECKM_DATA_PATH={params.checkm_data}
        checkm qa -f {output.file} {input.lineage} {input.dir}
    """


rule checkm_to_drep:
    input:
        checkm_results=join(PROJECT_DIR, "04_dRep/checkm.tsv"),
    output:
        checkm_results=join(PROJECT_DIR, "04_dRep/checkm_for_drep.csv"),
    shell:
        """
        python lib/convert_checkm_to_drep.py -i {input} -o {output}
        """


rule dRep_genome_collection_strains:
    input:
        genome_list = join(PROJECT_DIR, "04_dRep/genome_list.txt"),
        combined_dir = join(PROJECT_DIR, "04_dRep/combined_bins/"),
        checkm_results=join(PROJECT_DIR, "04_dRep/checkm_for_drep.csv"),
    output:
        directory(join(PROJECT_DIR, "04_dRep/dereplicated_genomes/"))
    singularity: "docker://quay.io/biocontainers/drep:3.2.0--py_0"
    threads: 16
    resources:
        mem = 128,
        time = 48
    params:
        drep_dir = join(PROJECT_DIR, "04_dRep/"),
        completeness = 75,
        contamination = 10,
        ani_primary = 0.95,
        ani_secondary = 0.999,
    shell: """
        rm -fr {output}
        dRep dereplicate {params.drep_dir} \
        -g {input.genome_list} \
        -comp {params.completeness} -con {params.contamination} \
        --genomeInfo {input.checkm_results} \
        -pa {params.ani_primary} -sa {params.ani_secondary} -p {threads}
    """



