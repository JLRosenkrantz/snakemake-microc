## Snakefile for processing Micro-C data
## following analysis recommendations from [Dovetail Genomics Micro-C Analysis Guide](https://micro-c.readthedocs.io/en/latest/index.html)


## conda activate microc_env

configfile: "config/config.yml"

import re
import glob
import os

# Auto-detect samples from raw data directory
raw_fastq_files = glob.glob("data/raw/*.fastq*")
print("Raw fastq files:", raw_fastq_files)

reads = [re.sub(r"\.fastq(\.gz)?$", "", os.path.basename(f)) for f in raw_fastq_files]
print("Detected reads:", reads)

samples = list(set(re.sub(r"_[Rr][12].*", "", os.path.basename(f)) for f in raw_fastq_files))
print("Detected samples:", samples)

#####################################
#### RULES
#####################################

rule all:
    input:
        expand("results/fastqc/{read}.html", read=reads),
        expand("results/mapped/{sample}/{sample}.mapped.PT.bam", sample=samples),
        expand("results/mapped/{sample}/{sample}.duplication_stats_summary.txt", sample=samples),
        expand("results/mapped/{sample}/{sample}.complexity.txt", sample=samples),
        expand("results/hic/{sample}.contact_map.hic", sample=samples),

rule fastqc:
    input:
        "data/raw/{read}.fastq.gz"
    output:
        html="results/fastqc/{read}.html",
        zip="results/fastqc/{read}_fastqc.zip"
    log:
        "results/logs/fastqc_{read}.log"
    params:
        "--threads 6"
    wrapper:
        "v1.5.0/bio/fastqc"

rule multiqc_raw:
    input:
       expand("results/fastqc/{read}_fastqc.zip", read=reads),
    output:
        "results/fastqc/multiqc_report.html",
        directory("results/fastqc/multiqc_data"),
    log:
        "results/logs/multiqc_raw.log"
    wrapper:
        "v5.8.3/bio/multiqc"


# Single-Step Pipeline: FASTQ → Aligned BAM → Filtered Pairs → Final BAM
rule fastq_to_valid_pairs_bam:
    input:
        ref=config["reference_genome"],
        r1="data/raw/{sample}_R1.fastq.gz",
        r2="data/raw/{sample}_R2.fastq.gz",
        genome_file=config["genome_file"],
    output:
        bam="results/mapped/{sample}/{sample}.mapped.PT.bam",
        pairs="results/mapped/{sample}/{sample}.mapped.pairs",
        stats="results/mapped/{sample}/{sample}.duplication_stats.txt"
    conda:
        "envs/microc_env.yaml"
    log:
        "results/logs/valid_pairs_{sample}.log"
    threads: 16
    shell:
        """
        bwa mem -5SP -T0 -t {threads} {input.ref} {input.r1} {input.r2} | \
        pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in {threads} --nproc-out {threads} --chroms-path {input.genome_file} | \
        pairtools sort --nproc {threads} | \
        pairtools dedup --nproc-in {threads} --nproc-out {threads} --mark-dups --output-stats {output.stats} | \
        pairtools split --nproc-in {threads} --nproc-out {threads} --output-pairs {output.pairs} --output-sam - | \
        samtools view -bS -@ {threads} | \
        samtools sort -@ {threads} -o {output.bam}
        
        # Index the final BAM file
        samtools index {output.bam}
        """

# Summarize stats
rule summarize_stats:
    input:
        stats="results/mapped/{sample}/{sample}.duplication_stats.txt"
    output:
        summary="results/mapped/{sample}/{sample}.duplication_stats_summary.txt"
    shell:
        """
        python ./scripts/get_qc.py -p {input.stats} > {output.summary}
        """
        
# Analyze Library Complexity - in this example the output file out.preseq will detail the extrapolated complexity curve of your library, with the number of reads in the first column and the expected distinct read value in the second column. For a typical experiment (human sample) check the expected complexity at 300M reads (to show the content of the file, type cat out.preseq). Expected unique pairs at 300M sequencing is at least ~ 120 million.
rule lib_complexity:
    input:
        bam="results/mapped/{sample}/{sample}.mapped.PT.bam",
    output:
        complexity="results/mapped/{sample}/{sample}.complexity.txt"
    conda:
        "envs/microc_env.yaml"
    shell:
        """
        preseq lc_extrap -bam -pe -extrap 2.1e9 -step 1e8 -seg_len 1000000000 -output {output.complexity} {input.bam}
        """


rule contact_matrix:
    input:
        pairs="results/mapped/{sample}/{sample}.mapped.pairs",
        genome_file=config["genome_file"],
    output:
        cmap="results/hic/{sample}.contact_map.hic",
    conda:
        "envs/microc_env.yaml"
    threads: 16
    shell:
        """
        java -Xmx32G  -Djava.awt.headless=true -jar ./scripts/juicer_tools_1.22.01.jar pre --threads {threads} {input.pairs} {output.cmap} {input.genome_file}
        """

        
rule index_mapped_pairs:
    input:
        pairs="results/mapped/{sample}/{sample}.mapped.pairs",
    output:
        compressed="results/mapped/{sample}/{sample}.mapped.pairs.gz",
        index="results/mapped/{sample}/{sample}.mapped.pairs.gz.px2",
    conda:
        "envs/microc_env.yaml"
    threads: 16
    shell:
        """
        bgzip -@ {threads} -c {input.pairs} > {output.compressed}
        pairix {output.compressed}
        """

# Rule to create .cool files from .mapped.pairs.gz using cooler cload
rule cooler_cload:
    input:
        pairs="results/mapped/{sample}/{sample}.mapped.pairs.gz"
    output:
        cool="results/cool/{sample}.cool"
    conda:
        "envs/hic.yml"
    log:
        "results/logs/cooler_cload.{sample}.log"
    threads: 16
    shell:
        """
        cooler cload pairix -p {threads} --assembly {config[ASMBLY]} {config[CHRSIZES]}:10000 {input.pairs} {output.cool} >{log} 2>&1
        """

# Rule to create multi-resolution .mcool files for visualization in HiGlass
rule zoomify:
    input:
        cool="results/cool/{sample}.cool"
    output:
        mcool="results/cool/{sample}.mcool"
    log:
        "results/logs/zoomify.{sample}.log"
    conda:
        "envs/hic.yml"
    threads: 16 
    shell:
        """
        cooler zoomify -p {threads} --balance -o {output.mcool} {input.cool} >{log} 2>&1
        """
        
        
rule hicFindTADs:
    input:
        mcool="results/cool/{sample}.mcool"
    output:
        boundaries="results/TADs/{sample}_min10_max60_fdr01_d01_boundaries.bed"
    conda:
        "envs/hicexplorer.yml"
    log:
        "results/logs/hicFindTADs_narrow.{sample}.log"
    threads: 16
    shell:
        """
        hicFindTADs -m {input.mcool}::resolutions/10000 \
                    --minDepth 100000 \
                    --maxDepth 600000 \
                    --outPrefix results/TADs/{wildcards.sample}_min10_max60_fdr01_d01 \
                    --correctForMultipleTesting fdr \
                    -p {threads} >{log} 2>&1
        """
        
        
