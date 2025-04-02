# Micro-C Workflow

This Snakemake workflow performs **quality control (QC) and processing of Micro-C sequencing data**, following analysis recommendations from [Dovetail Genomics Micro-C Analysis Guide](https://micro-c.readthedocs.io/en/latest/index.html)

***
## Requirements
This workflow requires:
-   **Snakemake** (for workflow automation)
-   **Conda** (for environment management)

***
## Installation & Setup
If you haven't already installed Conda and Snakemake:
1.  **Install Conda**:\
    Follow instructions at Miniconda installation.

2.  **Install Snakemake**:
```bash
    conda install -c conda-forge -c bioconda snakemake
```

3.  **Activate the Snakemake environment** before running the workflow:
```bash
    conda activate snakemake
```

***
## Input Data Requirements

### 1. Raw Data
- Create new directory for storing raw data `mkdir -p data/raw`
- Copy raw sequencing files into the **`data/raw/`** directory. If you have **write access** to the raw data files, create symbolic links to avoid duplication. If you **do not have write access**, copy the files instead.
-   Only **`.fastq.gz`** files are accepted.
-   If your files are not .gz compressed, use `gzip` or `pigz` tool to compress

### 2. Sample Naming Convention
- Files must be named using the following format: `Sample1_R1.fastq.gz`, `Sample1_R2.fastq.gz`
-   `_R1` → Read 1 of paired-end sequencing
-   `_R2` → Read 2 of paired-end sequencing
- **Incorrect names will cause the workflow to fail.**
- **Manually rename files if needed** before running the workflow.

***
## Configuration Setup
Before running the workflow, update the configuration file: `config/config.yml` and set the correct file paths. See example below:

```
reference_genome: "/home/groups/hoolock2/u0/genomes/ucsc/hg38/indexes/bwa/hg38.fa.gz"
genome_file: "/home/groups/hoolock2/u0/genomes/ucsc/hg38/hg38.genome"
chrsizes: "/home/groups/hoolock2/u0/genomes/ucsc/hg38/hg38.chrom.sizes"
ASMBLY: "hg38"
```
- **`reference_genome`** → Path to the BWA index file (`.fa.gz`). If this is the first time running for a specific genome, you **must generate index file** (see instructions below).
- **`genome_file`** → Path to the **`.genome`** file. If this is the first time running for a specific genome, you **must generate genome file** (see instructions below).
-   **`chrsizes`** → Path to the chromosome sizes (`.chrom.sizes`) file.


### Generate BWA index file
Micro-C data is aligned using **Burrows-Wheeler Aligner (BWA)**.\
Before running the workflow, create a **BWA index** for your reference genome.
Run this command in the directory containing the `.fasta` genome sequence:
```
bwa index hg38.fasta
```
This generates the necessary BWA index files.


### Generating a Genome File
A genome file (`.genome`) is required for downstream analysis. If you don't have one, create it from the `.fai` index file:

```
cut -f1,2 hg38.fa.fai > hg38.genome
```

This file contains chromosome names and sizes in **tab-separated format**.

***
## Running the Workflow
Once everything is set up, execute the Snakemake workflow:

### 1. Dry-Run to Check for Issues
Before running, test for missing files or errors:
```
snakemake --use-conda -np
```

### 2. Run the Workflow (on hoolock2)

If running on hoolock2, start a screen session (so workflow will continue even after exiting session):
```
screen -S microc
```

Reactivate snakemake conda environment:
```
conda activate snakemake
```

Execute the full workflow with the desired number of CPU cores:
```
snakemake --use-conda --cores 32 > $(date +"%y%m%d%H%M%S")_snakemake.out 2>&1
```

***
## Results
all results are located in newly generated results/ folder, including:
- fastqc/multiqc reports
- Mapping results
- Summary stats of mapping results
- .hic files for analyzing and visualizing data
 


