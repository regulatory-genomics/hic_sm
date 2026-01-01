# HiC-sm

A Snakemake workflow for processing Hi-C sequencing data from raw FASTQ files to analysis-ready contact matrices and downstream analyses.

## Overview

HiC-sm is a comprehensive, modular Snakemake pipeline designed for Hi-C data processing. It supports multiple alignment strategies, quality control, contact matrix generation, and downstream analyses including loop detection and distance-dependent contact analysis.

## Features

- **Multiple Alignment Strategies**: Supports BWA-MEM, BWA-MEM2, BWA-MEME, Bowtie2 (with rescue mapping), and Chromap
- **Quality Control**: Integrated fastp for read trimming and MultiQC for comprehensive QC reporting
- **Flexible Processing**: 
  - Automatic aggregation of multiple runs per sample
  - Support for paired-end and single-end reads
  - Configurable ligation site handling
- **Contact Matrix Generation**: 
  - Cooler format (.mcool) with multiple resolutions
  - Optional .hic format output
  - Multiple filtering strategies
- **Downstream Analysis**:
  - Mustache loop detection
  - Distance-dependent contact analysis
  - Scaling analysis
- **Modular Architecture**: Clean separation of concerns with domain-specific rule files
- **PEP Integration**: Uses Portable Encapsulated Project (PEP) for sample metadata management

## Requirements

- **Snakemake** >= 8.20.1
- **Conda** or **Mamba** for environment management
- **Python** 3.11+
- Genome reference files (FASTA, BWA/Bowtie2 indices, chromosome sizes)

## Installation

1. Clone the repository:
```bash
git clone <repository-url>
cd hic_sm
```

2. Install Snakemake (if not already installed):
```bash
conda install -c conda-forge -c bioconda snakemake
# or
mamba install -c conda-forge -c bioconda snakemake
```

3. The workflow uses conda environments defined in `workflow/envs/`. These will be automatically created when you run the workflow.

## Quick Start

### 1. Prepare Sample Sheet

Create a CSV file with sample information (see `config/hg38_test_sample.csv` for example):

```csv
sample_name,run,R1,R2,ligation_site,skip_ligation
sample1,1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,GATCGATC,false
sample2,1,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz,GATCGATC,false
```

Required columns:
- `sample_name`: Unique sample identifier
- `run`: Run number (integer)
- `R1`, `R2`: Paths to forward and reverse FASTQ files (or `fastq1`, `fastq2`)

Optional columns:
- `ligation_site`: Ligation site sequence (overrides config default)
- `skip_ligation`: Boolean to skip ligation site processing
- `passqc`: QC flag (1 = pass, 0 = fail)

### 2. Configure Workflow

Edit `config/test_config.yml` to set:
- Genome assembly and reference paths
- Mapper selection (`bwa-mem`, `bwa-mem2`, `bwa-meme`, `bowtie2`, `chromap`)
- Output directories
- Resolution bins
- Filtering options

### 3. Set Up PEP Configuration

Configure your PEP project file (default: `pep/project_config.yaml`) to define:
- Sample table path
- Path variables (e.g., `database_dir`, `process_dir`)

### 4. Run the Workflow

Dry-run to check the workflow:
```bash
snakemake -n
```

Execute the workflow:
```bash
snakemake --use-conda --cores <number_of_cores>
```

For cluster execution:
```bash
snakemake --use-conda --profile <profile> --cores <number_of_cores>
```


## Contributing

Contributions are welcome! Please ensure code follows the modular structure and includes appropriate documentation.
