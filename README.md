# Step 1: VCF Processing Workflow

This step performs quality control and sample subsetting on VCF files using Snakemake workflow management system.

## Overview

The workflow processes VCF files by:
1. Subsetting samples based on a provided sample list
2. Applying quality control filters (MAF and depth)
3. Generating QC reports for each chromosome
4. Creating a final summary report

## Prerequisites

- Miniconda3 (version 24.1.2-py310)
- Active conda environment with bcftools
- SLURM
- Sufficient storage space for temporary files and outputs

## Directory Structure

```
workflow/step1/
├── cluster_config.yaml  # SLURM cluster configuration
├── config.yaml         # Pipeline configuration
├── submit_snakemake.sh # Job submission script
└── vcf_processing.snakemake # Main workflow file
```

## Configuration

### Pipeline Configuration (config.yaml)

Key parameters in `config.yaml`:
- `input_dir`: Directory containing input VCF files
- `output_dir`: Directory for processed VCF files
- `samples_file`: File containing list of samples to extract
- `min_maf`: Minimum minor allele frequency (default: 0.01)
- `min_dp`: Minimum depth threshold (default: 7)
- `threads`: Number of threads for parallel processing
- `memory_gb`: Memory allocation in GB
- `temp_dir`: Directory for temporary files

### Cluster Configuration (cluster_config.yaml)

SLURM resource allocations:
- Default job settings:
  - Time: 4 hours
  - Memory: 32GB
  - CPUs: 24
- Specific rules:
  - `subset_vcf`: 12 hours, 128GB RAM
  - `create_summary`: 2 hours, 16GB RAM

## Usage

1. Ensure all configuration files are properly set up
2. Create required directories:
   ```bash
   mkdir -p logs qc_reports
   ```
3. Submit the workflow:
   ```bash
   sbatch submit_snakemake.sh
   ```

## Output Files

The workflow generates:
- Processed VCF files: `{output_dir}/chr{N}.subset.vcf.gz`
- VCF indexes: `{output_dir}/chr{N}.subset.vcf.gz.tbi`
- QC reports: `qc_reports/chr{N}.qc_report.txt`
- Summary report: `summary_report.txt`

## Quality Control

The workflow applies the following QC filters:
- Minor Allele Frequency (MAF) ≥ 0.01
- Depth (DP) ≥ 7
- Generates comprehensive QC reports for each chromosome

## Error Handling

- The workflow includes error checking and logging
- Failed jobs are recorded in SLURM error logs: `logs/slurm_{rule}_{wildcards}.err`
- Uses `--keep-going` flag to continue with independent jobs if one fails
- Implements `--rerun-incomplete` to handle interrupted jobs

## Resource Management

- Implements efficient resource allocation through SLURM
- Uses temporary directory for intermediate files
- Configurable thread and memory usage per rule
- Includes 180-second latency wait for filesystem sync

## Notes

- Input VCF files must follow the naming convention: `ALL.chr{N}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.vcf.gz`
- Requires indexed VCF files (.tbi)
