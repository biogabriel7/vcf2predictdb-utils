#!/bin/bash
#SBATCH --job-name=vcf_process
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=128GB
#SBATCH --account=PAS2598
#SBATCH -o logs/slurm_%j.out
#SBATCH -e logs/slurm_%j.err

# Exit on error
set -e

# Create logs directory if it doesn't exist
mkdir -p logs

# Load conda
module load miniconda3/24.1.2-py310

# Activate conda environment
source activate bcftools_env

# Print environment information
echo "=== Environment Information ==="
echo "Python version: $(python --version)"
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Working directory: $(pwd)"
echo "=========================="

# Run snakemake
snakemake --snakefile vcf_processing.snakemake \
    --cores ${SLURM_CPUS_PER_TASK} \
    --use-conda \
    --rerun-incomplete \
    --keep-going \
    --latency-wait 180 \
    --jobs 23 \
    --verbose

# Check exit status
if [ $? -ne 0 ]; then
    echo "Error: Snakemake pipeline failed"
    exit 1
fi

# Deactivate conda environment
conda deactivate