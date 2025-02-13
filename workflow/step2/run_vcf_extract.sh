#!/bin/bash
#SBATCH --job-name=vcf_process
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=128GB
#SBATCH --account=PAS2598
#SBATCH --output=process_vcfs_%j.log
#SBATCH --error=process_vcfs_%j.err

# Exit on error
set -e

# Set directories
INPUT_DIR="${1:-/fs/ess/PAS2598/h5/Genomic_Data/vcf_snps_only/vcf_processed}"  # Can be overridden by command line argument
OUTPUT_DIR="/fs/ess/PAS2598/h5/Genomic_Data/vcf_snps_only/predictdb"

echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"

# Verify input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# Load required modules for Cardinal cluster
module load miniconda3/24.1.2-py310

# Activate conda environment
source activate bcftools_env

# Create necessary directories
echo "Creating output directories..."
mkdir -p predictdb_files
mkdir -p "$OUTPUT_DIR"

# Extract sample IDs from the first VCF file
echo "Extracting sample IDs..."
if [ -f "${INPUT_DIR}/chr1.subset.vcf.gz" ]; then
    bcftools query -l "${INPUT_DIR}/chr1.subset.vcf.gz" > predictdb_files/samples.txt || {
        echo "Error: Failed to extract sample IDs from chr1.subset.vcf.gz"
        exit 1
    }
    echo "Sample IDs extracted and saved to predictdb_files/samples.txt"
else
    echo "Error: chr1.subset.vcf.gz not found in ${INPUT_DIR}"
    exit 1
fi

# Verify samples file was created and is not empty
if [ ! -s predictdb_files/samples.txt ]; then
    echo "Error: Sample ID file is empty or was not created"
    exit 1
fi

# Process each VCF file
echo "Starting VCF processing..."
for vcf_file in "${INPUT_DIR}"/chr{1..22}.subset.vcf.gz "${INPUT_DIR}"/chrX.subset.vcf.gz; do
    if [ ! -f "$vcf_file" ]; then
        echo "Warning: $vcf_file not found, skipping..."
        continue
    fi
    
    echo "Processing $vcf_file..."
    
    # Extract chromosome number from filename
    chr=$(echo $vcf_file | sed 's/.*chr\([^.]*\).subset.vcf.gz/\1/')
    
    # Run the extraction script
    python vcf_extract.py \
        --vcf "$vcf_file" \
        --sample-list predictdb_files/samples.txt \
        --out-prefix predictdb_files/chr${chr}
    
    # Check if the processing was successful
    if [ $? -eq 0 ]; then
        echo "Successfully processed $vcf_file"
    else
        echo "Error processing $vcf_file"
        exit 1
    fi
done

# Combine all chromosome files
echo "Combining chromosome files..."

# Initialize combined files with headers
echo "Creating combined files..."
if [ -f predictdb_files/chr1.snp_annotation.txt ]; then
    head -n 1 predictdb_files/chr1.snp_annotation.txt > predictdb_files/combined.snp_annotation.txt
    head -n 1 predictdb_files/chr1.genotype.txt > predictdb_files/combined.genotype.txt
else
    echo "Error: Chr1 output files not found. Cannot create combined files."
    exit 1
fi

# Combine files for each chromosome
for chr in {1..22} X; do
    echo "Processing chromosome $chr files..."
    
    if [ -f "predictdb_files/chr${chr}.snp_annotation.txt" ]; then
        echo "Adding chr${chr} to SNP annotation file..."
        tail -n +2 "predictdb_files/chr${chr}.snp_annotation.txt" >> predictdb_files/combined.snp_annotation.txt
    fi
    
    if [ -f "predictdb_files/chr${chr}.genotype.txt" ]; then
        echo "Adding chr${chr} to genotype file..."
        tail -n +2 "predictdb_files/chr${chr}.genotype.txt" >> predictdb_files/combined.genotype.txt
    fi
done

# Verify combined files were created and are not empty
if [ ! -s predictdb_files/combined.snp_annotation.txt ] || [ ! -s predictdb_files/combined.genotype.txt ]; then
    echo "Error: Combined files are empty or were not created properly"
    exit 1
fi

# Print summary of files processed
echo -e "\nProcessed files:"
ls -lh predictdb_files/

# Print some basic stats
echo -e "\nFinal file sizes:"
ls -lh predictdb_files/combined*

# Count number of variants and samples
echo -e "\nFinal statistics:"
echo "Number of variants: $(tail -n +2 predictdb_files/combined.snp_annotation.txt | wc -l)"
echo "Number of samples: $(($(head -n1 predictdb_files/combined.genotype.txt | tr '\t' '\n' | wc -l) - 1))"

echo "Job completed successfully"