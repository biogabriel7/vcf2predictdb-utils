import argparse
import gzip
from pathlib import Path
import pandas as pd
import numpy as np

def parse_vcf_line(line):
    """Parse a VCF data line into components."""
    fields = line.strip().split('\t')
    # Remove 'chr' prefix if present and ensure X/Y are preserved
    chrom = fields[0].replace('chr', '')
    pos = fields[1]
    rsid = fields[2]
    ref = fields[3]
    alt = fields[4]
    info = dict(x.split('=') for x in fields[7].split(';') if '=' in x)
    
    # Extract AF (allele frequency) from INFO field
    maf = float(info.get('AF', '0'))
    # Make sure MAF is <= 0.5
    maf = min(maf, 1 - maf)
    
    # Handle R2 (imputation quality)
    r2 = info.get('R2', info.get('DR2', info.get('AR2', '1.0')))  # Check multiple common R2 field names
    try:
        r2 = float(r2)
        # Ensure R2 is in valid range
        r2 = max(0.0, min(1.0, r2))
    except ValueError:
        r2 = 1.0  # Default to 1.0 if parsing fails
    
    # Create variant ID in required format: chr_pos_ref_alt_b37
    var_id = f"{chrom}_{pos}_{ref}_{alt}_b37"
    
    # Get genotype data - skip format field
    genotypes = fields[9:]
    
    return {
        'chromosome': chrom,
        'pos': pos,
        'varID': var_id,
        'ref_vcf': ref,
        'alt_vcf': alt,
        'R2': r2,  # Use processed R2 value
        'MAF': maf,
        'rsid': rsid,
        'genotypes': genotypes
    }

def process_vcf_file(vcf_path, sample_ids):
    """Process VCF file and return SNP annotations and genotype data."""
    snp_annotations = []
    genotype_data = []
    
    # Open VCF file (handles both .gz and regular files)
    opener = gzip.open if str(vcf_path).endswith('.gz') else open
    with opener(vcf_path, 'rt') as vcf:
        # Skip header lines and get sample IDs
        for line in vcf:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                file_samples = line.strip().split('\t')[9:]
                # Create mapping from file column to desired sample order
                sample_map = {i: file_samples.index(sample) 
                            for i, sample in enumerate(sample_ids)
                            if sample in file_samples}
                break
        
        # Process variant lines
        for line in vcf:
            variant = parse_vcf_line(line)
            
            # Add to SNP annotations
            snp_annotations.append({
                'chromosome': variant['chromosome'],
                'pos': variant['pos'],
                'varID': variant['varID'],
                'ref_vcf': variant['ref_vcf'],
                'alt_vcf': variant['alt_vcf'],
                'R2': variant['R2'],
                'MAF': variant['MAF'],
                'rsid': variant['rsid']
            })
            
            # Process genotypes into dosages (0-2)
            dosages = []
            for sample_idx in range(len(sample_ids)):
                if sample_idx in sample_map:
                    geno = variant['genotypes'][sample_map[sample_idx]].split(':')[0]
                    if geno == '0/0' or geno == '0|0':
                        dosages.append('0')
                    elif geno == '0/1' or geno == '1/0' or geno == '0|1' or geno == '1|0':
                        dosages.append('1')
                    elif geno == '1/1' or geno == '1|1':
                        dosages.append('2')
                    else:
                        dosages.append('NA')
                else:
                    dosages.append('NA')
            
            genotype_data.append([variant['varID']] + dosages)
    
    return pd.DataFrame(snp_annotations), pd.DataFrame(genotype_data)

def validate_outputs(snp_annot_df, genotype_df, expected_samples):
    """Validate the output files against expected metrics."""
    print("\nValidation Results:")
    
    # Check number of variants
    print(f"Number of variants: {len(snp_annot_df)}")
    
    # Check sample count
    sample_count = len(genotype_df.columns) - 1  # -1 for varID column
    print(f"Number of samples: {sample_count} (Expected: {expected_samples})")
    
    # Check MAF distribution
    maf_dist = snp_annot_df['MAF'].describe()
    print("\nMAF Distribution:")
    print(f"Min: {maf_dist['min']:.4f}")
    print(f"Max: {maf_dist['max']:.4f}")
    print(f"Mean: {maf_dist['mean']:.4f}")
    
    # Check R2 distribution
    r2_dist = snp_annot_df['R2'].describe()
    print("\nR2 Distribution:")
    print(f"Min: {r2_dist['min']:.4f}")
    print(f"Max: {r2_dist['max']:.4f}")
    print(f"Mean: {r2_dist['mean']:.4f}")
    
    # Verify chromosome format
    invalid_chroms = snp_annot_df[~snp_annot_df['chromosome'].str.match(r'^([1-9]|1\d|2[0-2]|X)$')]
    if not invalid_chroms.empty:
        print("\nWARNING: Found invalid chromosome formats:")
        print(invalid_chroms['chromosome'].unique())

def main():
    parser = argparse.ArgumentParser(description='Extract SNP annotations and genotype data from VCF')
    parser.add_argument('--vcf', required=True, help='Input VCF file')
    parser.add_argument('--sample-list', required=True, help='File containing sample IDs')
    parser.add_argument('--out-prefix', required=True, help='Output file prefix')
    args = parser.parse_args()
    
    # Read sample IDs
    with open(args.sample_list) as f:
        sample_ids = [line.strip() for line in f if line.strip()]
    
    print(f"Processing VCF file: {args.vcf}")
    print(f"Number of samples in list: {len(sample_ids)}")
    
    try:
        # Process VCF file
        snp_annot_df, genotype_df = process_vcf_file(args.vcf, sample_ids)
        
        # Set column names for genotype data
        genotype_df.columns = ['varID'] + sample_ids
        
        # Validate outputs
        validate_outputs(snp_annot_df, genotype_df, len(sample_ids))
        
        # Save outputs
        out_snp = f"{args.out_prefix}.snp_annotation.txt"
        out_geno = f"{args.out_prefix}.genotype.txt"
        
        print(f"\nSaving outputs:")
        print(f"SNP annotation file: {out_snp}")
        print(f"Genotype file: {out_geno}")
        
        snp_annot_df.to_csv(out_snp, sep='\t', index=False)
        genotype_df.to_csv(out_geno, sep='\t', index=False)
        
        print("Processing completed successfully")
        
    except Exception as e:
        print(f"Error processing VCF file: {str(e)}")
        raise

if __name__ == '__main__':
    main()