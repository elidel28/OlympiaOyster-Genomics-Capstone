#!/usr/bin/env python3
"""
PLINK Filtering Pipeline for SNP Data

This script performs:
1. Converts a SNP-only VCF to PLINK binary format
2. Computes MAF and missingness
3. Filters SNPs by MAF and genotype missingness
4. Performs LD pruning
5. Deduplicates SNP IDs from prune list
6. Extracts high-quality SNPs and exports to VCF
7. Creates a BED file of final SNP positions
"""

import os
import subprocess
from collections import Counter
import yaml

def load_config(config_path="config/config.yaml"):
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

def run_plink_conversion(vcf_input, prefix):
    subprocess.run([
        "plink", "--vcf", vcf_input,
        "--allow-extra-chr",
        "--make-bed",
        "--out", prefix
    ], check=True)

def run_plink_stats(bfile_prefix, stats_prefix):
    subprocess.run([
        "plink", "--bfile", bfile_prefix,
        "--allow-extra-chr",
        "--freq", "--missing",
        "--out", stats_prefix
    ], check=True)

def run_plink_filtering(bfile_prefix, filtered_prefix, maf, geno):
    subprocess.run([
        "plink", "--bfile", bfile_prefix,
        "--allow-extra-chr",
        "--maf", str(maf),
        "--geno", str(geno),
        "--make-bed",
        "--out", filtered_prefix
    ], check=True)

def run_ld_pruning(filtered_prefix, prune_prefix, window, step, threshold):
    subprocess.run([
        "plink", "--bfile", filtered_prefix,
        "--allow-extra-chr",
        "--indep-pairwise", str(window), str(step), str(threshold),
        "--out", prune_prefix
    ], check=True)

def deduplicate_snps(prune_in_file, deduped_extract):
    seen = set()
    with open(prune_in_file, "r") as infile, open(deduped_extract, "w") as outfile:
        for line in infile:
            snp = line.strip()
            if snp not in seen:
                seen.add(snp)
                outfile.write(snp + "\n")
    snps = [line.strip() for line in open(prune_in_file)]
    dupes = [snp for snp, count in Counter(snps).items() if count > 1]
    print(f"Found {len(dupes)} duplicated SNP IDs (removed).")

def extract_snps_to_vcf(filtered_prefix, deduped_extract, vcf_output_prefix):
    subprocess.run([
        "plink", "--bfile", filtered_prefix,
        "--allow-extra-chr",
        "--extract", deduped_extract,
        "--recode", "vcf",
        "--out", vcf_output_prefix
    ], check=True)

def create_bed_file(bim_file, bed_output):
    with open(bim_file) as bim, open(bed_output, 'w') as bed:
        for line in bim:
            fields = line.strip().split()
            chrom = fields[0]
            pos = int(fields[3])
            bed.write(f"{chrom}\t{pos - 1}\t{pos}\n")

def main():
    config = load_config()
    vcf_input = config["vcf_input"]
    work_dir = config["work_dir"]
    os.makedirs(work_dir, exist_ok=True)

    maf = config.get("maf_threshold", 0.1)
    geno = config.get("geno_threshold", 0.05)
    window = config.get("ld_window", 50)
    step = config.get("ld_step", 5)
    threshold = config.get("ld_threshold", 0.2)

    plink_prefix = f"{work_dir}/cohort_plink"
    stats_prefix = f"{work_dir}/cohort_stats"
    filtered_prefix = f"{work_dir}/cohort_filtered"
    prune_prefix = f"{work_dir}/cohort_prune"
    deduped_extract = f"{prune_prefix}.prune.in.dedup"
    final_vcf_prefix = f"{work_dir}/cohort_filtered_snps"
    bed_output = f"{work_dir}/cohort_filtered.snps.bed"
    bim_file = f"{filtered_prefix}.bim"

    run_plink_conversion(vcf_input, plink_prefix)
    run_plink_stats(plink_prefix, stats_prefix)
    run_plink_filtering(plink_prefix, filtered_prefix, maf, geno)
    run_ld_pruning(filtered_prefix, prune_prefix, window, step, threshold)
    deduplicate_snps(f"{prune_prefix}.prune.in", deduped_extract)
    extract_snps_to_vcf(filtered_prefix, deduped_extract, final_vcf_prefix)
    create_bed_file(bim_file, bed_output)
    print(f"Done. Final VCF: {final_vcf_prefix}.vcf and BED: {bed_output}")

if __name__ == "__main__":
    main()
