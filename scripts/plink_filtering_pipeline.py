"""
PLINK filtering pipeline for SNP data.

Steps:
1. Convert VCF to PLINK format
2. Calculate allele frequencies and missingness
3. Filter variants by MAF and missingness thresholds
4. Export filtered dataset back to VCF format
"""

import os
import subprocess
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
        "--freq",
        "--missing",
        "--out", stats_prefix
    ], check=True)


def run_plink_filtering(bfile_prefix, filtered_prefix, maf="0.1", geno="0.05"):
    subprocess.run([
        "plink", "--bfile", bfile_prefix,
        "--allow-extra-chr",
        "--maf", maf,
        "--geno", geno,
        "--make-bed",
        "--out", filtered_prefix
    ], check=True)


def run_plink_export_to_vcf(filtered_prefix, vcf_output_prefix):
    subprocess.run([
        "plink", "--bfile", filtered_prefix,
        "--allow-extra-chr",
        "--recode", "vcf",
        "--out", vcf_output_prefix
    ], check=True)


def main():
    config = load_config()

    genotyped_dir = config["genotyped_dir"]
    plink_dir = config["plink_dir"]
    os.makedirs(plink_dir, exist_ok=True)

    vcf_input = os.path.join(genotyped_dir, "cohort_snps.vcf")

    plink_prefix = os.path.join(plink_dir, "cohort_plink")
    stats_prefix = os.path.join(plink_dir, "cohort_stats")
    filtered_prefix = os.path.join(plink_dir, "cohort_filtered")
    final_vcf_prefix = os.path.join(plink_dir, "cohort_filtered_snps")

    run_plink_conversion(vcf_input, plink_prefix)
    run_plink_stats(plink_prefix, stats_prefix)
    run_plink_filtering(plink_prefix, filtered_prefix)
    run_plink_export_to_vcf(filtered_prefix, final_vcf_prefix)

    print(f"Done. Final filtered VCF written to: {final_vcf_prefix}.vcf")


if __name__ == "__main__":
    main()
