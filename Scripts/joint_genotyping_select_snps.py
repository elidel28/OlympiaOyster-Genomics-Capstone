#!/usr/bin/env python3
"""
Joint Genotyping and SNP Selection using GATK.

This script performs:
1. Joint genotyping on a combined GVCF
2. SNP-only selection from the joint-genotyped VCF
"""

import os
import subprocess
import yaml
import csv


def load_config(config_path="config/config.yaml"):
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def read_sample_sheet(sample_sheet_path):
    samples = []
    with open(sample_sheet_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            samples.append({"sample_id": row["sample_id"]})
    return samples


def joint_genotyping(config):
    """Run GenotypeGVCFs to produce a joint VCF."""
    gatk = config["gatk_path"]
    ref = config["reference"]
    combined_dir = config["combined_gvcf_dir"]
    genotyped_dir = config["genotyped_dir"]
    ram = config.get("ram", "8g")

    os.makedirs(genotyped_dir, exist_ok=True)

    input_gvcf = os.path.join(combined_dir, "cohort_combined.g.vcf")
    joint_vcf = os.path.join(genotyped_dir, "cohort_genotyped.vcf")

    print(f"Running joint genotyping → {joint_vcf}")
    subprocess.run([
        gatk, "--java-options", f"-Xmx{ram}", "GenotypeGVCFs",
        "-R", ref,
        "-V", input_gvcf,
        "-O", joint_vcf
    ], check=True)

    return joint_vcf


def select_snps(config, joint_vcf):
    """Extract SNPs from the joint-genotyped VCF."""
    gatk = config["gatk_path"]
    ref = config["reference"]
    genotyped_dir = config["genotyped_dir"]
    ram = config.get("ram", "8g")

    snps_vcf = os.path.join(genotyped_dir, "cohort_snps.vcf")

    print(f"Selecting SNPs only → {snps_vcf}")
    subprocess.run([
        gatk, "--java-options", f"-Xmx{ram}", "SelectVariants",
        "-R", ref,
        "-V", joint_vcf,
        "--select-type-to-include", "SNP",
        "-O", snps_vcf
    ], check=True)


def main():
    config = load_config()
    _ = read_sample_sheet(config["sample_sheet"])  # not used but kept for consistency
    joint_vcf = joint_genotyping(config)
    select_snps(config, joint_vcf)


if __name__ == "__main__":
    main()
