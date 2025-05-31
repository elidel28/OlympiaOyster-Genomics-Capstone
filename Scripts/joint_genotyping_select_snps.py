"""
Joint Genotyping and SNP Selection using GATK.

This script performs:
1. Joint genotyping on a combined GVCF
2. SNP-only selection from the joint-genotyped VCF
"""

import os
import subprocess
import yaml


def load_config(config_path="config/config.yaml"):
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def joint_genotyping(config):
    """Run GenotypeGVCFs to produce a joint VCF."""
    gatk = config["gatk_path"]
    ref = config["reference"]
    gvcf = config["combined_gvcf"]
    output = config["joint_vcf"]
    ram = config.get("ram", "8g")

    print(f"Running joint genotyping → {output}")
    subprocess.run([
        gatk, "--java-options", f"-Xmx{ram}", "GenotypeGVCFs",
        "-R", ref,
        "-V", gvcf,
        "-O", output
    ], check=True)


def select_snps(config):
    """Extract SNPs from a joint-genotyped VCF."""
    gatk = config["gatk_path"]
    ref = config["reference"]
    joint_vcf = config["joint_vcf"]
    snps_vcf = config["snps_vcf"]
    ram = config.get("ram", "8g")

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
    joint_genotyping(config)
    select_snps(config)


if __name__ == "__main__":
    main()
