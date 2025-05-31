"""
Combine individual GVCF files into a single multi-sample GVCF using GATK's CombineGVCFs.

Steps:
1. Ensure all per-sample GVCFs exist and are indexed.
2. Combine them into a single cohort-level .g.vcf file.
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
            samples.append({
                "sample_id": row["sample_id"]
            })
    return samples


def validate_and_index_gvcfs(sample_ids, gvcf_dir, gatk_path):
    valid_gvcfs = []

    for sample in sample_ids:
        sample_id = sample["sample_id"]
        gvcf_path = os.path.join(gvcf_dir, f"{sample_id}.final.g.vcf")
        idx_path = gvcf_path + ".idx"

        if not os.path.exists(gvcf_path):
            print(f"Missing GVCF: {gvcf_path}")
            continue

        if not os.path.exists(idx_path):
            print(f"Index missing for {sample_id} — indexing now.")
            subprocess.run([
                gatk_path, "IndexFeatureFile", "-I", gvcf_path
            ], check=True)
        else:
            print(f"Index exists for {sample_id}")

        valid_gvcfs.append(gvcf_path)

    return valid_gvcfs


def combine_gvcfs(reference, gvcfs, output_path, gatk_path):
    print(f"Combining {len(gvcfs)} GVCFs into: {output_path}")

    variant_args = []
    for gvcf in gvcfs:
        variant_args.extend(["-V", gvcf])

    subprocess.run([
        gatk_path, "CombineGVCFs",
        "-R", reference,
        "-O", output_path
    ] + variant_args, check=True)


def main():
    config = load_config()
    samples = read_sample_sheet(config["sample_sheet"])

    reference = config["reference"]
    gatk_path = config["gatk_path"]
    gvcf_dir = config["gvcf_dir"]
    combined_dir = config["combined_gvcf_dir"]

    os.makedirs(combined_dir, exist_ok=True)

    gvcfs = validate_and_index_gvcfs(samples, gvcf_dir, gatk_path)

    combined_gvcf_path = os.path.join(combined_dir, "cohort_combined.g.vcf")
    combine_gvcfs(reference, gvcfs, combined_gvcf_path, gatk_path)


if __name__ == "__main__":
    main()
