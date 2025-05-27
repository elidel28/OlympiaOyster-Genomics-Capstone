"""
Combine individual GVCF files into a single multi-sample GVCF using GATK's CombineGVCFs.

This script:
1. Validates and indexes per-sample .final.g.vcf files if necessary.
2. Combines them into one cohort-level .g.vcf file.
"""

import os
import subprocess
import yaml

# === Load config ===
def load_config(config_path="config/config.yaml"):
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


# === Ensure .idx exists for each GVCF ===
def validate_and_index_gvcfs(sample_ids, gvcf_dir, gatk_path):
    valid_gvcfs = []

    for sample in sample_ids:
        gvcf_path = os.path.join(gvcf_dir, f"{sample}.final.g.vcf")
        idx_path = gvcf_path + ".idx"

        if not os.path.exists(gvcf_path):
            print(f" Missing GVCF: {gvcf_path}")
            continue

        if not os.path.exists(idx_path):
            print(f"Indexing missing for {sample} — indexing now.")
            subprocess.run([
                gatk_path, "IndexFeatureFile", "-I", gvcf_path
            ], check=True)
        else:
            print(f"Index exists for {sample}")
        
        valid_gvcfs.append(gvcf_path)

    return valid_gvcfs


# === Combine GVCFs into a cohort file ===
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


# === Main ===
def main():
    config = load_config()

    reference = config["reference"]
    output_dir = config["output_dir"]
    gatk_path = config["gatk_path"]
    samples = config["samples"]  # list of sample IDs with _RG suffixes

    os.makedirs(output_dir, exist_ok=True)

    # Validate and collect GVCFs
    gvcfs = validate_and_index_gvcfs(samples, output_dir, gatk_path)

    # Output path for cohort GVCF
    combined_gvcf = os.path.join(output_dir, "cohort_combined.g.vcf")
    combine_gvcfs(reference, gvcfs, combined_gvcf, gatk_path)


if __name__ == "__main__":
    main()
