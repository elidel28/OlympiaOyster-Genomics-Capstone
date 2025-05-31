"""
SAM to Sorted and Indexed BAM Pipeline using SAMtools.

Steps:
1. Convert SAM to BAM
2. Sort BAM by coordinates
3. Index sorted BAM
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


def sam_to_bam(sam_path, bam_path):
    subprocess.run([
        "samtools", "view", "-bS", sam_path, "-o", bam_path
    ], check=True)


def sort_bam(bam_path, sorted_bam_path, threads, ram):
    subprocess.run([
        "samtools", "sort",
        "-@", str(threads),
        "-m", ram,
        "-o", sorted_bam_path,
        bam_path
    ], check=True)


def index_bam(sorted_bam_path, threads):
    subprocess.run([
        "samtools", "index",
        "-@", str(threads),
        sorted_bam_path
    ], check=True)


def process_sample(sample_id, config):
    sam_path = os.path.join(config["sam_dir"], f"{sample_id}.sam")
    bam_path = os.path.join(config["bam_dir"], f"{sample_id}.bam")
    sorted_bam_path = os.path.join(config["sorted_bam_dir"], f"{sample_id}.sorted.bam")

    os.makedirs(config["bam_dir"], exist_ok=True)
    os.makedirs(config["sorted_bam_dir"], exist_ok=True)

    print(f"Processing sample: {sample_id}")
    sam_to_bam(sam_path, bam_path)
    sort_bam(bam_path, sorted_bam_path, config["threads"], config["samtools_ram"])
    index_bam(sorted_bam_path, config["threads"])


def main():
    config = load_config()
    samples = read_sample_sheet(config["sample_sheet"])
    for sample in samples:
        process_sample(sample["sample_id"], config)


if __name__ == "__main__":
    main()