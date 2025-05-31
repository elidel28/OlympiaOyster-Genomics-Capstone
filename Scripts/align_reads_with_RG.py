"""
BWA alignment script with Read Group (RG) extraction.

Aligns paired-end FASTQ files to a reference genome using BWA MEM, while assigning read groups
based on FASTQ headers. Outputs SAM files per sample.
"""

import os
import csv
import gzip
import subprocess
import yaml


def load_config(config_path="config/config.yaml"):
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def read_sample_sheet(sample_sheet_path):
    """Read samples.csv and return list of sample dictionaries."""
    samples = []
    with open(sample_sheet_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            samples.append({
                "sample_id": row["sample_id"],
                "forward_read": row["forward_read"],
                "reverse_read": row["reverse_read"]
            })
    return samples


def extract_readgroup_info(r1_path, sample_id):
    """Extract @RG string from the first header of a gzipped FASTQ file."""
    with gzip.open(r1_path, "rt") as f:
        header = f.readline().strip()
    fields = header.split(":")
    instrument = fields[0].replace("@", "")
    flowcell = fields[2]
    lane = fields[3]
    pu = f"{flowcell}.{lane}"

    rg_dict = {
        'ID': instrument,
        'SM': sample_id,
        'PL': 'ILLUMINA',
        'LB': 'LCWG_Olympia',
        'PU': pu
    }

    return "@RG\\t" + "\\t".join([f"{k}:{v}" for k, v in rg_dict.items()])


def align_sample(sample, config):
    """Run BWA alignment for a sample."""
    r1_path = os.path.join(config["fastq_dir"], sample["forward_read"])
    r2_path = os.path.join(config["fastq_dir"], sample["reverse_read"])
    sample_id = sample["sample_id"]
    output_dir = config["sam_dir"]
    ref = config["reference"]
    threads = str(config["threads"])
    bwa_path = config["bwa_path"]

    os.makedirs(output_dir, exist_ok=True)
    output_sam = os.path.join(output_dir, f"{sample_id}.sam")

    print(f"Aligning {sample_id} → {output_sam}")
    rg_string = extract_readgroup_info(r1_path, sample_id)

    cmd = [
        bwa_path, "mem",
        "-t", threads,
        "-R", rg_string,
        ref, r1_path, r2_path
    ]

    with open(output_sam, "w") as out_f:
        subprocess.run(cmd, stdout=out_f)


def main():
    config = load_config()
    samples = read_sample_sheet(config["sample_sheet"])

    for sample in samples:
        if sample.get("forward_read") and sample.get("reverse_read"):
            align_sample(sample, config)
        else:
            print(f"Incomplete FASTQ pair for {sample['sample_id']} — skipping.")


if __name__ == "__main__":
    main()
