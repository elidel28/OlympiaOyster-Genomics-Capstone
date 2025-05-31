"""
Variant Calling Pipeline using GATK.

Per sample, this script performs:
1. BAM indexing
2. Initial variant calling (HaplotypeCaller)
3. Variant filtration
4. Base Recalibration (BQSR)
5. Final GVCF generation
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

def index_bam(bam_path):
    subprocess.run(["samtools", "index", bam_path], check=True)


def haplotype_caller(gatk, ref, input_bam, output_vcf, threads, java_mem=None, emit_gvcf=False):
    cmd = [gatk]
    if java_mem:
        cmd += ["--java-options", f"-Xmx{java_mem}"]

    cmd += [
        "HaplotypeCaller",
        "-R", ref,
        "-I", input_bam,
        "-O", output_vcf,
        "--native-pair-hmm-threads", str(threads)
    ]
    if emit_gvcf:
        cmd += ["-ERC", "GVCF"]

    subprocess.run(cmd, check=True)


def hard_filter_variants(gatk, input_vcf, output_filtered, output_passonly):
    subprocess.run([
        gatk, "VariantFiltration",
        "-V", input_vcf,
        "--filter-expression", "QD < 2.0 || FS > 60.0 || MQ < 40.0",
        "--filter-name", "FAIL",
        "-O", output_filtered
    ], check=True)

    subprocess.run([
        gatk, "SelectVariants",
        "-V", output_filtered,
        "--exclude-filtered",
        "-O", output_passonly
    ], check=True)


def run_bqsr(gatk, ref, input_bam, known_sites, out_table):
    subprocess.run([
        gatk, "BaseRecalibrator",
        "-I", input_bam,
        "-R", ref,
        "--known-sites", known_sites,
        "-O", out_table
    ], check=True)


def apply_bqsr(gatk, ref, input_bam, recal_table, output_bam):
    subprocess.run([
        gatk, "ApplyBQSR",
        "-R", ref,
        "-I", input_bam,
        "--bqsr-recal-file", recal_table,
        "-O", output_bam
    ], check=True)


def process_sample(sample_id, config):
    ref = config["reference"]
    gatk = config["gatk_path"]
    threads = config["threads"]
    java_mem = config.get("gatk_mem", "8g")

    bam_dir = config["sorted_bam_dir"]
    variant_dir = config["variant_dir"]
    recal_dir = config["recal_dir"]
    gvcf_dir = config["gvcf_dir"]

    os.makedirs(variant_dir, exist_ok=True)
    os.makedirs(recal_dir, exist_ok=True)
    os.makedirs(gvcf_dir, exist_ok=True)

    sorted_bam = os.path.join(bam_dir, f"{sample_id}.sorted.bam")
    raw_vcf = os.path.join(variant_dir, f"{sample_id}.raw.vcf")
    filtered_vcf = os.path.join(variant_dir, f"{sample_id}.filtered.vcf")
    passonly_vcf = os.path.join(variant_dir, f"{sample_id}.passonly.vcf")
    recal_table = os.path.join(recal_dir, f"{sample_id}.recal.table")
    recal_bam = os.path.join(recal_dir, f"{sample_id}.recalibrated.bam")
    final_gvcf = os.path.join(gvcf_dir, f"{sample_id}.final.g.vcf")

    print(f"Processing sample: {sample_id}")
    index_bam(sorted_bam)
    haplotype_caller(gatk, ref, sorted_bam, raw_vcf, threads, java_mem=java_mem)
    hard_filter_variants(gatk, raw_vcf, filtered_vcf, passonly_vcf)
    run_bqsr(gatk, ref, sorted_bam, passonly_vcf, recal_table)
    apply_bqsr(gatk, ref, sorted_bam, recal_table, recal_bam)
    haplotype_caller(gatk, ref, recal_bam, final_gvcf, threads, java_mem=java_mem, emit_gvcf=True)


def main():
    config = load_config()
    samples = read_sample_sheet(config["sample_sheet"])
    for sample in samples:
        process_sample(sample["sample_id"], config)

if __name__ == "__main__":
    main()