"""
Single-sample variant calling, SNP filtering, and window-based SNP density filtering.

Steps:
1. Variant calling (GATK)
2. Hard filtering of variants
3. Base Quality Score Recalibration (BQSR)
4. Final GVCF calling
5. Genotyping (placeholder for single-sample)
6. SNP selection
7. Annotate with MAF and F_MISSING, filter
8. Convert to BED
9. Create 30bp windows
10. Count SNPs per window
11. Keep SNPs with no neighbors within 30bp
"""

import os
import subprocess
import yaml


def load_config(path="config/config.yaml"):
    with open(path, "r") as f:
        return yaml.safe_load(f)


def run_variant_workflow(sample_id, config):
    ref = config["reference"]
    out_dir = config["output_dir"]
    bam_dir = config["sorted_bam_dir"]
    genome_file = config["genome_file"]
    gatk = config["gatk_path"]
    threads = str(config["threads"])

    sorted_bam = os.path.join(bam_dir, f"{sample_id}.sorted.bam")
    raw_vcf = os.path.join(out_dir, f"{sample_id}.raw.vcf")
    filtered_vcf = os.path.join(out_dir, f"{sample_id}.filtered.vcf")
    passonly_vcf = os.path.join(out_dir, f"{sample_id}.passonly.vcf")
    recal_table = os.path.join(out_dir, f"{sample_id}.recal.table")
    recal_bam = os.path.join(out_dir, f"{sample_id}.recalibrated.bam")
    final_gvcf = os.path.join(out_dir, f"{sample_id}.final.g.vcf")
    joint_vcf = os.path.join(out_dir, f"{sample_id}.joint_genotyped.vcf")
    snps_vcf = os.path.join(out_dir, f"{sample_id}.joint_snps.vcf")
    annotated_vcf = os.path.join(out_dir, f"{sample_id}.annotated_snps.vcf")
    maf_filtered_vcf = os.path.join(out_dir, f"{sample_id}.selected_snps.vcf")
    snp_bed = os.path.join(out_dir, f"{sample_id}.snps.bed")
    snp_windows = os.path.join(out_dir, f"{sample_id}.snp_windows_30bp.bed")
    snp_counts = os.path.join(out_dir, f"{sample_id}.snps_with_counts.txt")
    final_snps_bed = os.path.join(out_dir, f"{sample_id}.final_snps.bed")

    os.makedirs(out_dir, exist_ok=True)

    subprocess.run([gatk, "HaplotypeCaller", "-R", ref, "-I", sorted_bam, "-O", raw_vcf, "--native-pair-hmm-threads", threads], check=True)

    subprocess.run([gatk, "VariantFiltration", "-V", raw_vcf,
                    "--filter-expression", "QD < 2.0 || FS > 60.0 || MQ < 40.0",
                    "--filter-name", "FAIL", "-O", filtered_vcf], check=True)

    subprocess.run([gatk, "SelectVariants", "-V", filtered_vcf,
                    "--exclude-filtered", "-O", passonly_vcf], check=True)

    subprocess.run([gatk, "BaseRecalibrator", "-I", sorted_bam, "-R", ref,
                    "--known-sites", passonly_vcf, "-O", recal_table,
                    "--num-threads", threads], check=True)

    subprocess.run([gatk, "ApplyBQSR", "-R", ref, "-I", sorted_bam,
                    "--bqsr-recal-file", recal_table, "-O", recal_bam,
                    "--num-threads", threads], check=True)

    subprocess.run([gatk, "HaplotypeCaller", "-R", ref, "-I", recal_bam,
                    "-O", final_gvcf, "-ERC", "GVCF", "--native-pair-hmm-threads", threads], check=True)

    subprocess.run([gatk, "GenotypeGVCFs", "-R", ref, "-V", final_gvcf, "-O", joint_vcf], check=True)

    subprocess.run([gatk, "SelectVariants", "-V", joint_vcf,
                    "--select-type-to-include", "SNP", "-O", snps_vcf], check=True)

    # Annotate with MAF and F_MISSING
    with open(annotated_vcf, 'w') as annot_out:
        fill_tags = subprocess.Popen([
            "bcftools", "+fill-tags", snps_vcf, "-Ou", "--", "-t", "MAF,F_MISSING"
        ], stdout=subprocess.PIPE)
        subprocess.run(["bcftools", "view", "-Ov"], stdin=fill_tags.stdout, stdout=annot_out, check=True)

    # Filter based on MAF > 0.1 and F_MISSING < 0.05
    with open(maf_filtered_vcf, 'w') as filt_out:
        subprocess.run([
            "bcftools", "view",
            "-i", "MAF>0.1 && F_MISSING<0.05",
            annotated_vcf
        ], stdout=filt_out, check=True)

    # Convert to BED
    with open(snp_bed, 'w') as bed_out:
        subprocess.run([
            "bcftools", "query",
            "-f", "%CHROM\t%POS0\t%POS\t%ID\n",
            maf_filtered_vcf
        ], stdout=bed_out, check=True)

    # Create 30bp windows
    with open(snp_windows, 'w') as win_out:
        subprocess.run([
            "bedtools", "slop", "-i", snp_bed, "-g", genome_file, "-b", "30"
        ], stdout=win_out, check=True)

    # Count SNPs per window
    with open(snp_counts, 'w') as count_out:
        subprocess.run([
            "bedtools", "intersect", "-a", snp_windows, "-b", snp_bed, "-c"
        ], stdout=count_out, check=True)

    # Filter SNPs with no neighbors within 30bp
    with open(snp_counts) as infile, open(final_snps_bed, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            if int(fields[4]) == 1:
                outfile.write('\t'.join(fields[:3]) + '\n')

    print(f"Finished processing sample: {sample_id}")


def main():
    config = load_config()
    for sample_id in config["samples"]:
        run_variant_workflow(sample_id, config)


if __name__ == "__main__":
    main()
