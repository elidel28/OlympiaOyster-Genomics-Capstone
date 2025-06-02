import subprocess
import os
from collections import Counter

def run_plink_filtering_pipeline(
    vcf_input,
    work_dir,
    maf_threshold=0.1,
    geno_threshold=0.05,
    ld_window=50,
    ld_step=5,
    ld_threshold=0.2,
):
    os.makedirs(work_dir, exist_ok=True)

    # Define file paths
    plink_prefix = f"{work_dir}/cohort_plink"
    stats_prefix = f"{work_dir}/cohort_stats"
    filtered_prefix = f"{work_dir}/cohort_filtered"
    prune_prefix = f"{work_dir}/cohort_prune"
    final_vcf = f"{work_dir}/cohort_filtered_snps.vcf"
    snp_bed = f"{work_dir}/cohort_filtered.snps.bed"
    plink_bim = f"{filtered_prefix}.bim"
    prune_in_file = f"{prune_prefix}.prune.in"
    deduped_extract = f"{prune_prefix}.prune.in.dedup"

    print("Converting VCF to PLINK binary format...")
    subprocess.run([
        "plink", "--vcf", vcf_input,
        "--allow-extra-chr",
        "--make-bed", "--out", plink_prefix
    ], check=True)

    print("Calculating MAF and missingness...")
    subprocess.run([
        "plink", "--bfile", plink_prefix,
        "--allow-extra-chr",
        "--freq", "--missing",
        "--out", stats_prefix
    ], check=True)

    print("Filtering SNPs with MAF > 0.1 and missingness < 0.05...")
    subprocess.run([
        "plink", "--bfile", plink_prefix,
        "--allow-extra-chr",
        "--maf", str(maf_threshold),
        "--geno", str(geno_threshold),
        "--make-bed",
        "--out", filtered_prefix
    ], check=True)

    print("Performing LD pruning...")
    subprocess.run([
        "plink", "--bfile", filtered_prefix,
        "--allow-extra-chr",
        "--indep-pairwise", str(ld_window), str(ld_step), str(ld_threshold),
        "--out", prune_prefix
    ], check=True)

    print("Deduplicating SNP IDs in prune list...")
    seen = set()
    with open(prune_in_file, "r") as infile, open(deduped_extract, "w") as outfile:
        for line in infile:
            snp = line.strip()
            if snp not in seen:
                seen.add(snp)
                outfile.write(snp + "\n")

    with open(prune_in_file) as f:
        snps = [line.strip() for line in f]
        dupes = [snp for snp, count in Counter(snps).items() if count > 1]
        print(f"Found {len(dupes)} duplicated SNP IDs (removed).")

    print("Extracting pruned SNPs and writing final VCF...")
    subprocess.run([
        "plink", "--bfile", filtered_prefix,
        "--allow-extra-chr",
        "--extract", deduped_extract,
        "--recode", "vcf",
        "--out", final_vcf.replace(".vcf", "")
    ], check=True)

    print("Generating BED file of final SNPs...")
    with open(plink_bim) as bim_file, open(snp_bed, 'w') as bed_file:
        for line in bim_file:
            fields = line.strip().split()
            chrom = fields[0]
            pos = int(fields[3])
            bed_file.write(f"{chrom}\t{pos - 1}\t{pos}\n")

    print("PLINK SNP filtering and LD pruning complete.")
