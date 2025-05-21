import os
import subprocess

# === Configuration ===

# Path to the reference genome FASTA file
REFERENCE_FASTA = "/path/to/reference.fasta"

# Directory for output files
OUTPUT_DIR = "/path/to/output"

# Directory containing sorted BAM files with read groups
BAM_DIR = "/path/to/bam_directory"

# Path to genome file (required for bedtools slop)
GENOME_FILE = "/path/to/genome.txt"

# Path to GATK executable
GATK_PATH = "/path/to/gatk"

# Number of threads to use for parallel operations
THREADS = "16"

# List of sample base filenames (without file extensions)
SAMPLES = [
    "sample1",
    # Add more sample names as needed
]

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# === Pipeline Execution ===

for sample in SAMPLES:
    print(f"🔄 Processing sample: {sample}")

    # Define paths for input and output files
    sorted_bam = f"{BAM_DIR}/{sample}.sorted.bam"
    raw_vcf = f"{OUTPUT_DIR}/{sample}.raw.vcf"
    filtered_vcf = f"{OUTPUT_DIR}/{sample}.filtered.vcf"
    passonly_vcf = f"{OUTPUT_DIR}/{sample}.passonly.vcf"
    recal_table = f"{OUTPUT_DIR}/{sample}.recal.table"
    recal_bam = f"{OUTPUT_DIR}/{sample}.recalibrated.bam"
    final_gvcf = f"{OUTPUT_DIR}/{sample}.final.g.vcf"
    joint_vcf = f"{OUTPUT_DIR}/{sample}.joint_genotyped.vcf"
    snps_vcf = f"{OUTPUT_DIR}/{sample}.joint_snps.vcf"
    annotated_vcf = f"{OUTPUT_DIR}/{sample}.annotated_snps.vcf"
    maf_filtered_vcf = f"{OUTPUT_DIR}/{sample}.selected_snps.vcf"
    snp_bed = f"{OUTPUT_DIR}/{sample}.snps.bed"
    snp_windows = f"{OUTPUT_DIR}/{sample}.snp_windows_30bp.bed"
    snp_counts = f"{OUTPUT_DIR}/{sample}.snps_with_counts.txt"
    final_snps_bed = f"{OUTPUT_DIR}/{sample}.final_snps.bed"

    # Step 1: Variant Calling
    subprocess.run([
        GATK_PATH, "HaplotypeCaller",
        "-R", REFERENCE_FASTA,
        "-I", sorted_bam,
        "-O", raw_vcf,
        "--native-pair-hmm-threads", THREADS
    ])

    # Step 2: Variant Filtering
    subprocess.run([
        GATK_PATH, "VariantFiltration",
        "-V", raw_vcf,
        "--filter-expression", "QD < 2.0 || FS > 60.0 || MQ < 40.0",
        "--filter-name", "FAIL",
        "-O", filtered_vcf
    ])
    subprocess.run([
        GATK_PATH, "SelectVariants",
        "-V", filtered_vcf,
        "--exclude-filtered",
        "-O", passonly_vcf
    ])

    # Step 3: Base Quality Score Recalibration (BQSR)
    subprocess.run([
        GATK_PATH, "BaseRecalibrator",
        "-I", sorted_bam,
        "-R", REFERENCE_FASTA,
        "--known-sites", passonly_vcf,
        "-O", recal_table,
        "--num-threads", THREADS
    ])
    subprocess.run([
        GATK_PATH, "ApplyBQSR",
        "-R", REFERENCE_FASTA,
        "-I", sorted_bam,
        "--bqsr-recal-file", recal_table,
        "-O", recal_bam,
        "--num-threads", THREADS
    ])

    # Step 4: Final Variant Calling
    subprocess.run([
        GATK_PATH, "HaplotypeCaller",
        "-R", REFERENCE_FASTA,
        "-I", recal_bam,
        "-O", final_gvcf,
        "-ERC", "GVCF",
        "--native-pair-hmm-threads", THREADS
    ])

    # Step 5: Joint Genotyping (single-sample case)
    subprocess.run([
        GATK_PATH, "GenotypeGVCFs",
        "-R", REFERENCE_FASTA,
        "-V", final_gvcf,
        "-O", joint_vcf
    ])

    # Step 6: Extract SNPs
    subprocess.run([
        GATK_PATH, "SelectVariants",
        "-V", joint_vcf,
        "--select-type-to-include", "SNP",
        "-O", snps_vcf
    ])

    # Step 7: Annotate and filter SNPs based on MAF and missingness
    with open(annotated_vcf, 'w') as annot_out:
        fill_tags = subprocess.Popen([
            "bcftools", "+fill-tags", snps_vcf, "-Ou", "--", "-t", "MAF,F_MISSING"
        ], stdout=subprocess.PIPE)
        subprocess.run([
            "bcftools", "view", "-Ov"
        ], stdin=fill_tags.stdout, stdout=annot_out)

    with open(maf_filtered_vcf, 'w') as filt_out:
        subprocess.run([
            "bcftools", "view",
            "-i", "MAF>0.1 && F_MISSING<0.05",
            annotated_vcf
        ], stdout=filt_out)

    # Step 8: Convert VCF to BED format
    with open(snp_bed, 'w') as bed_out:
        subprocess.run([
            "bcftools", "query",
            "-f", "%CHROM\t%POS0\t%POS\t%ID\n",
            maf_filtered_vcf
        ], stdout=bed_out)

    # Step 9: Create 30 bp windows around SNPs
    with open(snp_windows, 'w') as out_win:
        subprocess.run([
            "bedtools", "slop",
            "-i", snp_bed,
            "-g", GENOME_FILE,
            "-b", "30"
        ], stdout=out_win)

    # Step 10: Count SNPs per window
    with open(snp_counts, 'w') as out_count:
        subprocess.run([
            "bedtools", "intersect",
            "-a", snp_windows,
            "-b", snp_bed,
            "-c"
        ], stdout=out_count)

    # Step 11: Filter SNPs without nearby neighbors
    with open(snp_counts) as infile, open(final_snps_bed, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            if int(fields[4]) == 1:
                outfile.write('\t'.join(fields[:3]) + '\n')

    print(f"Finished processing: {sample}")
