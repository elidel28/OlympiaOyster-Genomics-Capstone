# Olympia Oyster GTSeq Panel Pipeline

## Project Description:
The project aims to develop a genomic toolset and analysis pipeline to assess the genetic relatedness between farmed and wild Olympia oysters (*Ostrea lurida*) using a categorized panel of location-specific single-nucleotide polymorphisms (SNPs) identified through GTseq. This project supports government efforts, such as those by the California Department of Fish and Wildlife, in evaluating the suitability of farmed oysters for reintroduction by comparing their genetic diversity to that of various unique wild populations along the West Coast.

## About the Data:
Data used for our pipeline consisted of raw genomic data from a larger dataset of 254 individual Olympia oysters across 14 locations on the West Coast.

## Tools Needed:
- Burrows-Wheeler Aligner (BWA): Align sequencing reads from individual oyster samples to the Olympia oyster reference genome
- SamTools: Interacting with SAM/BAM files
- *GATK*: HaplotypeCaller: Call germline SNPs and indels via local re-assembly of haplotypes
- *GATK*: VariantFiltration: Hard filter variant calls based on given criteria and label failed calls
- *GATK*: SelectVariants: Generate a vcf file with all the failed calls removed
- *GATK*: BaseRecalibrator: Generates recalibration table for Base Quality Score Recalibration (BQSR)
- *GATK*: ApplyBQSR: applies recalibration table to adjust the quality scores of sequencing reads
- *GATK*: CombineGVCFs: Merge individual gVCFs to single multi-sample gVFC
- *GATK*: GenotypeGVCFs: Joint genotyping on combined gVCF to produce a fully genotyped multi-sample VCF
- Plink: SNP-level filtering, linkage disequilibrium pruning, and population based quality control
- bcftools: VCF file operations (annotating, filtering, format conversion)
- bedtools: Genomiv interval operations to great targeted genomics regions and count SNPs per region

## Scripts Overview:

### align_reads_with_RG: 
This script aligns paired-end FASTQ reads to a reference genome using BWA-MEM, automatically extracts read group (@RG) information from each R1 file, and generates SAM files for each sample.

Inputs:
- Reference genome FASTA file (.fasta)
- Paired-end FASTQ files (R1 and R2) for each sample (gzipped)
- Output directory path
- Path to bwa executable
- Genome-specific thread count (user-defined)

Outputs:
- SAM files for each sample saved in the output directory (<sample>_test.sam)


### sam_to_sorted_bam: 
  Converts SAM files to BAM (Binary Alignment Map) files, sorts the BAMs, and indexes them for further use in the pipeline.

Inputs:
- Sequence Alignment Map files (.sam)

Outputs:
- Binary Alignment Map files (.bam)
- Sorted BAM files (.sorted.bam)
- Sorted BAM Index files (.sorted.bam.bai)


### variant_calling_bqsr: 
  This script runs a GATK-based pipeline for initial variant calling, filtering unwanted reads, base quality score recalibration (BQSR), and producing final GVCF files for downstream analysis.

Inputs: 
- Sorted and read-group-added BAM files (.sorted.bam)
- Reference genome FASTA (.fasta)
- Known variant set for BQSR (generated from initial variant calls)
- Sample file names list (defined in files list)

Outputs: 
- Indexed BAM files (.bam.bai) (only indexed if indexes are not present from SAMtools_script)
- Initial raw VCF (.raw.vcf)
- Filtered VCF (.filtered.vcf)
- Passing variants only VCF (.passonly.vcf)
- Recalibration table (.recal.table)
- Recalibrated BAM file (.recalibrated.bam)
- Final GVCF file (.final.g.vcf) — used for joint genotyping

Step-by-Step Overview:

Step 1: Index BAM
Input: <sample>.sorted.bam
Output: <sample>.sorted.bam.bai

Step 2: Initial Variant Calling
Input: <sample>.sorted.bam, reference
Output: <sample>.raw.vcf

Step 3: Hard Filtering of Variants
Input: <sample>.raw.vcf
Output: <sample>.filtered.vcf

Step 4: Select Passing Variants
Input: <sample>.filtered.vcf
Output: <sample>.passonly.vcf

Step 5: Base Recalibration (BQSR Table)
Input: <sample>.sorted.bam, <sample>.passonly.vcf, reference
Output: <sample>.recal.table

Step 6: Apply BQSR
Input: <sample>.sorted.bam, recal table
Output: <sample>.recalibrated.bam

Step 7: Final Variant Calling (GVCF mode)
Input: <sample>.recalibrated.bam, reference
Output: <sample>.final.g.vcf


### combine_gvcfs:
Checks for, and indexes per-sample GVCF files generated after BQSR and final variant calling. Once all GVCFs are confirmed and indexed, it combines them into a single multi-sample GVCF using GATK’s CombineGVCFs, in preparation for joint genotyping.

Inputs:
- Final per-sample GVCF files (<sample>.final.g.vcf)
- Reference genome FASTA file (.fasta)

Outputs:
- Combined multi-sample GVCF (cohort_combined.g.vcf)
- Indexed .idx files for each input GVCF

### joint_genotyping_select_snps:
Performs joint genotyping on a combined multi-sample GVCF using GATK’s GenotypeGVCFs, then extracts only SNPs from the resulting VCF using SelectVariants. It is a core step in transitioning from raw variant discovery to marker selection.

Inputs:
- Combined GVCF file (cohort_combined.g.vcf)
- Reference genome FASTA file (.fasta)

Outputs:
- Joint-genotyped VCF (cohort_genotyped.vcf)
- SNP-only VCF (cohort_snps.vcf)

### plink_filtering_pipeline:
Performs SNP-level filtering using PLINK. It converts a SNP-only VCF into PLINK format, calculates minor allele frequency (MAF) and missingness statistics, filters variants based on quality thresholds, and outputs a high-confidence VCF file for downstream analysis or GT-seq panel selection.

Inputs:
- SNP-only VCF file (cohort_snps.vcf)

Outputs:
- PLINK binary files (.bed, .bim, .fam)
- Summary statistics (.frq, .imiss, .lmiss)
- Filtered PLINK binary dataset
- Final filtered VCF (cohort_filtered_snps.vcf)

Step-by-Step Overview:

Step 1: Convert VCF to PLINK Format (Translates VCF into PLINK binary format for efficient processing)

Input: cohort_snps.vcf
Output: cohort_plink.bed/.bim/.fam

Step 2: Compute SNP Statistics (Calculates MAF and per-site/per-sample missingness)

Input: PLINK files from step 1
Output: cohort_stats.frq, .imiss, .lmiss

Step 3: Filter SNPs (Removes SNPs with MAF < 0.1 or missingness > 5%)

Input: PLINK files from step 1
Output: cohort_filtered.bed/.bim/.fam


Step 4: Export Filtered SNPs to VCF (Converts final filtered dataset back to VCF format)

Input: Filtered PLINK files
Output: cohort_filtered_snps.vcf

