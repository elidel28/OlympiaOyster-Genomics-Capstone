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

### Align_Reads: 
  This script aligns paired-end FASTQ reads to a reference genome using BWA-MEM and generates SAM files for each sample.

Inputs: 
- Reference genome FASTA file (.fasta)
- Paired-end FASTQ files (R1 and R2 for each sample)

Outputs:
- Sequence Alignment Map (SAM) files for each sample saved in the results/ directory (<sample>.sam)


### SAMtools_Script: 
  Converts SAM files to BAM (Binary Alignment Map) files, sorts the BAMs, and indexes them for further use in the pipeline.

Inputs:
- Sequence Alignment Map files (.sam)

Outputs:
- Binary Alignment Map files (.bam)
- Sorted BAM files (.sorted.bam)
- Sorted BAM Index files (.sorted.bam.bai)


### VariantCallingScripts: 
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
