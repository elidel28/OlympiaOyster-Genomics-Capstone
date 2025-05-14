# Olympia Oyster GTSeq Panel Pipeline

## Project Description:

## About the Data:

## Tools Needed:

## Scipts Overview:

align_reads : This script aligns paired-end FASTQ reads to a reference genome using BWA-MEM and generates SAM files for each sample.
Inputs: Reference genome FASTA file (.fasta) and paired-end FASTQ files (R1 and R2 for each sample)
Outputs: SAM files for each sample saved in the results/ directory (<sample>.sam)

BQSRscripts : This script runs a GATK-based pipeline for base quality score recalibration (BQSR), variant calling, and filtering on multiple BAM files, producing final GVCF files for downstream analysis.

Inputs: 
- Sorted and read-group-added BAM files (.sorted.bam)
- Reference genome FASTA (.fasta)
- Known variant set for BQSR (generated from initial variant calls)
- Sample file names list (defined in files list)

Outputs: 
- Indexed BAM files (.bam.bai)
- Initial raw VCF (.raw.vcf)
- Filtered VCF (.filtered.vcf)
- Passing variants only VCF (.passonly.vcf)
- Recalibration table (.recal.table)
- Recalibrated BAM file (.recalibrated.bam)
- Final GVCF file (.final.g.vcf) — used for joint genotyping
