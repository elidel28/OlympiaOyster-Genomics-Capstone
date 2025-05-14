# Olympia Oyster GTSeq Panel Pipeline

## Project Description:

## About the Data:

## Tools Needed:

## Scipts Overview:

align_reads : This script aligns paired-end FASTQ reads to a reference genome using BWA-MEM and generates SAM files for each sample.

Inputs: Reference genome FASTA file (.fasta) and paired-end FASTQ files (R1 and R2 for each sample)

Outputs: SAM files for each sample saved in the results/ directory (<sample>.sam)
