#!/bin/bash
set -e  # The pipeline will exit if any step fails.
# echo "Activating conda environment..."
# conda activate olympia-pipeline

echo "========== STEP 1: Align Reads with BWA =========="
python scripts/align_reads_with_RG.py

echo "========== STEP 2: Convert SAM to Sorted and Indexed BAM =========="
python scripts/sam_to_sorted_bam.py

echo "========== STEP 3: Run Variant Calling and BQSR =========="
python scripts/variant_calling_bqsr.py

echo "========== STEP 4: Combine GVCFs =========="
python scripts/combine_gvcfs.py

echo "========== STEP 5: Joint Genotyping and SNP Selection =========="
python scripts/joint_genotyping_select_snps.py

echo "========== STEP 6: PLINK SNP Filtering =========="
python scripts/plink_filtering_pipeline.py

echo "========== PIPELINE COMPLETE =========="
