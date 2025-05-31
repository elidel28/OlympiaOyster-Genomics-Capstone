#!/bin/bash
set -e  # Exit immediately on any error

# === Clean BOM (Byte Order Mark) from samples.csv if it exists ===
echo "Cleaning BOM from samples.csv if present..."
sed -i '1s/^\xEF\xBB\xBF//' config/samples.csv

=== Activate environment ===
echo "Activating conda environment..."
conda activate olympia-pipeline

# === Check and Create Reference Indices ===
REFERENCE_FASTA=$(grep '^reference:' config/config.yaml | awk '{print $2}')
FAI_FILE="${REFERENCE_FASTA}.fai"
DICT_FILE="${REFERENCE_FASTA%.fasta}.dict"
BWA_INDEX_PREFIX="${REFERENCE_FASTA}"

echo "Checking reference genome indices..."

# GATK .fai
if [ ! -f "$FAI_FILE" ]; then
    echo "Creating FASTA index (.fai)..."
    samtools faidx "$REFERENCE_FASTA"
else
    echo "Found FASTA index: $FAI_FILE"
fi

# GATK .dict
if [ ! -f "$DICT_FILE" ]; then
    echo "Creating sequence dictionary (.dict)..."
    gatk CreateSequenceDictionary -R "$REFERENCE_FASTA"
else
    echo "Found sequence dictionary: $DICT_FILE"
fi

# BWA index files
BWA_EXTS=(.amb .ann .bwt .pac .sa)
MISSING_BWA_INDEX=false
for ext in "${BWA_EXTS[@]}"; do
    if [ ! -f "${BWA_INDEX_PREFIX}${ext}" ]; then
        MISSING_BWA_INDEX=true
        break
    fi
done

if [ "$MISSING_BWA_INDEX" = true ]; then
    echo "Creating BWA index..."
    bwa index "$REFERENCE_FASTA"
else
    echo "Found BWA index files."
fi

echo "Reference genome indexing check complete."

# === Run Pipeline ===
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
