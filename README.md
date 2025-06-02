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

---

## Repository Structure

```
OlympiaOyster-Genomics-Capstone/
├── config/                 # Pipeline configuration and sample metadata
│   ├── config.yaml
│   └── samples.csv
│
├── results/                # All pipeline outputs organized by step
│   ├── sam/                # SAM files from BWA
│   ├── bam/                # BAM files
│   ├── bam/sorted/         # Sorted and indexed BAMs
│   ├── variants/           # Raw + filtered VCFs
│   ├── recalibration/      # BQSR tables and recalibrated BAMs
│   ├── gvcfs/              # Per-sample GVCFs
│   ├── combined/           # Combined GVCF for joint calling
│   ├── genotyped/          # Joint-genotyped VCFs and SNPs
│   └── plink_analysis/     # Filtered SNPs and PLINK outputs
│
├── scripts/                # Individual pipeline steps
│   ├── align_reads_with_RG.py
│   ├── sam_to_sorted_bam.py
│   ├── variant_calling_bqsr.py
│   ├── combine_gvcfs.py
│   ├── joint_genotyping_select_snps.py
│   ├── plink_filtering_pipeline.py
│   └── snp_annotation_and_filtering.py
│
├── environment.yml         # Conda environment file
├── run_pipeline.sh         # One-command pipeline runner
└── README.md             
```

---

## Configuration

### `config.yaml`

* Defines paths for input/output directories, reference genome, and tools.
* Set `reference` and `fastq_dir` to where your files are stored.

### `samples.csv`

* List your sample IDs and FASTQ filenames:

```csv
sample_id,forward_read,reverse_read
10KP,10KP_R1.fastq.gz,10KP_R2.fastq.gz
14KP,14KP_R1.fastq.gz,14KP_R2.fastq.gz
```

---

## Getting Started

### Step 1: Clone the Repository

```bash
git clone https://github.com/elidel28/OlympiaOyster-Genomics-Capstone.git
cd OlympiaOyster-Genomics-Capstone
```

### Step 2: Create the Conda Environment

```bash
conda env create -f environment.yml
conda activate olympia-pipeline
```

### Step 3: Set Up Data Paths

* Create a symlink to your data folder:

```bash
ln -s /mnt/jupiter/johnsonlab johnsonlab_data
```

* Update `config.yaml`:

```yaml
reference: johnsonlab_data/Final_assembly_Olympia/your_reference.fasta
fastq_dir: johnsonlab_data/lcwg_Olympia/
```

### Step 4: Configure Samples

Edit `config/samples.csv` to list all samples and their FASTQ files names.

---

## Run the Pipeline

### Preferred:

```bash
bash run_pipeline.sh
```

### Manual (step-by-step):

```bash
# Clean BOM if needed
sed -i '1s/^\xEF\xBB\xBF//' config/samples.csv

# Generate index files required for BWA MEM
bwa index your_reference.fasta

# Run each step
python scripts/align_reads_with_RG.py
python scripts/sam_to_sorted_bam.py
python scripts/variant_calling_bqsr.py
python scripts/combine_gvcfs.py
python scripts/joint_genotyping_select_snps.py
python scripts/plink_filtering_pipeline.py
```

## Scripts Overview:

| Script                            | Reads From                    | Writes To                       |
| --------------------------------- | ----------------------------- | ------------------------------- |
| `align_reads_with_RG.py`          | `fastq_dir`, `samples.csv`    | `sam_dir`                       |
| `sam_to_sorted_bam.py`            | `sam_dir`                     | `bam_dir`, `sorted_bam_dir`     |
| `variant_calling_bqsr.py`         | `sorted_bam_dir`, `reference` | `variants/`, `recal/`, `gvcfs/` |
| `combine_gvcfs.py`                | `gvcfs/`, `reference`         | `combined/`                     |
| `joint_genotyping_select_snps.py` | `combined/`, `reference`      | `genotyped/`                    |
| `plink_filtering_pipeline.py`     | `genotyped/`                  | `plink_analysis/`               |

### align_reads_with_RG: 
This script aligns paired-end FASTQ reads to a reference genome using BWA-MEM, automatically extracts read group (@RG) information from each R1 file, and generates SAM files for each sample.

Inputs:
- Reference genome FASTA file (.fasta)
- Paired-end FASTQ files (R1 and R2) for each sample (gzipped)
- Output directory path
- Path to bwa executable
- Genome-specific thread count (user-defined)

Outputs:
- SAM files for each sample saved in the output directory


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
- Sample file names list

Outputs: 
- Indexed BAM files (.bam.bai) (only indexed if indexes are not present from SAMtools_script)
- Initial raw VCF (.raw.vcf)
- Filtered VCF (.filtered.vcf)
- Passing variants only VCF (.passonly.vcf)
- Recalibration table (.recal.table)
- Recalibrated BAM file (.recalibrated.bam)
- Final GVCF file (.final.g.vcf) — used for joint genotyping

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
- BED file with genomic coordinates of final SNPs

Step-by-Step Overview:

#### Step 1: Convert VCF to PLINK Format (Translates VCF into PLINK binary format for efficient processing)

Input: cohort_snps.vcf
Output: cohort_plink.bed/.bim/.fam

Step 2: Compute SNP Statistics (Calculates MAF and per-site/per-sample missingness)

Input: PLINK files from step 1
Output: cohort_stats.frq, .imiss, .lmiss

Step 3: Filter SNPs (Removes SNPs with MAF < 0.1 or missingness > 5%)

Input: PLINK files from step 1
Output: cohort_filtered.bed/.bim/.fam


Step 4: LD Pruning - Identify/removes SNPs in strong linkage disequilibrium (R^2 > 0.2)

Input: Filtered PLINK dataset
Output: cohort_prune.prune.in


Step 5: Renives duplicate entries frm the LD pune list for error handling

Input: cohort_prune.prune.in
Output: cohort_prune.prune.in.dedup


Step 6: SNP selection and exports to VCF

Input: Deduplicated SNP list and filtered PLINK dataset
Output: cohort_filtered_snps.vcf


Step 7: Generate BED file of SNP coordinates

Input: cohort_filtered.bim
Output: cohort_filtered.snps.bed


















