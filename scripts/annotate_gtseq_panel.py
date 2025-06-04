"""
annotate_gtseq_panel.py

This script enriches a GT-seq panel CSV (from Primer3) by adding genomic annotations:
- Chromosome and position (from metadata)
- Start and end coordinates (from BED file)
- Minor allele frequency (MAF)

Inputs:
- gtseq_panel.csv: Primer3 output with SNP_ID, primer sequences, and amplicon size
- snp_metadata.tsv: SNP metadata (CHROM, POS, REF, ALT, MAF)
- cohort_final_snps_100bp.bed: BED file with CHROM, START, END

Output:
- gtseq_panel_annotated.csv: Final annotated panel ready for synthesis or analysis
"""

import csv

PANEL_CSV = "/mnt/jupiter/johnsonlab/Capstone_proj/results/gtseq_design/gtseq_panel.csv"
METADATA_TSV = "/mnt/jupiter/johnsonlab/Capstone_proj/results/gtseq_design/snp_metadata.tsv"
BED_FILE = "/mnt/jupiter/johnsonlab/Capstone_proj/results/gtseq_design/cohort_final_snps_100bp.bed"
OUTPUT_CSV = "/mnt/jupiter/johnsonlab/Capstone_proj/results/gtseq_design/gtseq_panel_annotated.csv"


def load_gtseq_panel(filepath):
    """Loads GTseq panel from CSV file."""
    with open(filepath, newline='') as f:
        return list(csv.DictReader(f))


def load_snp_metadata(filepath):
    """Loads SNP metadata from TSV (CHROM, POS, REF, ALT, MAF)."""
    metadata = []
    with open(filepath, newline='') as f:
        for line in f:
            chrom, pos, ref, alt, maf = line.strip().split('\t')
            metadata.append({
                "CHROM": chrom,
                "POS": pos,
                "REF": ref,
                "ALT": alt,
                "MAF": maf
            })
    return metadata


def load_bed_coordinates(filepath):
    """Loads BED file (CHROM, START, END)."""
    bed = []
    with open(filepath, newline='') as f:
        for line in f:
            chrom, start, end = line.strip().split('\t')
            bed.append({
                "BED_CHROM": chrom,
                "START": start,
                "END": end
            })
    return bed


def combine_annotations(panel, metadata, bed):
    """Combines GTseq panel, SNP metadata, and BED info."""
    combined = []
    for i in range(len(panel)):
        combined_row = {
            "SNP_ID": panel[i]["SNP_ID"],
            "CHROM": metadata[i]["CHROM"],
            "POS": metadata[i]["POS"],
            "START": bed[i]["START"],
            "END": bed[i]["END"],
            "MAF": metadata[i]["MAF"],
            "Forward_Primer": panel[i]["Forward_Primer"],
            "Reverse_Primer": panel[i]["Reverse_Primer"],
            "Amplicon_Size": panel[i]["Amplicon_Size"]
        }
        combined.append(combined_row)
    return combined


def write_annotated_panel(rows, output_file):
    """Writes annotated panel to CSV."""
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main():
    panel = load_gtseq_panel(PANEL_CSV)
    metadata = load_snp_metadata(METADATA_TSV)
    bed = load_bed_coordinates(BED_FILE)
    final_rows = combine_annotations(panel, metadata, bed)
    write_annotated_panel(final_rows, OUTPUT_CSV)


if __name__ == "__main__":
    main()
