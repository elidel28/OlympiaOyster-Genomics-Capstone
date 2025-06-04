"""
GTseq Primer Design Script Using Primer3

This script performs three core tasks:
1. Converts SNP-centered FASTA sequences to Primer3 input format
2. Runs Primer3 to design primers
3. Parses Primer3 output and writes a GTseq-compatible CSV file

Each FASTA entry should be ~100 bp, centered on a SNP of interest.
"""

import os
import csv
import subprocess

# === CONFIGURATION ===
INPUT_FASTA = "results/gtseq_design/cohort_snps_100bp.fa"
OUTPUT_DIR = "results/gtseq_design"
PRIMER3_INPUT = os.path.join(OUTPUT_DIR, "input.pr3")
PRIMER3_OUTPUT = os.path.join(OUTPUT_DIR, "output.txt")
GTSEQ_CSV = os.path.join(OUTPUT_DIR, "gtseq_panel.csv")


def fasta_to_primer3_input(fasta_path, output_path):
    """
    Convert SNP-centered FASTA entries to Primer3 input format.
    Each entry includes a SEQUENCE_TARGET at the SNP site.
    """
    def write_block(seq_id, sequence, out):
        center = len(sequence) // 2
        out.write(f"SEQUENCE_ID={seq_id}\n")
        out.write(f"SEQUENCE_TEMPLATE={sequence}\n")
        out.write(f"SEQUENCE_TARGET={center},1\n")
        out.write(f"SEQUENCE_INTERNAL_EXCLUDED_REGION={center},1\n")
        out.write("PRIMER_TASK=generic\n")
        out.write("PRIMER_PICK_LEFT_PRIMER=1\n")
        out.write("PRIMER_PICK_RIGHT_PRIMER=1\n")
        out.write("PRIMER_OPT_SIZE=20\n")
        out.write("PRIMER_MIN_SIZE=18\n")
        out.write("PRIMER_MAX_SIZE=25\n")
        out.write("PRIMER_OPT_TM=60.0\n")
        out.write("PRIMER_MIN_TM=57.0\n")
        out.write("PRIMER_MAX_TM=63.0\n")
        out.write("PRIMER_MIN_GC=20.0\n")
        out.write("PRIMER_MAX_GC=80.0\n")
        out.write("PRIMER_PRODUCT_SIZE_RANGE=80-120\n")
        out.write("PRIMER_NUM_RETURN=1\n")
        out.write("PRIMER_EXPLAIN_FLAG=1\n")
        out.write("=\n")

    with open(fasta_path) as fasta, open(output_path, "w") as out:
        seq_id = None
        seq = []
        for line in fasta:
            line = line.strip()
            if line.startswith(">"):
                if seq_id:
                    write_block(seq_id, "".join(seq), out)
                seq_id = line[1:]  # drop ">"
                seq = []
            else:
                seq.append(line)
        if seq_id:
            write_block(seq_id, "".join(seq), out)


def run_primer3(input_file, output_file):
    """
    Runs Primer3 on the provided input file and writes raw output.
    """
    try:
        with open(output_file, "w") as out:
            subprocess.run(["primer3_core", input_file], stdout=out, check=True)
    except subprocess.CalledProcessError as e:


def parse_primer3_output(input_file, csv_output):
    """
    Parses Primer3 output and writes a CSV with SNP ID, primers, and amplicon size.
    """
    with open(input_file) as infile, open(csv_output, "w", newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(["SNP_ID", "Forward_Primer", "Reverse_Primer", "Amplicon_Size", "Notes"])

        block = {}
        for line in infile:
            line = line.strip()
            if not line:
                continue
            if line == "=":
                # Process a completed block
                if all(k in block for k in ["SEQUENCE_ID", "PRIMER_LEFT_0", "PRIMER_RIGHT_0",
                                            "PRIMER_LEFT_0_SEQUENCE", "PRIMER_RIGHT_0_SEQUENCE"]):
                    left_pos = int(block["PRIMER_LEFT_0"].split(",")[0])
                    right_pos = int(block["PRIMER_RIGHT_0"].split(",")[0])
                    amplicon_size = right_pos + 1 - left_pos
                    writer.writerow([
                        block["SEQUENCE_ID"],
                        block["PRIMER_LEFT_0_SEQUENCE"],
                        block["PRIMER_RIGHT_0_SEQUENCE"],
                        amplicon_size,
                        ""
                    ])
                block = {}
            elif "=" in line:
                k, v = line.split("=", 1)
                block[k] = v


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    fasta_to_primer3_input(INPUT_FASTA, PRIMER3_INPUT)
    run_primer3(PRIMER3_INPUT, PRIMER3_OUTPUT)
    parse_primer3_output(PRIMER3_OUTPUT, GTSEQ_CSV)


if __name__ == "__main__":
    main()
