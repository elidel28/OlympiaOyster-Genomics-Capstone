"""
GTseq Primer Design Pipeline using Primer3

This script performs:
1. Converts SNP-flanking FASTA sequences to Primer3 input format
2. Runs Primer3 to design primers
3. Parses Primer3 output and writes a GTseq-compatible CSV file
"""

import os
import subprocess
import csv

def fasta_to_primer3_input(fasta_path, output_path):

    def write_primer3_block(seq_id, sequence, out):
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
        out.write("PRIMER_EXPLAIN_FLAG=1\n")
        out.write("PRIMER_PRODUCT_SIZE_RANGE=80-120\n")
        out.write("PRIMER_NUM_RETURN=1\n")
        out.write("=\n")

    with open(fasta_path) as fasta, open(output_path, "w") as out:
        seq_id = None
        seq = []
        for line in fasta:
            if line.startswith(">"):
                if seq_id:
                    write_primer3_block(seq_id, "".join(seq), out)
                seq_id = line.strip()[1:]
                seq = []
            else:
                seq.append(line.strip())
        if seq_id:
            write_primer3_block(seq_id, "".join(seq), out)

def run_primer3(primer3_input_file, primer3_output_file):
    subprocess.run(["primer3_core", primer3_input_file], stdout=open(primer3_output_file, "w"), check=True)

def parse_primer3_output(primer3_output, csv_output):
    with open(primer3_output) as infile, open(csv_output, "w", newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(["SNP_ID", "Forward_Primer", "Reverse_Primer", "Amplicon_Size", "Notes"])
        block = {}
        count = 0
        for line in infile:
            line = line.strip()
            if not line:
                continue
            if line == "=":
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
                    count += 1
                block = {}
            else:
                if "=" in line:
                    k, v = line.split("=", 1)
                    block[k] = v

def main():
    input_fasta = "results/gtseq_design/cohort_snps_100bp.fa"
    output_dir = "results/gtseq_design"
    os.makedirs(output_dir, exist_ok=True)

    primer3_input_file = os.path.join(output_dir, "input.pr3")
    primer3_output_file = os.path.join(output_dir, "output.txt")
    gtseq_csv = os.path.join(output_dir, "gtseq_panel.csv")

    fasta_to_primer3_input(input_fasta, primer3_input_file)
    run_primer3(primer3_input_file, primer3_output_file)
    parse_primer3_output(primer3_output_file, gtseq_csv)

if __name__ == "__main__":
    main()
