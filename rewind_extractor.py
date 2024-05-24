#!/usr/bin/env python3
import argparse
import gzip
from Bio import SeqIO

CAPTURE_SEQ2 = "GCTCACCTATTAGCGGCTAAGG"

# Define the sets for W, S, and N conditions
W = set('AT')
S = set('GC')
N = set('ATGC')  # Though we won't actually check N as all are allowed

# Function to compute mismatches against the expected WSN pattern
def count_wsn_mismatches(sequence):
    mismatches = 0
    for i in range(0, len(sequence), 3):
        if i < len(sequence) and sequence[i] not in W:
            mismatches += 1
        if i + 1 < len(sequence) and sequence[i + 1] not in S:
            mismatches += 1
        # N is skipped as it can be any nucleotide
    return mismatches

def rewind_extractor(fastq1_path, fastq2_path, output_path, instrument_run_flowcell_ID):
    """
    Extracts cell barcodes, UMIs, and 84 nt lineage barcodes from read1 and read2 FASTQ files, respectively,
    and writes them to an output file along with the number of mismatches to the WSN pattern.

    Parameters
    ----------
    fastq1_path : str
        Path to the read1 FASTQ file which contains cell barcodes and UMIs.
    fastq2_path : str
        Path to the read2 FASTQ file which contains the 84 nt lineage barcodes.
    output_path : str
        Path to the output file where results will be saved.
    num_lines : int
        Number of lines to process for testing.

    Output
    ------
    The output file will contain the following tab-separated fields:
    - readid: The identifier for the read.
    - cell_barcode: The cell barcode extracted from read1.
    - umi: The UMI extracted from read1.
    - lineage_barcode: The 84 nt lineage barcode extracted from read2.
    - mismatches: The number of mismatches between the lineage barcode and the expected WSN pattern.

    Example
    -------
    Given an input read1 FASTQ file and read2 FASTQ file, the function extracts relevant data and writes it
    to the output file with the following format:
    readid    cell_barcode    umi    lineage_barcode    mismatches
    """
    cell_barcode_len = 16
    umi_len = 12

    # Lists to store cell barcodes, UMIs, and read IDs
    cell_barcodes = []
    umis = []
    read_ids = []

    print("Extract from Read1 ------------------------------")
    # Extract cell barcodes, UMIs, and read IDs from read1
    with gzip.open(fastq1_path, 'rt') as f1:
        for record in SeqIO.parse(f1, "fastq"):
            read_id = record.id.split(instrument_run_flowcell_ID)[1]
            cell_barcode = str(record.seq[:cell_barcode_len])
            umi = str(record.seq[cell_barcode_len:cell_barcode_len + umi_len])
            cell_barcodes.append(cell_barcode)
            umis.append(umi)
            read_ids.append(read_id)

    print("Extract from Read2 ------------------------------")
    # Extract lineage barcodes from read2 and write to the output file
    number_of_readids_unmatched = 0
    number_of_lineage_barcodes_mismatched = 0
    number_of_reads_with_captureseq2 = 0
    number_reads_in_address = 0
    with gzip.open(fastq2_path, 'rt') as f2, open(output_path, 'w') as out:
        for i, record in enumerate(SeqIO.parse(f2, "fastq")):
            read_id = record.id.split(instrument_run_flowcell_ID)[1]
            read2_sequence = str(record.seq)
            capture_seq2_index = read2_sequence.find(CAPTURE_SEQ2)
            if capture_seq2_index != -1:
                number_of_reads_with_captureseq2 += 1
                lineage_barcode_start = capture_seq2_index - 84
                if lineage_barcode_start >= 0:
                    lineage_barcode = read2_sequence[lineage_barcode_start:capture_seq2_index]
                    mismatches = count_wsn_mismatches(lineage_barcode)
                    if mismatches > 0:
                        number_of_lineage_barcodes_mismatched += 1
                    if read_id == read_ids[i]:
                        out.write(f"{read_id}\t{cell_barcodes[i]}\t{umis[i]}\t{lineage_barcode}\t{mismatches}\n")
                        number_reads_in_address += 1
                    else:
                        number_of_readids_unmatched += 1

    print(f'Number of lineage barcodes with at least one mismatch to the WSN pattern: {number_of_lineage_barcodes_mismatched}')
    print(f'Number of reads in Address: {number_reads_in_address}')
    print(f'Fraction of mismatched lineage barcode: {number_of_lineage_barcodes_mismatched/number_reads_in_address}')
    print(f'Number of reads containing CaptureSeq2: {number_of_reads_with_captureseq2}')
    print(f'Number of readids did not match between Read1 and Read2: {number_of_readids_unmatched}')
    print(f"Extraction complete. Results saved to {output_path}")

def main():
    parser = argparse.ArgumentParser(description='Extract cell barcodes, UMIs, and lineage barcodes from FASTQ files.')
    parser.add_argument('--fastq1', type=str, required=True, help='Path to the read1 FASTQ file')
    parser.add_argument('--fastq2', type=str, required=True, help='Path to the read2 FASTQ file')
    parser.add_argument('--output', type=str, required=True, help='Path to the output file')
    parser.add_argument('--instrument-run-flowcell-ID', type=str, required=True, help='instrument, run, and flowcell IDs of the fastq read header separated by :, ej AV100007:PY055:2336402118:')
    args = parser.parse_args()

    rewind_extractor(args.fastq1, args.fastq2, args.output, args.instrument_run_flowcell_ID)

if __name__ == "__main__":
    main()

