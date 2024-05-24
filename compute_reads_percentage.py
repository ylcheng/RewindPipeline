import gzip
import bz2
import subprocess
import argparse
import time

def count_reads_in_fastq_gz(file_path):
    """
    Count the number of reads in a .fastq.gz file.

    Parameters
    ----------
    file_path : str
        Path to the .fastq.gz file.

    Returns
    -------
    int
        Number of reads in the .fastq.gz file.
    """
    with gzip.open(file_path, 'rt') as f:
        line_count = sum(1 for _ in f)
    return line_count // 4

def count_lines(file_path):
    """
    Count the number of lines in a file (plain text, gzipped, or bzipped).

    Parameters
    ----------
    file_path : str
        Path to the file.

    Returns
    -------
    int
        Number of lines in the file.
    """
    if file_path.endswith('.gz'):
        with gzip.open(file_path, 'rt') as f:
            line_count = sum(1 for _ in f)
    elif file_path.endswith('.bz2'):
        with bz2.open(file_path, 'rt') as f:
            line_count = sum(1 for _ in f)
    else:
        with open(file_path, 'r') as f:
            line_count = sum(1 for _ in f)
    return line_count

def compute_percentage(original_count, processed_count):
    """
    Compute the percentage of processed reads relative to original reads.

    Parameters
    ----------
    original_count : int
        Number of reads in the original file.
    processed_count : int
        Number of reads in the processed file.

    Returns
    -------
    float
        Percentage of processed reads relative to original reads.
    """
    if original_count == 0:
        return 0.0
    return (processed_count / original_count) * 100

def main(original_file_path, processed_file_path):
    """
    Main function to compute the percentage of reads in the processed file.

    Parameters
    ----------
    original_file_path : str
        Path to the original file.
    processed_file_path : str
        Path to the processed file.
    """
    if original_file_path.endswith('.fastq.gz'):
        original_count = count_reads_in_fastq_gz(original_file_path)
    else:
        original_count = count_lines(original_file_path)
    processed_count = count_lines(processed_file_path)
    percentage = compute_percentage(original_count, processed_count)

    print(f"Number of reads in original file: {original_count}")
    print(f"Number of reads in processed file: {processed_count}")
    print(f"Percentage of processed reads: {percentage:.2f}%")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute the percentage of processed reads relative to the original reads in a file.")
    parser.add_argument('--original', required=True, help="Path to the original file (fastq.gz or text file with reads per line).")
    parser.add_argument('--processed', required=True, help="Path to the processed file (text, gzipped, or bzipped).")

    args = parser.parse_args()

    starttime = time.time()

    main(args.original, args.processed)

    endtime = time.time()
    seconds = endtime - starttime
    print(f'---------- {seconds} seconds ----------')

