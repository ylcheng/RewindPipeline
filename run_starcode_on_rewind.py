#!/usr/bin/env python3
import os
import subprocess
import tempfile
import argparse

def run_starcode(input_data, starcode_path, max_distance):
    """
    Runs Starcode on the provided input data and returns the output lines.

    Parameters
    ----------
    input_data : str
        The input data to be processed by Starcode.
    starcode_path : str
        The path to the Starcode executable.
    max_distance : int
        The maximum Levenshtein distance for clustering.

    Returns
    -------
    list
        The lines of output from Starcode.
    """
    # Create temporary input and output files for Starcode
    with tempfile.NamedTemporaryFile(delete=False) as temp_input_file:
        temp_input_file.write(input_data.encode())
        temp_input_name = temp_input_file.name

    temp_output_name = temp_input_name + "_output"

    # Create temporary log file
    with tempfile.NamedTemporaryFile(delete=False) as temp_log_file:
        log_file_name = temp_log_file.name


    try:
        # 1. Run Starcode with the specified parameters and redirect output to the log file
        with open(log_file_name, "w") as log_output:
            subprocess.run(
                [starcode_path, "-i", temp_input_name, "-o", temp_output_name, "-d", str(max_distance), "--seq-id"],
                stdout=log_output,
                stderr=log_output,
                check=True  # This will raise a CalledProcessError if the command fails
            )

        # 2. Read Starcode output
        with open(temp_output_name, "r") as temp_output_file:
            output_lines = temp_output_file.readlines()

        # 3. Clean up temporary files
        os.remove(temp_input_name)
        os.remove(temp_output_name)

        # 4. If no exception was raised, delete the log file
        os.remove(log_file_name)

    except subprocess.CalledProcessError as e:
        # Handle the error: print an error message and do not delete the log file
        print(f"Starcode failed with error: {e}")
        print(f"Log file: {log_file_name}")
        raise  # Re-raise the exception to indicate failure

    finally:
        # Ensure temporary files are cleaned up if an exception occurs
        if os.path.exists(temp_input_name):
            os.remove(temp_input_name)
        if os.path.exists(temp_output_name):
            os.remove(temp_output_name)

    return output_lines

def parse_starcode_output(line):
    """
    Parses a line of Starcode output to extract the canonical sequence,
    cluster size, and sequence indices.

    Parameters
    ----------
    line : str
        A single line of Starcode output.

    Returns
    -------
    dict
        A dictionary containing the canonical sequence, cluster size, and sequence indices.
    """
    # Split the line into components
    parts = line.strip().split('\t')

    # Extract the canonical sequence, cluster size, and sequence indices
    canonical_sequence = parts[0]
    cluster_size = int(parts[1])
    sequence_indices = list(map(int, parts[2].split(',')))

    return {
        "canonical_sequence": canonical_sequence,
        "cluster_size": cluster_size,
        "sequence_indices": sequence_indices
    }

def process_umi_counts(input_file, output_file, starcode_path, max_distance):
    """
    Processes the UMI counts, runs Starcode for lineage barcodes within each cell barcode,
    and stores the output in a combined results file.

    Parameters
    ----------
    input_file : str
        Path to the input file containing cell barcode, UMI, lineage barcode, and read counts.
    output_file : str
        Path to the output file where the results will be saved.
    starcode_path : str
        Path to the Starcode executable.
    max_distance : int
        The maximum Levenshtein distance for clustering.

    Returns
    -------
    None
    """
    # Ensure the output file is empty before starting and write the column names
    with open(output_file, "w") as f:
        f.write("cell_barcode\tcanonical_sequence\sequence_in_cluster\n")

    # Read the input file and process each cell barcode
    with open(input_file, "r") as file:
        lines = file.readlines()

    # Dictionary to hold lineage barcodes and UMI counts for each cell barcode
    cell_data = {}

    # Process each line to populate the cell_data dictionary
    for line in lines:
        cell_barcode, lineage_barcode, umi_count, molecules_per_cell_count, lineage_barcode_per_cell_count = line.strip().split()
        umi_count = int(umi_count)
        if cell_barcode not in cell_data:
            cell_data[cell_barcode] = []
        cell_data[cell_barcode].extend([lineage_barcode] * umi_count)

    # Process each cell barcode
    for cell_barcode, lineage_barcodes in cell_data.items():
        # Prepare input data for Starcode by joining lineage barcodes with newlines
        input_data = "\n".join(lineage_barcodes) + "\n"

        # Run Starcode on the prepared input data
        starcode_output = run_starcode(input_data, starcode_path, max_distance)

        # Parse Starcode output and append to the final output file
        with open(output_file, "a") as f:
            for line in starcode_output:
                parsed_output = parse_starcode_output(line)
                canonical_sequence = parsed_output["canonical_sequence"]
                sequence_indices = parsed_output["sequence_indices"]

                for index in sequence_indices:
                    lineage_barcode = lineage_barcodes[index - 1]
                    f.write(f"{cell_barcode}\t{canonical_sequence}\t{lineage_barcode}\n")

    print(f"Processing complete. Results saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Process UMI counts and run Starcode for lineage barcode clustering within each cell barcode.')
    parser.add_argument('--input', type=str, required=True, help='Path to the input file containing cell barcode, UMI, lineage barcode, and read counts.')
    parser.add_argument('--output', type=str, required=True, help='Path to the output file.')
    parser.add_argument('--starcode_path', type=str, required=True, help='Path to the Starcode executable.')
    parser.add_argument('--max_distance', type=int, required=True, help='Maximum Levenshtein distance for clustering.')

    args = parser.parse_args()

    process_umi_counts(args.input, args.output, args.starcode_path, args.max_distance)

if __name__ == "__main__":
    main()

