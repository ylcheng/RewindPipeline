import argparse

def filter_by_mismatch(input_path, output_path, max_mismatches):
    """
    Filters records based on the number of mismatches, keeping only those with mismatches less than
    or equal to the specified maximum.

    Parameters
    ----------
    input_path : str
        Path to the file containing records with mismatch counts.
    output_path : str
        Path where the filtered records will be saved.
    max_mismatches : int
        Maximum allowed mismatches for the records to be kept.

    Returns
    -------
    None
    """
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            parts = line.strip().split()
            if len(parts) < 5:
                continue  # Skip lines that don't have enough data
            mismatches = int(parts[4])  # Mismatch count is the fourth element
            if mismatches <= max_mismatches:
                outfile.write(line)  # Write the line to the output file if it passes the filter

    print(f"Records filtered by mismatches. Results saved to {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Filter records based on mismatch counts.")
    parser.add_argument('--input', type=str, required=True, help="Path to the input file.")
    parser.add_argument('--output', type=str, required=True, help="Path to the output file.")
    parser.add_argument('--max_mismatches', type=int, required=True, help="Maximum allowed mismatches.")
    
    args = parser.parse_args()
    filter_by_mismatch(args.input, args.output, args.max_mismatches)

if __name__ == "__main__":
    main()
