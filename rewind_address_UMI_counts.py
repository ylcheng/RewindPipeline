#!/usr/bin/env python3
import argparse
import pandas as pd

def process_rewind_address(input_file, output_file):
    """
    Processes the input file to generate an output file with the following columns:
    cell barcode, lineage barcode, the number of unique UMIs for each cell barcode-lineage barcode pair,
    total unique UMIs per cell barcode, and total unique lineage barcodes per cell barcode.
    The output is sorted by the cell barcode with the highest total number of unique UMIs,
    and then within each cell barcode group, it is sorted by the lineage barcode with the highest number of unique UMIs.

    Parameters
    ----------
    input_file : str
        Path to the input file containing cell barcode, UMI, lineage barcode, and read counts.
    output_file : str
        Path to the output file where the results will be saved.

    Returns
    -------
    None
    """
    # Read the input file into a DataFrame
    df = pd.read_csv(input_file, sep='\t', header=None, names=['cellbc', 'umi', 'lineage_barcode', 'read_counts'])

    # Group by cell barcode and lineage barcode, and count unique UMIs
    grouped = df.groupby(['cellbc', 'lineage_barcode']).agg(unique_umi_count=('umi', 'nunique')).reset_index()

    # Sort by unique UMIs per lineage barcode within each cell barcode
    grouped = grouped.sort_values(by=['unique_umi_count'], ascending=False)

    # Calculate the total number of unique lineage barcodes per cell barcode
    cellbc_total_lineage_barcodes = grouped.groupby('cellbc')['lineage_barcode'].agg("nunique").reset_index().rename(columns={"lineage_barcode": "total_unique_lineage_barcodes"})

    # Calculate the total number of unique UMIs per cell barcode
    cellbc_total_umis = grouped.groupby('cellbc')['unique_umi_count'].sum().reset_index().rename(columns={'unique_umi_count': 'total_unique_umis'})

    # Merge the totals with the grouped data
    totals_df = pd.merge(cellbc_total_umis, cellbc_total_lineage_barcodes, how="outer", on="cellbc")

    # Sort the totals DataFrame by total unique UMIs and total unique lineage barcodes
    totals_df = totals_df.sort_values(by=totals_df.columns[1:].to_list(), ascending=False)

    # Merge the totals DataFrame back with the grouped data
    grouped = grouped.merge(totals_df, on='cellbc')

    # Check if any NaN values exist in the entire DataFrame
    # If no NaN values are found, it means all cell barcodes are accounted for in the processing
    has_nan = grouped.isnull().values.any()
    print(f"NaN values found in DataFrame: {has_nan}")  # Output: True if any NaN values are found, otherwise False


    # Write the output to a new file
    grouped.to_csv(output_file, sep='\t', index=False, header=False)

    print(f"Processing complete. Output saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Process input file to generate a sorted output file with cell barcode, lineage barcode, unique UMI count, total unique UMIs per cell barcode, and total unique lineage barcodes per cell barcode.')
    parser.add_argument('--input', type=str, required=True, help='Path to the input file containing cell barcode, UMI, lineage barcode, and read counts.')
    parser.add_argument('--output', type=str, required=True, help='Path to the output file.')
    args = parser.parse_args()

    process_rewind_address(args.input, args.output)

if __name__ == "__main__":
    main()

