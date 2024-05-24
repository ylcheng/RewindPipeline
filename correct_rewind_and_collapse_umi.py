#!/usr/bin/env python3

import pandas as pd
import argparse
import numpy as np

def read_starcode_file(file_path):
    """
    Read the starcode.txt file into a DataFrame.

    Parameters
    ----------
    file_path : str
        Path to the starcode.txt file.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the data from the file.
    """
    # Read the first line of the file to check for the header
    with open(file_path, 'r') as file:
        first_line = file.readline().strip()

    # Define the expected header
    expected_header = ['cell_barcode', 'canonical_sequence', 'sequence_in_cluster']

    # Check if the first line matches the expected header
    if first_line.split('\t') == expected_header:
        return pd.read_csv(file_path, sep='\t')
    else:
        return pd.read_csv(file_path, sep='\t', header=None, names=expected_header)

def read_umi_counts_file(file_path):
    """
    Read the UMI_counts.txt file into a DataFrame, keeping only the first three columns.

    Parameters
    ----------
    file_path : str
        Path to the UMI_counts.txt file.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the data from the file.
    """
    return pd.read_csv(file_path, sep='\t', header=None, names=['cell_barcode', 'lineage_barcode', 'umi_count'], usecols=[0, 1, 2])

def correct_lineage_barcodes(umi_df, starcode_df):
    """
    Correct lineage barcodes in the UMI_counts DataFrame using the starcode DataFrame.

    Parameters
    ----------
    umi_df : pd.DataFrame
        DataFrame containing the UMI counts data.
    starcode_df : pd.DataFrame
        DataFrame containing the starcode data.

    Returns
    -------
    pd.DataFrame
        DataFrame with corrected lineage barcodes.
    """
    # Create a dictionary for correction using (cell_barcode, sequence_in_cluster) as key
    correction_dict = dict(zip(zip(starcode_df['cell_barcode'], starcode_df['sequence_in_cluster']), starcode_df['canonical_sequence']))

    # Apply correction to the lineage barcodes
    umi_df['corrected_lineage_barcode'] = umi_df.apply(
        lambda row: correction_dict.get((row['cell_barcode'], row['lineage_barcode']), row['lineage_barcode']),
        axis=1
    )

    return umi_df

def process_corrected_umi_counts(corrected_umi_df):
    """
    Process the corrected UMI counts DataFrame to combine and sort by UMI counts.

    Parameters
    ----------
    corrected_umi_df : pd.DataFrame
        DataFrame with corrected lineage barcodes.

    Returns
    -------
    pd.DataFrame
        DataFrame with combined and sorted UMI counts.
    pd.DataFrame
        DataFrame with total molecules and lineage barcodes per cell.
    """
    # Combine rows with the same cell barcode and corrected lineage barcode
    combined_df = corrected_umi_df.groupby(['cell_barcode', 'corrected_lineage_barcode'], as_index=False).agg({
        'umi_count': 'sum'
    })

    # Sort by UMI counts from large to small
    sorted_df = combined_df.sort_values(by='umi_count', ascending=False)

    # Compute total molecules per cell
    total_molecules_df = sorted_df.groupby('cell_barcode', as_index=False)['umi_count'].sum().rename(columns={'umi_count': 'total_molecules_per_cell'})

    # Compute total unique lineage barcodes per cell
    total_lineages_df = sorted_df.groupby('cell_barcode', as_index=False)['corrected_lineage_barcode'].nunique().rename(columns={'corrected_lineage_barcode': 'total_lineages_per_cell'})

    # Create total_df by merging total molecules and total lineage barcodes per cell
    total_df = pd.merge(total_molecules_df, total_lineages_df, on='cell_barcode')

    # Sort total_df by total molecules per cell and then by lineage barcodes per cell from large to small
    total_df = total_df.sort_values(by=['total_molecules_per_cell', 'total_lineages_per_cell'], ascending=False)

    # Merge sorted_df with total_df by cell barcode
    final_df = pd.merge(sorted_df, total_df, on='cell_barcode')

    # Check for cell barcodes loss
    initial_cell_barcodes = set(corrected_umi_df['cell_barcode'])
    final_cell_barcodes = set(final_df['cell_barcode'])
    lost_barcodes = initial_cell_barcodes - final_cell_barcodes

    if lost_barcodes:
        print(f"Warning: The following cell barcodes were lost during the merge process: {lost_barcodes}")

    return final_df

def main(starcode_file_path, umi_counts_file_path, output_file_path):
    """
    Main function to correct lineage barcodes and process UMI counts.

    Parameters
    ----------
    starcode_file_path : str
        Path to the starcode.txt file.
    umi_counts_file_path : str
        Path to the UMI_counts.txt file.
    output_file_path : str
        Path to the output file where the corrected data will be saved.
    """
    # Read input files
    starcode_df = read_starcode_file(starcode_file_path)
    umi_counts_df = read_umi_counts_file(umi_counts_file_path)

    # Correct lineage barcodes
    corrected_umi_counts_df = correct_lineage_barcodes(umi_counts_df, starcode_df)

    # Process the corrected UMI counts
    final_df = process_corrected_umi_counts(corrected_umi_counts_df)

    # percent of cells with only one lineage barcode
    fraction_of_one_event = np.sum(final_df['total_lineages_per_cell'] == 1) / final_df.shape[0]
    print(f'Percentage of cells with one lineage barcode: {fraction_of_one_event * 100:.2f}%')

    # Save the corrected data to the output file
    final_df.to_csv(output_file_path, sep='\t', index=False)
    print(f"Corrected data saved to {output_file_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Correct lineage barcodes using starcode data and process UMI counts.")
    parser.add_argument("--starcode_file", required=True, help="Path to the starcode.txt file.")
    parser.add_argument("--umi_counts_file", required=True, help="Path to the UMI_counts.txt file.")
    parser.add_argument("--output_file", required=True, help="Path to the output file where the corrected data will be saved.")

    args = parser.parse_args()
    main(args.starcode_file, args.umi_counts_file, args.output_file)

