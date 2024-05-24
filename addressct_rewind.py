#!/usr/bin/env python3
import sys
import argparse
import gzip
import bz2

def open_file(filename):
    """
    Opens a file with support for .gz, .bz2, and .tsv files.
    
    Parameters
    ----------
    filename : str
        Path to the input file.
    
    Returns
    -------
    file object
        Opened file object.
    """
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    elif filename.endswith('.bz2'):
        return bz2.open(filename, 'rt')
    else:
        return open(filename, 'r')

def collapse_counts(input_file, ct_outfile):
    """
    Collapses read addresses based on UMI sequences without accounting for errors.
    
    Parameters
    ----------
    input_file : str
        Path to the input file containing read addresses.
    ct_outfile : str
        Path to the output file where collapsed counts will be saved.
    
    Returns
    -------
    None
    """
    address_dict = {}
    
    with open_file(input_file) as infile:
        for line in infile:
            llist = line.strip().split()
            if len(llist) == 5:
                address = "\t".join(llist[1:4])  # Construct the address from cell_barcode, umi, and lineage_barcode
                mismatches = int(llist[4])
                
                if address in address_dict:
                    address_dict[address] += 1
                else:
                    address_dict[address] = 1
    
    with open(ct_outfile, 'w') as output:
        for address, count in address_dict.items():
            output.write(f"{address}\t{count}\n")

def main():
    parser = argparse.ArgumentParser(description='Collapse read addresses based on UMI sequences without accounting for errors.')
    parser.add_argument('--input', type=str, required=True, help='Path to the input file (.gz, .bz2, or .tsv)')
    parser.add_argument('--output', type=str, required=True, help='Path to the output file')
    args = parser.parse_args()
    
    collapse_counts(args.input, args.output)

if __name__ == "__main__":
    main()
