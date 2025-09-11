#!/bin/bash

# Check if exactly two arguments are passed
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <bedfile> <outputfile>"
    exit 1
fi

# Assign the arguments to variables
bedfile="$1"
outputfile="$2"

# Check if the provided BED files exist
if [ ! -f "$bedfile" ]; then
    echo "Error: File $bedfile not found."
    exit 1
fi

if [ ! -f "$outputfile" ]; then
    echo "Error: File $outputfile not found."
    exit 1
fi

# Run bedtools intersect to find overlapping regions and process the data
bedtools intersect -a "$bedfile" -b "$bedfile" -s -wo | \
awk '{
    len_a = $3 - $2;  # Length of the feature in the first BED file
    len_b = $9 - $8;  # Length of the feature in the second BED file
    if (len_a == 0 || len_b == 0) next;  # Skip lines with zero-length features
    shorter = (len_a < len_b) ? len_a : len_b;  # Select the shorter length
    if ($NF / shorter >= 0.5 && $4 != $10)  # If overlap ratio is >= 50% and genes are not the same
        print $4, $10;  # Print the overlapping gene names
}' > "$outputfile"  
