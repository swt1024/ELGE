#!/bin/bash

# Ensure the script stops if any command fails
set -e

# Input BED files
BED_LncBook="./lncbook_lnc_0based.bed"  # Replace with the path to your first BED file
BED_NPInter="./gencode/npinter_lnc_0based.bed"  # Replace with the path to your second BED file

# Output files
MAPPED_OUTPUT="mapped_lncRNA_temp.bed"
EXTRACTED_COLUMNS_OUTPUT="./gencode/mapped_lncRNA.txt"

# Perform 100% mapping using bedtools intersect with -f 1.0 -r for reciprocal full overlap
echo "Performing 100% mapping using bedtools..."
bedtools intersect -a $BED_LncBook -b $BED_NPInter -f 1.0 -r -s -wo > $MAPPED_OUTPUT

echo "Mapped results saved to $MAPPED_OUTPUT."

# Extract specific columns from the mapped output
echo "Extracting specific columns..."
cut -f7,14 $MAPPED_OUTPUT > $EXTRACTED_COLUMNS_OUTPUT

rm mapped_lncRNA_temp.bed

echo "Extracted columns saved to $EXTRACTED_COLUMNS_OUTPUT."

echo "Process completed successfully!"
