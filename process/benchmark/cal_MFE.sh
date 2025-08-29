#!/usr/bin/env bash
set -euo pipefail

# -------------------------------
# Usage check
# -------------------------------
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_fasta_file> <output_MFE_file>"
    exit 1
fi

input_fasta="$1"
output_file="$2"

# -------------------------------
# Preconditions
# -------------------------------
if [ ! -f "$input_fasta" ]; then
    echo "Error: Input FASTA file does not exist: $input_fasta"
    exit 1
fi

if ! command -v RNAfold >/dev/null 2>&1; then
    echo "Error: RNAfold not found in PATH."
    exit 1
fi

if ! command -v parallel >/dev/null 2>&1; then
    echo "Error: GNU parallel not found in PATH."
    exit 1
fi

# -------------------------------
# Split FASTA into per-transcript files
# Each output file: temp_fasta/<id>.fasta
# -------------------------------
rm -rf temp_fasta
mkdir -p temp_fasta

awk '
BEGIN { fname="" }
/^>/ {
    # when hitting a new header, flush previous sequence if any
    if (fname != "" && seq != "") {
        print seq >> fname
        close(fname)
    }
    # derive id (strip leading ">" and anything after first whitespace; also strip suffix after first dot)
    header=$0
    id=substr($1,2)
    sub(/\..*/, "", id)
    fname=sprintf("temp_fasta/%s.fasta", id)
    print header > fname
    seq=""
    next
}
{
    # accumulate sequence lines verbatim
    seq = (seq == "" ? $0 : seq "\n" $0)
}
END {
    if (fname != "" && seq != "") {
        print seq >> fname
        close(fname)
    }
}' "$input_fasta"

# -------------------------------
# Write CSV header
# -------------------------------
echo "transcript_id,MFE" > "$output_file"

# -------------------------------
# Parallel RNAfold MFE computation
# For each per-transcript FASTA:
# - Extract transcript id from filename
# - Run RNAfold --noPS
# - Parse the MFE from the structure line: ... (-12.30)
# -------------------------------
export LC_ALL=C

find "temp_fasta" -type f -name "*.fasta" | parallel --jobs 16 --bar '
    f="{}"
    id="$(basename "$f" .fasta)"
    # RNAfold prints sequence and a structure+energy line; grab the energy field (last field) and strip parentheses
    mfe="$(RNAfold --noPS < "$f" | awk "/\\)/{print \$NF; exit}" | tr -d "()")"
    echo "${id},${mfe}"
' >> "$output_file"

# -------------------------------
# Cleanup
# -------------------------------
rm -rf temp_fasta

echo "Final MFE calculations have been saved to: $output_file"
