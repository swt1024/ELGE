#!/bin/bash
# Download, decompress and clean up Ensembl GTF files using version numbers

# List all release numbers you want to download
versions=(
84
85
86
87
88
89
90
91
92
93
94
95
96
97
98
99
100
)

# URL template: just replace version number
base_url="https://ftp.ensembl.org/pub/release-%s/gtf/mus_musculus/Mus_musculus.GRCm38.%s.gtf.gz"

for ver in "${versions[@]}"
do
    url=$(printf "$base_url" "$ver" "$ver")
    echo "Downloading $url"
    wget -c "$url"
done

# Decompress all .gtf.gz files in current folder
for gzfile in *.gtf.gz
do
    echo "Decompressing $gzfile"
    gunzip "$gzfile"
done

# Remove all original .gtf.gz files (should already be gone after gunzip, but for safety)
rm -f *.gtf.gz

echo "All files downloaded and extracted."
