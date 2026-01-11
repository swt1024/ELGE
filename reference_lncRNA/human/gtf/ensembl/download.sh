#!/bin/bash
# Download multiple Ensembl GTF files using version numbers

# List all release numbers you want to download
versions=(
76
78
80
84
87
93
104
106
107
108
109
110
111
112
113
)

# URL template: just replace version number
# %s will be replaced by version number in the loop
base_url="https://ftp.ensembl.org/pub/release-%s/gtf/homo_sapiens/Homo_sapiens.GRCh38.%s.gtf.gz"

for ver in "${versions[@]}"
do
    url=$(printf "$base_url" "$ver" "$ver")
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