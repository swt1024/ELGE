#!/bin/bash
# Download, decompress and clean up Ensembl FASTA files using version numbers

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
base_url="https://ftp.ensembl.org/pub/release-%s/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz"

for ver in "${versions[@]}"
do
    url=$(printf "$base_url" "$ver")
    gz_name="Homo_sapiens.GRCh38.v${ver}.ncrna.fa.gz"
    fa_name="Homo_sapiens.GRCh38.v${ver}.ncrna.fa"

    echo "Downloading $url as $gz_name"
    wget -c "$url" -O "$gz_name"

    echo "Decompressing $gz_name to $fa_name"
    gunzip -c "$gz_name" > "$fa_name"

    echo "Deleting $gz_name"
    rm -f "$gz_name"
done

echo "All files downloaded, renamed, and extracted."
