import os
import csv
import sys
import re

# Pre-compile regular expressions for speed
re_gene_id = re.compile(r'gene_id "([^"]+)"')
re_gene_name = re.compile(r'gene_name "([^"]+)"')

def gtf_to_bed(gtf_file, bed_file):
    """Convert a GTF file to BED format and save the output (optimized)."""
    buffer = []
    buffer_size = 10000  # Batch size for writing, you can adjust

    with open(gtf_file, 'r') as infile, open(bed_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')

        for line in infile:
            if line[0] == '#':
                continue

            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue

            if fields[2] != 'gene':
                continue

            # Use 0-based BED format
            chrom = 'chr' + fields[0]
            start = str(int(fields[3]) - 1)
            end = fields[4]
            strand = fields[6]
            attributes = fields[8]

            # Fast regex extraction
            gene_id_match = re_gene_id.search(attributes)
            gene_name_match = re_gene_name.search(attributes)
            gene_id = gene_id_match.group(1).split('.')[0] if gene_id_match else ''
            gene_name = gene_name_match.group(1) if gene_name_match else ''

            buffer.append([chrom, start, end, gene_name, gene_id, strand])

            # Write in batches
            if len(buffer) >= buffer_size:
                writer.writerows(buffer)
                buffer.clear()

        if buffer:
            writer.writerows(buffer)

def process_all_gtf_files(input_folder, output_folder):
    """Find all GTF files in the given folder and convert them to BED format."""
    os.makedirs(output_folder, exist_ok=True)
    gtf_files = [f for f in os.listdir(input_folder) if f.endswith('.gtf')]

    if not gtf_files:
        print("No GTF files found in the input directory.")
        return

    for gtf_file in gtf_files:
        input_gtf_path = os.path.join(input_folder, gtf_file)
        output_bed_path = os.path.join(output_folder, gtf_file.replace('.gtf', '.bed'))
        print(f"Processing: {gtf_file} → {output_bed_path}")
        gtf_to_bed(input_gtf_path, output_bed_path)

    print(f"Conversion complete. BED files saved to {output_folder}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python gtf2bed.py <input_directory> <output_directory>")
        sys.exit(1)

    input_directory = sys.argv[1]
    output_directory = sys.argv[2]
    process_all_gtf_files(input_directory, output_directory)
