import os
import sys
from Bio import SeqIO

def process_fasta(input_file, output_file):
    """Process a single FASTA file: clean transcript IDs and save to output."""
    with open(output_file, "w") as output_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            # Remove the transcript version (keep only the ID before the first '.')
            cleaned_id = record.id.split('.')[0]
            record.id = cleaned_id
            record.description = ''  # Clear description to keep only ID
            SeqIO.write(record, output_handle, "fasta")
    
    print(f"Processed: {input_file} → Saved to: {output_file}")

def process_all_fasta_files(input_folder, output_folder):
    """Find all FASTA (.fa) files in input_folder and process them into output_folder."""
    
    # Ensure output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # List all FASTA files in the input directory
    fasta_files = [f for f in os.listdir(input_folder) if f.endswith('.fa')]

    if not fasta_files:
        print(f"No FASTA files found in {input_folder}")
        return

    for fasta_file in fasta_files:
        input_fasta_path = os.path.join(input_folder, fasta_file)
        output_fasta_path = os.path.join(output_folder, fasta_file.replace('.fa', '_processed.fa'))

        process_fasta(input_fasta_path, output_fasta_path)

    print(f"\nAll FASTA files processed. Output saved to: {output_folder}")

if __name__ == "__main__":

    if len(sys.argv) != 3:
        print("Usage: python process_fasta.py <input_directory> <output_directory>")
        sys.exit(1)

    input_directory = os.path.abspath(sys.argv[1])
    output_directory = os.path.abspath(sys.argv[2])

    # Process all FASTA files in the input directory
    process_all_fasta_files(input_directory, output_directory)
