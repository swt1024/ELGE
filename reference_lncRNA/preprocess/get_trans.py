import pandas as pd
import os
import sys

def process_gtf(input_file, output_dir):
    # Read the GTF file and extract required columns
    try:
        df = pd.read_csv(input_file, sep='\t', comment='#', header=None,
                         usecols=[2, 8], names=['type', 'attributes'])
    except Exception as e:
        print(f"Error reading {input_file}: {e}")
        return

    # Filter rows where 'type' is 'transcript'
    transcripts = df[df['type'] == 'transcript'].copy()

    # Extract 'gene_id' and 'transcript_id' from the attributes column
    transcripts['gene_id'] = transcripts['attributes'].str.extract('gene_id "([^"]+)"')
    transcripts['transcript_id'] = transcripts['attributes'].str.extract('transcript_id "([^"]+)"')

    # Remove version numbers from IDs (keep only IDs before the first '.')
    transcripts['gene_id'] = transcripts['gene_id'].str.split('.').str[0]
    transcripts['transcript_id'] = transcripts['transcript_id'].str.split('.').str[0]

    # Select only the 'gene_id' and 'transcript_id' columns
    result = transcripts[['gene_id', 'transcript_id']]

    # Define the output filename based on the input filename
    base_name = os.path.basename(input_file)
    output_file = os.path.join(output_dir, f"{base_name.rsplit('.', 1)[0]}_trans.txt")

    # Save the result to a tab-separated file without headers
    result.to_csv(output_file, sep='\t', index=False, header=False)
    print(f"Processed {input_file}, output saved to {output_file}")

def process_all_gtf_files(gtf_folder, output_dir):
    # Create output directory if it does not exist
    os.makedirs(output_dir, exist_ok=True)

    # Find all .gtf files in the specified folder
    gtf_files = [f for f in os.listdir(gtf_folder) if f.endswith('.gtf')]

    if not gtf_files:
        print(f"No .gtf files found in {gtf_folder}")
        return

    # Process each GTF file
    for gtf_file in gtf_files:
        input_gtf_path = os.path.join(gtf_folder, gtf_file)
        process_gtf(input_gtf_path, output_dir)

if __name__ == "__main__":

    # Check if the correct number of arguments is provided
    if len(sys.argv) != 3:
        print("Usage: python get_trans.py <input_directory> <output_directory>")
        sys.exit(1)

    # Retrieve the input and output directory paths from command line arguments
    gtf_folder = sys.argv[1]
    trans_dir = sys.argv[2]

    # Process all GTF files in the specified folder
    process_all_gtf_files(gtf_folder, trans_dir)
