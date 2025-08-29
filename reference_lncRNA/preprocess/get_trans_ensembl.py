import pandas as pd
import os
import sys

def process_gtf(input_file, output_folder):
    """Process a single GTF file and extract gene_id, gene_name, and transcript_id."""
    try:
        # Read GTF file, ignoring comment lines (#)
        df = pd.read_csv(input_file, sep='\t', comment='#', header=None,
                         usecols=[2, 8], names=['type', 'attributes'])
    except Exception as e:
        print(f"Error reading {input_file}: {e}")
        return

    # Filter for transcript entries
    transcripts = df[df['type'] == 'transcript'].copy()

    # Extract gene_id, transcript_id, and gene_name (symbol)
    transcripts.loc[:, 'gene_id'] = transcripts['attributes'].str.extract('gene_id "([^"]+)"')
    transcripts.loc[:, 'symbol'] = transcripts['attributes'].str.extract('gene_name "([^"]+)"')
    transcripts.loc[:, 'transcript_id'] = transcripts['attributes'].str.extract('transcript_id "([^"]+)"')

    transcripts['gene_id'] = transcripts['gene_id'].str.split('.').str[0]
    transcripts['symbol'] = transcripts['symbol'].str.split('.').str[0]
    transcripts['transcript_id'] = transcripts['transcript_id'].str.split('.').str[0]

    # Select relevant columns
    result = transcripts[['gene_id', 'symbol', 'transcript_id']]

    # Define output file name
    base_name = os.path.basename(input_file).rsplit('.', 1)[0]  # Extract file name without extension
    output_file = os.path.join(output_folder, f"{base_name}_trans.txt")

    # Save result to output folder
    result.to_csv(output_file, sep='\t', index=False, header=False)
    print(f"Processed: {input_file} → Saved to: {output_file}")

def process_all_gtf_files(input_folder, output_folder):
    """Find all GTF files in input_folder and convert them to text files in output_folder."""
    
    # Ensure output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # List all GTF files in the input directory
    gtf_files = [f for f in os.listdir(input_folder) if f.endswith('.gtf')]

    if not gtf_files:
        print(f"No GTF files found in {input_folder}")
        return

    for gtf_file in gtf_files:
        input_gtf_path = os.path.join(input_folder, gtf_file)
        process_gtf(input_gtf_path, output_folder)

    print(f"\nAll GTF files processed. Output files saved to: {output_folder}")

if __name__ == "__main__":

    # Check if the correct number of arguments is provided
    if len(sys.argv) != 3:
        print("Usage: python get_trans_ensembl.py <input_directory> <output_directory>")
        sys.exit(1)

    # Retrieve the input and output directory paths from command line arguments
    gtf_folder = sys.argv[1]
    trans_dir = sys.argv[2]

    process_all_gtf_files(gtf_folder, trans_dir)
