from Bio import SeqIO
import pandas as pd
import os
import re

# Define sequence length limit
MAX_SEQUENCE_LENGTH = 20000  # Remove genes if any transcript exceeds this length

# Load the mapping file and initialize dictionaries
def load_transcript_ids(mapping_file):
    df = pd.read_csv(mapping_file, header=None, sep='\t', names=['lncRNA_ID', 'transcript_ID'])
    
    gene_to_transcripts = {}  # Store mapping from gene ID to transcript ID list
    transcript_to_gene = {}   # Store mapping from transcript ID to gene ID

    for _, row in df.iterrows():
        gene_id = row['lncRNA_ID'].strip()
        transcript_id = row['transcript_ID'].strip()

        if gene_id not in gene_to_transcripts:
            gene_to_transcripts[gene_id] = set()
        gene_to_transcripts[gene_id].add(transcript_id)
        transcript_to_gene[transcript_id] = gene_id

    return gene_to_transcripts, transcript_to_gene

# Fetch sequences for transcript IDs from a list of FASTA files
def fetch_sequences(fasta_files, gene_to_transcripts):
    found_transcripts = {}  # Store found transcript sequences
    missing_genes = set()   # Store genes whose sequences are missing
    oversized_genes = set() # Store genes that contain transcripts exceeding MAX_SEQUENCE_LENGTH

    # Create a set of all transcript IDs that need to be found
    all_transcript_ids = {tid for tids in gene_to_transcripts.values() for tid in tids}
    needed_ids = all_transcript_ids.copy()  # Copy set to track missing IDs

    for file in fasta_files:
        print(f"Checking file: {file}")  # Debugging output
        if not needed_ids:  # Stop early if all sequences have been found
            break
        for record in SeqIO.parse(file, "fasta"):
            if record.id in needed_ids:
                if len(record.seq) > MAX_SEQUENCE_LENGTH:
                    oversized_genes.add(transcript_to_gene[record.id])  # Mark gene for removal
                else:
                    found_transcripts[record.id] = record.seq  # Store valid transcript
                needed_ids.remove(record.id)  # Remove found ID from the set

    # Identify genes for which all transcripts are missing
    for gene_id, transcript_ids in gene_to_transcripts.items():
        if not any(tid in found_transcripts for tid in transcript_ids):  # If all transcripts of a gene are missing
            missing_genes.add(gene_id)

    # Combine genes to remove: those missing and those with oversized transcripts
    genes_to_remove = missing_genes.union(oversized_genes)

    # Remove the affected genes and their transcripts
    for gene_id in genes_to_remove:
        for tid in gene_to_transcripts[gene_id]:  # Iterate through all transcripts of the gene
            found_transcripts.pop(tid, None)  # Ensure all associated transcripts are removed

    # Filter out removed genes from the gene-to-transcript dictionary
    filtered_gene_to_transcripts = {gene: trans for gene, trans in gene_to_transcripts.items() if gene not in genes_to_remove}

    print(f"Total found sequences: {len(found_transcripts)}")
    print(f"Total missing genes removed: {len(missing_genes)}")
    print(f"Total oversized genes removed: {len(oversized_genes)}")

    return found_transcripts, filtered_gene_to_transcripts

# Write the found sequences to a new FASTA file
def write_fasta(sequences, output_file):
    with open(output_file, "w") as f:
        for trans_id, seq in sequences.items():
            SeqIO.write(SeqIO.SeqRecord(seq, id=trans_id, description=""), f, "fasta")
        print(f"Written {len(sequences)} sequences to {output_file}")  # Debugging output

# Write the filtered gene-to-transcript mappings to a new file
def write_filtered_mapping(filtered_gene_to_transcripts, output_mapping_file):
    with open(output_mapping_file, "w") as f:
        for gene_id, transcript_ids in filtered_gene_to_transcripts.items():
            for transcript_id in transcript_ids:
                f.write(f"{gene_id}\t{transcript_id}\n")
    print(f"Written filtered gene-transcript mappings to {output_mapping_file}")  # Debugging output

# File paths and execution
mapping_file = "lnc_trans.txt"
base_directory = "../../reference_lncRNA/human/fasta/"
processed_ensembl_dir = "../../reference_lncRNA/human/fasta/processed_ensembl/"

def extract_version(filename):
    #match = re.search(r'v(\d+)', filename)
    match = re.search(r'GRCh38\.(\d+).ncrna_processed\.fa', filename)
    return int(match.group(1)) if match else -1  # Extract version number, default to -1 if no match

# Get all .fa files in processedensembl directory
processed_ensembl_fa_files = [f for f in os.listdir(processed_ensembl_dir) if f.endswith(".fa")]
sorted_ensembl_files = [os.path.join(processed_ensembl_dir, f) for f in sorted(processed_ensembl_fa_files, key=extract_version, reverse=True)]

fasta_files = [
    os.path.join(base_directory, "LncBookv2_OnlyLnc.fa"),
    os.path.join(base_directory, "NONCODEv6_human_processed.fa"),
    os.path.join(base_directory, "NONCODEv5_human_processed.fa")
] + sorted_ensembl_files  # Add dynamically found files

output_fasta_file = "transcript_sequences.fasta"
output_mapping_file = "filtered_lnc_trans.txt"

# Load required transcript IDs and gene mappings
gene_to_transcripts, transcript_to_gene = load_transcript_ids(mapping_file)

# Fetch sequences and remove missing/oversized genes
sequences, filtered_gene_to_transcripts = fetch_sequences(fasta_files, gene_to_transcripts)

# Write the filtered sequences to output FASTA file
write_fasta(sequences, output_fasta_file)

# Write the filtered gene-transcript mappings to a new file
write_filtered_mapping(filtered_gene_to_transcripts, output_mapping_file)
