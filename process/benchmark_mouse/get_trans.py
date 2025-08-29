import pandas as pd
import os
import re

# Set the directory containing ensembl transcript data
ensembl_dir = "../../reference_lncRNA/mouse/transcript/ensembl/"

# Read LncBook and NONCODE transcript files
cols = ['gene_id', 'transcript_id']
noncodev5_trans = pd.read_csv('../../reference_lncRNA/mouse/transcript/NONCODEv5_mouse_mm10_lncRNA_trans.txt', sep='\t', header=None, names=cols)

# Read lncRNA ID data
lncRNA = pd.read_csv('../../data/LPI/mouse/lncRNA.csv')
lncRNA = lncRNA[['lncRNA_ID', 'identifier']]

# Initialize remaining lncRNA list
remained_lncRNA = lncRNA[['lncRNA_ID', 'identifier']].copy()

# Function to remove matched rows
def update_remained_lncRNA(df, matched_df):
    """Remove matched lncRNA_ID from the remaining lncRNA list"""
    updated_df = df[~df['lncRNA_ID'].isin(matched_df['lncRNA_ID'])]
    return updated_df[['lncRNA_ID', 'identifier']]

# Store results
results = []

trans_lnc_noncodev5 = pd.merge(remained_lncRNA, noncodev5_trans, left_on='identifier', right_on='gene_id', how='inner')
results.append(trans_lnc_noncodev5)
remained_lncRNA = update_remained_lncRNA(remained_lncRNA, trans_lnc_noncodev5)

# Extract version numbers and sort filenames in descending order
def extract_version(filename):
    match = re.search(r'GRCm38\.(\d+)_trans\.txt', filename)
    return int(match.group(1)) if match else -1  # Extract version number, default to -1 if no match

ensembl_files = [f for f in os.listdir(ensembl_dir) if f.endswith(".txt")]
sorted_ensembl_files = sorted(ensembl_files, key=extract_version, reverse=True)  # Sort by version number (descending)

# Iterate through sorted ensembl transcript files
for txt_file in sorted_ensembl_files:
    file_path = os.path.join(ensembl_dir, txt_file)
    print(f"Processing {txt_file}...")

    # Read ensembl transcript file
    ensembl_trans = pd.read_csv(file_path, sep='\t', header=None, names=['gene_id', 'symbol', 'transcript_id'])

    # Match by gene_id
    trans_lnc_ensembl = pd.merge(remained_lncRNA, ensembl_trans[['gene_id', 'transcript_id']], left_on='identifier', right_on='gene_id', how='inner')
    results.append(trans_lnc_ensembl)
    remained_lncRNA = update_remained_lncRNA(remained_lncRNA, trans_lnc_ensembl)

    # Match by symbol
    trans_lnc_symbol = pd.merge(remained_lncRNA, ensembl_trans[['symbol', 'transcript_id']], left_on='identifier', right_on='symbol', how='inner')
    results.append(trans_lnc_symbol)
    remained_lncRNA = update_remained_lncRNA(remained_lncRNA, trans_lnc_ensembl)

remained_lncRNA.to_csv("no_trans.csv", index=False)

# Combine all results into a single DataFrame
trans_lnc = pd.concat(results, ignore_index=True)
trans_lnc = trans_lnc[['lncRNA_ID', 'transcript_id']]

# Save results to file
trans_lnc.to_csv('lnc_trans.txt', index=False, sep='\t', header=None)

print("Processing complete. Results saved.")
