import pandas as pd
import os
import re

# Set the directory containing Ensembl data
ensembl_dir = "../../reference_lncRNA/mouse/bed/ensembl/"

inter_npinter = pd.read_csv('../../data/NPInter_LPI/mouse/correct_NPInter_LPI.csv')

# Read NONCODE and LncBook BED files
noncodev5_bed = pd.read_csv('../../reference_lncRNA/mouse/bed/NONCODEv5_mm10.lncRNAGene.bed', sep='\t', header=None, names=['chr', 'start', 'end', 'gene_id', 'score', 'strand'])
noncodev6_bed = pd.read_csv('../../reference_lncRNA/mouse/bed/NONCODEv6_mm10.lncRNAGene.bed', sep='\t', header=None, names=['chr', 'start', 'end', 'gene_id', 'score', 'strand'])

noncodev5_bed['gene_id'] = noncodev5_bed['gene_id'].str.split('.').str[0]
noncodev6_bed['gene_id'] = noncodev6_bed['gene_id'].str.split('.').str[0]

# Initialize remaining lncRNA list
remained_inter_npinter = inter_npinter.copy()
remained_inter_npinter['identifier'] = remained_inter_npinter['gene_id']
results = []

# Get genomic position by noncode_id (noncodev6)
pos_lnc_noncodev6 = pd.merge(
    remained_inter_npinter,
    noncodev6_bed[['chr', 'start', 'end', 'strand', 'gene_id']],
    on='gene_id',
    how='inner'
)
results.append(pos_lnc_noncodev6)
remained_inter_npinter = remained_inter_npinter[~remained_inter_npinter['gene_id'].isin(pos_lnc_noncodev6['gene_id'])]

# Get genomic position by noncode_id (noncodev5)
pos_lnc_noncodev5 = pd.merge(
    remained_inter_npinter,
    noncodev5_bed[['chr', 'start', 'end', 'strand', 'gene_id']],
    on='gene_id',
    how='inner'
)
results.append(pos_lnc_noncodev5)
remained_inter_npinter = remained_inter_npinter[~remained_inter_npinter['gene_id'].isin(pos_lnc_noncodev5['gene_id'])]

# **STEP 1: Extract BED file version numbers**
def extract_version(filename):
    match = re.search(r'GRCm38\.(\d+)', filename)  # ensembl
    return int(match.group(1)) if match else -1

# **STEP 2: Retrieve all BED files and sort by version number**
bed_files = [f for f in os.listdir(ensembl_dir) if f.endswith(".bed")]
bed_files_sorted = sorted(bed_files, key=extract_version, reverse=True)  # Sort by version number descending

# **STEP 3: Iterate through BED files in sorted order**
for bed_file in bed_files_sorted:
    bed_path = os.path.join(ensembl_dir, bed_file)
    print(f"Processing {bed_file}...")

    # Read Ensembl BED file
    ensembl_bed = pd.read_csv(bed_path, sep='\t', header=None, names=['chr', 'start', 'end', 'symbol', 'gene_id', 'strand'])

    # Match by gene_id
    pos_lnc_ensembl = pd.merge(remained_inter_npinter, ensembl_bed[['chr', 'start', 'end', 'strand', 'gene_id']], on='gene_id', how='inner')
    results.append(pos_lnc_ensembl)
    remained_inter_npinter = remained_inter_npinter[~remained_inter_npinter['gene_id'].isin(pos_lnc_ensembl['gene_id'])]

    # Match by symbol
    pos_lnc_symbol = pd.merge(remained_inter_npinter, ensembl_bed[['chr', 'start', 'end', 'strand', 'symbol']], on='symbol', how='inner')
    remained_inter_npinter = remained_inter_npinter[~remained_inter_npinter['symbol'].isin(pos_lnc_symbol['symbol'])]
    pos_lnc_symbol['identifier'] = pos_lnc_symbol['symbol']
    results.append(pos_lnc_symbol)

# Combine all results
pos_inter_npinter = pd.concat(results, ignore_index=True)
pos_inter_npinter['lncRNA_ID'] = pos_inter_npinter['identifier'].apply(lambda x: 'l' + str(x))
pos_inter_npinter['protein_ID'] = pos_inter_npinter['tarName'].apply(lambda x: 'p' + str(x))

lncRNA_cols = ['lncRNA_ID', 'identifier', 'gene_id', 'symbol', 'chr', 'start', 'end', 'strand']
lnc_npinter = pos_inter_npinter[lncRNA_cols].drop_duplicates(subset=['identifier'])
lnc_npinter.to_csv('../../data/LPI/mouse/lncRNA.csv', index=False)

protein_cols = ['protein_ID', 'tarName']
protein_npinter = pos_inter_npinter[protein_cols].drop_duplicates(subset=['protein_ID'])
protein_npinter.columns = ['protein_ID', 'protein']
protein_npinter.to_csv('../../data/LPI/mouse/protein.csv', index=False)

# Save remaining lncRNAs without genomic positions
remained_inter_npinter[['identifier', 'gene_id', 'symbol']].drop_duplicates().to_csv('lnc_no_pos.csv', index=False)

lpi = pos_inter_npinter[['lncRNA_ID', 'protein_ID']]
lpi.to_csv('../../data/LPI/mouse/LPI.csv', index=False)

print("Processing complete. Results saved.")
