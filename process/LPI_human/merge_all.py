import numpy as np
from tqdm import tqdm
import pandas as pd

# Load the CSV files
NPInter_inter = pd.read_csv("../../data/NPInter_LPI/human/correct_NPInter_LPI.csv")
LncBook_inter = pd.read_csv("../../data/LncBook_LPI/LncBook_LPI.csv")
mapping_file = pd.read_csv("mapped_lncRNA.txt", sep='\t', header=None, names=['lncbook_id', 'npinter_id'])  # Contains gene_id and symbol mappings

# Load the BED files
NPInter_Lnc = pd.read_csv("npinter_lnc_0based.bed", sep='\t', header=None, names=['chr', 'start', 'end', 'gene_id', 'symbol', 'strand', 'identifier'])
LncBook_Lnc = pd.read_csv("lncbook_lnc_0based.bed", sep='\t', header=None, names=['chr', 'start', 'end', 'gene_id', 'symbol', 'strand', 'identifier'])

inter_cols = ['gene_id', 'symbol', 'tarName','tissueOrCell']

# Select relevant columns from the NPInter data
NPInter_inter = NPInter_inter[inter_cols]

# Rename the column of LncBook_inter
LncBook_inter.columns = inter_cols

# Replace commas with semicolons in the 'tissueOrCell' column of NPInter_inter
NPInter_inter['tissueOrCell'] = NPInter_inter['tissueOrCell'].str.replace(',', ';', regex=False)

# Merge coordination of lncRNAs.
NPInter_inter = pd.merge(NPInter_inter, NPInter_Lnc, on=['gene_id', 'symbol'], how='inner')
LncBook_inter = pd.merge(LncBook_inter, LncBook_Lnc, on=['gene_id', 'symbol'], how='inner')

# Define functions to get the most common values
get_common_gene_id = lambda x: x.mode()[0] if not x.isna().all() else None

# Apply transform to keep original structure
NPInter_inter['identifier'] = NPInter_inter.groupby(['chr', 'start', 'end'])['identifier'].transform(get_common_gene_id)

# Create a mapping dictionary from the mapping file for quick lookup
# Each lncbook_id points to a set of symbols it maps to
mapping_dict = {}
for _, row in mapping_file.iterrows():
    mapping_dict[row['npinter_id']] = row['lncbook_id']

NPInter_inter['identifier'] = NPInter_inter['identifier'].replace(mapping_dict)

# Initialize progress bar for the LncBook_inter apply function
tqdm.pandas()

# Apply the update logic to LncBook_inter and add a column 'Keep' to indicate unmatched rows
LncBook_inter['Keep'] = False

Lnc_map = {}
def update_row(row):

    inter_matches = pd.DataFrame()

    if row['identifier'] in Lnc_map:
        inter_matches = Lnc_map[row['identifier']]

    elif row['identifier'] in mapping_dict:
        inter_matches = NPInter_inter[NPInter_inter['identifier'] == row['identifier']]  

    if not inter_matches.empty: 
        Lnc_map[row['identifier']] = inter_matches

        # Filter rows with the same 'tarName' value
        protein_matches = inter_matches[inter_matches['tarName'] == row['tarName']]

        # Check if 'tissueOrCell' value of the row is contained in the 'tissueOrCell' column of filtered rows
        if not protein_matches.empty:
            tissue_matches = protein_matches[protein_matches['tissueOrCell'].apply(lambda x: row['tissueOrCell'] in x if x != '-' else False)]
            if not tissue_matches.empty:
                return row
    row['Keep'] = True
    return row

# Apply the update_symbol function to update the 'symbol' column in LncBook_inter with a progress bar
LncBook_inter = LncBook_inter.progress_apply(lambda row: update_row(row), axis=1)
LncBook_inter = LncBook_inter[LncBook_inter['Keep']]

NPInter_inter['lncRNA_ID'] = NPInter_inter['identifier'].apply(lambda x: 'l' + str(x))
LncBook_inter['lncRNA_ID'] = LncBook_inter['identifier'].apply(lambda x: 'l' + str(x))

# Combine rows from NPInter_inter with unmatched rows from LncBook_inter
result_df = pd.concat([NPInter_inter, LncBook_inter], ignore_index=True)

lnc_columns = ['lncRNA_ID', 'identifier', 'gene_id', 'symbol', 'chr', 'start', 'end', 'strand']
lnc = result_df[lnc_columns]
lnc = lnc.drop_duplicates(subset=['lncRNA_ID'])

# Save the lncRNA
lnc.to_csv("../../data/LPI/human/lncRNA.csv", index=False)

result_df['protein_ID'] = result_df['tarName'].apply(lambda x: 'p' + str(x))

pro_columns = ['protein_ID', 'tarName']
protein = result_df[pro_columns].drop_duplicates()
protein.columns = ['protein_ID', 'protein']
# Save the protein
protein.to_csv("../../data/LPI/human/protein.csv", index=False)

result_df = result_df[['lncRNA_ID', 'protein_ID', 'tissueOrCell']]

# Save the final result
result_df.to_csv("../../data/LPI/human/LPI.csv", index=False)
