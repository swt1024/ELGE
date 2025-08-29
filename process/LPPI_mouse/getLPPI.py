import pandas as pd

# File paths
lpi_file = '../../data/LPI/mouse/LPI.csv'
lpi_protein_file = '../../data/LPI/mouse/protein.csv'
ppi_file = '../../data/PPI/mouse/PPI.csv'
output_file = '../../data/LPPI/mouse/LPPI.csv'

# Load LPI and PPI files
lpi_data = pd.read_csv(lpi_file)
lpi_proteins = pd.read_csv(lpi_protein_file)
ppi_data = pd.read_csv(ppi_file)

# Extract unique protein names from and PPI
ppi_proteins = pd.concat([ppi_data['Protein A'], ppi_data['Protein B']]).drop_duplicates().reset_index(drop=True).to_frame(name='protein')
# Prefix protein IDs with 'p'
ppi_proteins['protein_ID'] = ['p' + str(x) for x in ppi_proteins['protein']]

# Save protein mappings to new files
proteins = pd.concat([ppi_proteins, lpi_proteins]).reset_index(drop=True)
proteins = proteins.drop_duplicates()
proteins.to_csv('../../data/LPPI/mouse/protein.csv', index=False)

ppi_data = pd.merge(ppi_data, ppi_proteins, left_on='Protein A', right_on='protein', how='left').drop(columns=['Protein A', 'protein'])
ppi_data = ppi_data.rename(columns={'protein_ID': 'protein_ID_A'})
ppi_data = pd.merge(ppi_data, ppi_proteins, left_on='Protein B', right_on='protein', how='left').drop(columns=['Protein B', 'protein'])
ppi_data = ppi_data.rename(columns={'protein_ID': 'protein_ID_B'})

# Process LPI data to standardize the format as two columns: Node_i, Node_j
lpi_edges = lpi_data.rename(columns={'lncRNA_ID': 'Node_i', 'protein_ID': 'Node_j'})[['Node_i', 'Node_j']]

# Process PPI data to standardize the format as two columns: Node_i, Node_j
ppi_edges = ppi_data.rename(columns={'protein_ID_A': 'Node_i', 'protein_ID_B': 'Node_j'})[['Node_i', 'Node_j']]

# Combine all edges into a single DataFrame
lppi = pd.concat([lpi_edges, ppi_edges]).reset_index(drop=True)

# Save the final LPPI file
lppi.to_csv(output_file, index=False)

print(f"Processing complete! The final LPPI file has been saved to {output_file}")
