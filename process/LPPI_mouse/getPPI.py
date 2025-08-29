import pandas as pd

# Define the file path
PPI_file = '../../data/raw/BIOGRID-Mus_musculus-4.4.241.tab3.txt'

# Read the file, assuming it is tab-separated
ppi = pd.read_csv(PPI_file, sep='\t')

# Select the required columns, including organism names for filtering
ppi = ppi[['Official Symbol Interactor A', 'Official Symbol Interactor B', 
                      'Organism Name Interactor A', 'Organism Name Interactor B']]

# Filter rows where both interactors belong to the species "Mus musculus"
ppi = ppi[
    (ppi['Organism Name Interactor A'] == 'Mus musculus') &
    (ppi['Organism Name Interactor B'] == 'Mus musculus')
]

# Retain only the columns with the official symbols of the interactors
ppi = ppi[['Official Symbol Interactor A', 'Official Symbol Interactor B']]
ppi.columns = ['Protein A', 'Protein B']

# Convert PPI dataframe columns to strings and strip whitespace
ppi['Protein A'] = ppi['Protein A'].astype(str).str.strip()
ppi['Protein B'] = ppi['Protein B'].astype(str).str.strip()

# Save the filtered interactions
#output_file = '../../data/PPI/human/PPI.csv'
output_file = '../../data/PPI/mouse/PPI.csv'
ppi.to_csv(output_file, index=False)
print(f"Filtered PPI interactions saved to {output_file}")

# Combine both columns into a single column, remove duplicates, and save
protein = pd.concat([ppi["Protein A"], ppi["Protein B"]]).dropna().drop_duplicates().reset_index(drop=True)
protein_df = pd.DataFrame({"protein": protein})

# Save the combined unique interactor list to another file
output_file_protein = "../../data/PPI/mouse/protein_in_ppi.csv"
protein_df.to_csv(output_file_protein, index=False)
print(f"Proteins in PPL saved to {output_file_protein}")
