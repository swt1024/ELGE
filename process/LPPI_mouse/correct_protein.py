import pandas as pd

# === Configuration file paths ===
mapping_path = "right_protein.csv"         
protein_path = "../../data/LPPI/mouse/protein.csv"        
lppi_path = "../../data/LPPI/mouse/LPPI.csv"               

# === Read mapping table ===
df_map = pd.read_csv(mapping_path, header=None, names=["old_symbol", "new_symbol"])
symbol_map = dict(zip(df_map["old_symbol"], df_map["new_symbol"]))

# === Process protein file ===
df_protein = pd.read_csv(protein_path)
df_protein["NewSymbol"] = df_protein["protein"].map(symbol_map).fillna(df_protein["protein"])
df_protein["NewProteinID"] = "p" + df_protein["NewSymbol"]

# Save the updated protein file
update_protein = df_protein[["NewSymbol", "NewProteinID"]]
update_protein.columns = ['protein', 'protein_ID']
update_protein = update_protein.drop_duplicates()
update_protein.to_csv("../../data/LPPI/mouse/protein_updated.csv", index=False)
print("Saved: protein_updated.csv")

# === Process LPPI file ===
df_lppi = pd.read_csv(lppi_path)

# Build protein ID map: p+old_symbol → p+new_symbol
protein_id_map = {
    "p" + old: "p" + new
    for old, new in zip(df_map["old_symbol"], df_map["new_symbol"])
}

# Replace IDs in LPPI where node IDs start with 'p'
def replace_protein_id(value):
    return protein_id_map.get(value, value) if value.startswith("p") else value

df_lppi["Node_i"] = df_lppi["Node_i"].apply(replace_protein_id)
df_lppi["Node_j"] = df_lppi["Node_j"].apply(replace_protein_id)

# Save the updated LPPI file
df_lppi.to_csv("../../data/LPPI/mouse/LPPI_updated.csv", index=False)
print("Saved: LPPI_updated.csv (includes LPI + PPI)")
