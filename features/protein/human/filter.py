import pandas as pd
import re

# List of 30 representative GTEx tissues
tissue_list = [
    'Brain_Cortex','Brain_Cerebellum','Brain_Hippocampus','Brain_Amygdala','Brain_Hypothalamus',
    'Heart_Left_Ventricle','Heart_Atrial_Appendage','Artery_Aorta','Artery_Coronary','Artery_Tibial',
    'Lung','Muscle_Skeletal','Skin_Sun_Exposed_Lower_leg','Skin_Not_Sun_Exposed_Suprapubic',
    'Whole_Blood','Spleen','Cells_Cultured_fibroblasts','Thyroid','Pituitary','Adrenal_Gland',
    'Liver','Stomach','Colon_Transverse','Small_Intestine_Terminal_Ileum','Esophagus_Mucosa',
    'Pancreas','Kidney_Cortex','Ovary','Testis','Prostate'
]

# === Step 1: Load GTEx GCT file (skip first 2 metadata rows) ===
exp_df = pd.read_csv(
    'GTEx_Analysis_2022-06-06_v10_RNASeQCv2.4.2_gene_median_tpm.gct',
    sep='\t', skiprows=2
)

# Keep only Ensembl ID, gene name, and selected tissues
exp_df = exp_df[['Name','Description'] + tissue_list]
exp_df = exp_df.rename(columns={'Name':'ensembl_id','Description':'gene_name'})

# === Step 2: Clean Ensembl IDs (remove version and _PAR_Y suffix) ===
def clean_id(x):
    if not isinstance(x, str):
        return x

    # Remove _PAR_Y suffix, e.g. ENSG00000123456.11_PAR_Y → ENSG00000123456
    x = x.replace("_PAR_Y", "")
    return x

exp_df["ensembl_id"] = exp_df["ensembl_id"].map(clean_id)

# === Step 3: Merge duplicated rows by cleaned Ensembl ID ===
# For expression columns: sum values (X + Y PAR copies)
expr_cols = tissue_list

agg_dict = {col: "sum" for col in expr_cols}
agg_dict["gene_name"] = "first"

# Group by cleaned Ensembl ID and aggregate
exp_df = (
    exp_df.groupby("ensembl_id", as_index=False)
    .agg(agg_dict)
)

# === Step 4: Save the merged expression matrix ===
ordered_cols = ["ensembl_id", "gene_name"] + tissue_list
exp_df = exp_df[ordered_cols]
exp_df.to_csv('filtered_exp.csv', index=False)

print(f"[DONE] Wrote {exp_df.shape[0]} genes to filtered_exp.csv")
