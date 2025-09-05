import pandas as pd

# import tsv file
df = pd.read_csv("gnomad.v4.1.constraint_metrics.tsv", sep='\t')

# 1. Keep line with canonical == True
df_filtered = df[df["canonical"] == True]

# 2. Group by gene_name
df_grouped = df_filtered.groupby("gene", as_index=False).mean(numeric_only=True)

# 3. Save
df_grouped.to_csv("pLoF_v4.txt", sep='\t', index=False)
