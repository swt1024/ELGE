import pandas as pd

# 读取原始CSV文件
df = pd.read_csv("gnomad.v4.1.constraint_metrics.tsv", sep='\t')

# 1. 仅保留 canonical == True 的行
df_filtered = df[df["canonical"] == True]

# 2. 按 gene 分组，并对其他列取均值
#    注意：这里会自动跳过非数值列，只对数值列取平均
df_grouped = df_filtered.groupby("gene", as_index=False).mean(numeric_only=True)

# 3. 保存到新文件
df_grouped.to_csv("pLoF_v4.txt", sep='\t', index=False)
