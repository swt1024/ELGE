import pandas as pd

# 读取 CSV 文件
df = pd.read_csv('./MLP_predictions_heart.csv')

# 计算第三列为 1 的数量（假设第三列是索引 2）
count = (df.iloc[:, 2] == 1).sum()

# 输出结果
print(f"The number of rows where the third column is 1: {count}")
