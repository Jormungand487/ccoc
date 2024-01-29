import pandas as pd
from sklearn.metrics import confusion_matrix

# 从两个CSV文件导入数据
file1 = 'label.csv'
file2 = 'pred.csv'

df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)

# 假设两个文件中都有一个名为'label'的列，表示标签
actual_labels = df1
predicted_labels = df2

# 计算混淆矩阵
conf_matrix = confusion_matrix(actual_labels, predicted_labels)

# 打印混淆矩阵
print("Confusion Matrix:")
print(conf_matrix)
