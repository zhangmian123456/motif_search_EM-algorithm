from algorithm_function import *

# 提取数据集
proteinList = preprocess(proteinsDataset)

# 读取已经训练好的权重
weightMatrixList = weight_matrix_read(weightFile)

for i in range(0, trainTime):
    # 存储不同训练motif对应最高分数片断
    motif_save(proteinList, weightMatrixList, i)
