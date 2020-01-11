from algorithm_function import *

# 提取数据集
proteinList = preprocess(proteinsDataset)

# 训练
maxSeqDict, maxSeqList, weightMatrixList = train(proteinList, trainTime)

# 存储理想最高分数
maxScoreSeqList, motifScoreList = max_score_save(resultFile, maxSeqDict, proteinList, weightMatrixList)

# 按照匹配得分排序并保存
match_score_save(matchScoreFile, maxScoreSeqList, motifScoreList)

# 存储训练得到的权重
weight_matrix_save(weightFile, weightMatrixList)

# 读取已经训练好的权重
weightMatrixList = weight_matrix_read(weightFile)

for i in range(0, trainTime):
    # 存储不同训练motif对应最高分数片断
    motif_save(proteinList, weightMatrixList, i)
