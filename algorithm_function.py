import csv
import json
from random import randint as randint
from math import log as log

proteinsDataset = "proteins.csv"

resultFile = "result.txt"

matchScoreFile = "match_score.txt"

weightFile = "weight.json"

trainNum = 2000

trainTime = 1000

motif_length = 7


def csv_reader(path):
    """读取csv文件"""
    origin_protein_list = []
    with open(path) as f:
        reader = csv.reader(f)
        for item in reader:
            origin_protein_list.append(item)
    return origin_protein_list


def preprocess(path):
    """数据预处理，去掉序号"""
    protein_order_list = csv_reader(path)[1:]
    protein_list = []
    for protein in protein_order_list:
        protein_list.append(protein[1:])
    return protein_list


def create_amino_list(protein_list):
    """找出所有氨基酸种类"""
    amino_species_list = []
    for protein in protein_list:
        for amino in protein:
            if amino not in amino_species_list:
                amino_species_list.append(amino)
    return amino_species_list


def transpose(matrix):
    """转置氨基酸序列矩阵"""
    new_matrix = []
    for i in range(len(matrix[0])):
        matrix1 = []
        for j in range(len(matrix)):
            matrix1.append(matrix[j][i])
        new_matrix.append(matrix1)
    return new_matrix


def pssm_algorithm(train_seqs):
    """计算权重矩阵"""
    # 氨基酸每个位置出现次数字典
    amino_count_dict = {}

    # 氨基酸出现频率字典
    amino_sum_fre_dict = {}

    # 氨基酸权重矩阵字典
    amino_weight_dict = {}

    amino_species_list = ['A', 'C', 'D', 'E', 'F',
                          'G', 'H', 'I', 'K', 'L',
                          'M', 'N', 'P', 'Q', 'R',
                          'S', 'T', 'V', 'W', 'Y']

    # 所有字典初始化
    for amino in amino_species_list:
        amino_count_dict[amino] = [0] * motif_length
        amino_sum_fre_dict[amino] = 0
        amino_weight_dict[amino] = [0] * motif_length

    # 转置氨基酸序列矩阵
    amino_list_trans = transpose(train_seqs)

    # 统计氨基酸每个位置出现次数
    t = 0
    for pos_seq in amino_list_trans:
        for amino in pos_seq:
            i = t
            amino_count_dict[amino][i] += 1
        t += 1

    # 统计氨基酸出现频率
    for amino in amino_count_dict:
        amino_sum = 0
        for aminoCount in amino_count_dict[amino]:
            amino_sum += aminoCount
        amino_sum_fre = amino_sum / (len(train_seqs) * motif_length)
        amino_sum_fre_dict[amino] = amino_sum_fre

    # 获得氨基酸权重矩阵
    for i in range(0, motif_length):
        for amino in amino_count_dict:
            amino_pos_fre = amino_count_dict[amino][i] / len(train_seqs)
            if amino_sum_fre_dict[amino] == 0:
                amino_weight = 0
            else:
                amino_weight = amino_pos_fre / amino_sum_fre_dict[amino]
            if amino_weight == 0:
                amino_weight = -10
            else:
                amino_weight = log(amino_weight, 10)
            amino_weight_dict[amino][i] += float(format(amino_weight, '.4f'))
    return amino_weight_dict


def random_seq_generation(protein_list):
    """从2000组蛋白质序列中，每组随机取一个motif_length长度的氨基酸链"""
    train_seq = []
    for i in range(0, trainNum):
        head = randint(0, 43)
        tail = head + motif_length
        seq = protein_list[i][head:tail]
        train_seq.append(seq)
    return train_seq


def score_cal(test_seq, weight_matrix):
    """根据测试序列和权重矩阵计算分数"""
    score = 0
    for i in range(0, motif_length):
        score += weight_matrix[test_seq[i]][i]
    return score


def train_seq_generation(protein_list, weight_matrix):
    """从每个蛋白质序列中计算所有motif_length氨基酸片断分数，选出分数最高的一个"""
    new_train_seqs = []
    max_seq_index_list = []
    for protein in protein_list:
        score_list = []
        for i in range(0, 44):
            test_seq = protein[i:i+motif_length]
            score_list.append(score_cal(test_seq, weight_matrix))
        max_seq_index = score_list.index(max(score_list))
        new_train_seq = protein[max_seq_index:max_seq_index+motif_length]
        new_train_seqs.append(new_train_seq)
        max_seq_index_list.append(max_seq_index)
    return new_train_seqs, max_seq_index_list


def max_score_seq_generation(weight_matrix):
    """返回最高分数片断，即理想化motif"""
    max_score_seq = []
    for i in range(0, motif_length):
        score_list = []
        for amino in weight_matrix:
            score_list.append(weight_matrix[amino][i])
        for amino in weight_matrix:
            if max(score_list) == weight_matrix[amino][i]:
                max_score_seq.append(amino)
    return max_score_seq


def bubble_sort_3(list1, list2, list3):
    """训练结果排序"""
    list_len = len(list1)-1
    for i in range(list_len):
        for j in range(list_len - i):
            if list1[j] < list1[j + 1]:
                list1[j], list1[j + 1] = list1[j + 1], list1[j]
                list2[j], list2[j + 1] = list2[j + 1], list2[j]
                list3[j], list3[j + 1] = list3[j + 1], list3[j]

    return list1, list2, list3


def result_check(protein_list, weight_matrix_list, seq_index):
    """返回motif和对应序号"""
    motif_list, max_seq_index = train_seq_generation(protein_list, weight_matrix_list[seq_index])
    return motif_list, max_seq_index


def train(protein_list, tt):
    """训练函数"""
    max_seq_list = []
    max_score_list = []
    weight_matrix_list = []
    max_seq_dict = {}

    for random_t in range(0, tt):
        # 训练矩阵
        train_seq = random_seq_generation(protein_list)
        weight_matrix = pssm_algorithm(train_seq)
        for train_t in range(0, 10):
            new_train_seqs, index = train_seq_generation(protein_list, weight_matrix)
            weight_matrix = pssm_algorithm(new_train_seqs)

        max_score_seq = max_score_seq_generation(weight_matrix)
        max_score = score_cal(max_score_seq, weight_matrix)
        max_score = float(format(max_score, '.4f'))

        max_score_list.append(max_score)
        max_seq_list.append(max_score_seq)
        weight_matrix_list.append(weight_matrix)

    max_score_list, max_seq_list, weight_matrix_list = bubble_sort_3(max_score_list, max_seq_list, weight_matrix_list)
    # 按照训练效果进行排序

    for i in range(0, tt):
        print(i)
        # 将得分与训练得到的motif对应成字典
        if max_score_list[i] in max_seq_dict:
            # 防止出现相同得分
            max_seq_dict[max_score_list[i]+0.000001] = max_seq_list[i]
        max_seq_dict[max_score_list[i]] = max_seq_list[i]
    return max_seq_dict, max_seq_list, weight_matrix_list


def motif_evaluation(motif_list, max_score_seq):
    evaluation_score = 0
    for motif in motif_list:
        for i in range(len(motif)):
            if motif[i] == max_score_seq[i]:
                evaluation_score += float(100 / motif_length)
    return evaluation_score / len(motif_list)


def motif_list_evaluation(max_seq_dict, protein_list, weight_matrix_list):
    index = -1
    max_score_seq_list = []
    motif_score_list = []
    for max_seq in max_seq_dict:
        index += 1
        motif_list, max_seq_index = train_seq_generation(protein_list, weight_matrix_list[index])
        max_score_seq = max_seq_dict[max_seq]
        motif_score = motif_evaluation(motif_list, max_score_seq)
        max_score_seq_list.append(max_score_seq)
        motif_score_list.append(motif_score)
    return max_score_seq_list, motif_score_list


def max_score_save(file_name, max_seq_dict, protein_list, weight_matrix_list, ):
    index = -1
    max_score_seq_list = []
    motif_score_list = []
    with open(file_name, "w") as max_score_txt:
        for max_seq in max_seq_dict:
            index += 1
            motif_list, max_seq_index = train_seq_generation(protein_list, weight_matrix_list[index])
            max_score_seq = max_score_seq_generation(weight_matrix_list[index])
            motif_score = motif_evaluation(motif_list, max_score_seq)
            max_score_seq_list.append(max_score_seq)
            motif_score_list.append(motif_score)
            max_score_txt.write("No." + str(format(index, '4d'))
                                + "  motif: " + str(max_seq_dict[max_seq])
                                + "  motif_score: " + str(format(max_seq, '.4f'))
                                + "  match_score: " + str(format(motif_score, '.4f')) + '\n')
    return max_score_seq_list, motif_score_list


def match_score_save(file_name, max_score_seq_list, motif_score_list):
    order_list = []
    for i in range(len(max_score_seq_list)):
        order_list.append(i)
    motif_score_list, max_score_seq_list, order_list = bubble_sort_3(motif_score_list, max_score_seq_list, order_list)
    index = -1
    with open(file_name, 'w') as match_score_txt:
        for item in max_score_seq_list:
            index += 1
            match_score_txt.write("seq" + str(format(order_list[index], '4d')) + ' ' + str(item) + ' '
                                  + 'match_score: ' + str(format(motif_score_list[index], '.4f')) + '\n')


def motif_save(protein_list, weight_matrix_list, motif_index):
    motif_list, max_seq_index = train_seq_generation(protein_list, weight_matrix_list[motif_index])
    max_score_seq = max_score_seq_generation(weight_matrix_list[motif_index])
    evaluation_score = motif_evaluation(motif_list, max_score_seq)
    file_name = 'motifs\\' + str(motif_index) + ' ' + str(max_score_seq) + '.txt'
    index = -1
    with open(file_name, "w") as motif_list_txt:
        motif_list_txt.write("motif No." + str(motif_index) + ' ')
        motif_list_txt.writelines(max_score_seq)
        motif_list_txt.write('\n' + "evaluation score: " + str(format(evaluation_score, '.4f')) + '\n')
        for motif in motif_list:
            index += 1
            motif_list_txt.write('\n' + "seq" + str(format(index, '4d')) + '            '
                                 + '(' + str(format(max_seq_index[index], '2d')) + ')' + '              ')
            motif_list_txt.writelines(motif)


def weight_matrix_save(file_name, weight_matrix_list):
    with open(file_name, 'w') as weight_matrix_file:
        json.dump(weight_matrix_list, weight_matrix_file)


def weight_matrix_read(file_name):
    with open(file_name, 'r') as weight_matrix_file:
        weight_matrix = json.load(weight_matrix_file)
    return weight_matrix
