from Bio import pairwise2
from Bio import SeqIO
from sklearn.cluster import AgglomerativeClustering
import matplotlib.pyplot as plt
import numpy as np

# 读取FASTA文件
fasta_file = "SELEX-13-LFJ14896_L3_1stastic.fasta.sort.fasta"
records = list(SeqIO.parse(fasta_file, "fasta"))

# 计算相似性矩阵
similarity_matrix = np.zeros((len(records), len(records)))

for i in range(len(records)):
    for j in range(i, len(records)):
        # 使用比对分数作为相似性度量
        alignments = pairwise2.align.globalxx(records[i].seq, records[j].seq, one_alignment_only=True)
        similarity_score = alignments[0].score
        similarity_matrix[i, j] = similarity_score
        similarity_matrix[j, i] = similarity_score

# 进行层次聚类
clustering = AgglomerativeClustering(n_clusters=3, affinity='precomputed', linkage='average')
clustering.fit(similarity_matrix)

# 绘制聚类图
plt.figure(figsize=(8, 6))
plt.scatter(range(len(records)), [0] * len(records), c=clustering.labels_, cmap='rainbow')
plt.title("Sequence Clustering")
plt.show()