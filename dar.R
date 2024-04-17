# 加载必要的库
library(seqinr)  # 用于处理fasta格式的序列
library(cluster)  # 用于K均值聚类
library(stringr)
library(ggplot2)


fasta_file_path <- "SELEX-13-LFJ14896_L3_1stastic.fasta.sort.fasta"

# 读取FASTA文件内容
fasta_content <- readLines(fasta_file_path)



# 解析fasta文件内容，提取序列部分
sequences <- fasta_content[grep("^[ACGTN]+$", fasta_content)]

# 创建特征矩阵，表示每种碱基的频率
feature_matrix <- matrix(0, nrow = length(sequences), ncol = 4)

for (i in 1:length(sequences)) {
  seq <- sequences[i]
  A_count <- sum(str_count(seq, "A"))
  C_count <- sum(str_count(seq, "C"))
  G_count <- sum(str_count(seq, "G"))
  T_count <- sum(str_count(seq, "T"))
  feature_matrix[i, ] <- c(A_count, C_count, G_count, T_count)
}





# 1. 肘部法则（Elbow Method）：计算误差平方和
wss <- numeric(length = 20)
for (i in 1:20) {
  kmeans_result <- kmeans(feature_matrix, centers = i)
  wss[i] <- kmeans_result$tot.withinss
}

# 绘制肘部法则图表
plot(1:20, wss, type = "b", xlab = "Number of Clusters", ylab = "Within-cluster Sum of Squares")











# 执行K均值聚类
k <- 13  # 聚类数
kmeans_result <- kmeans(feature_matrix, centers = k)




# 输出每个序列的聚类结果
cluster_results <- data.frame(Sequence = sequences, Cluster = kmeans_result$cluster)
print(cluster_results)

file_path <- "data.txt"
write.table(cluster_results, file = file_path, row.names = TRUE)

# 创建一个数据框，用于绘制聚类图
plot_data <- data.frame(PC1 = feature_matrix[, 1], PC2 = feature_matrix[, 2], Cluster = factor(kmeans_result$cluster))

# 使用ggplot2创建散点图
shape_mapping <- data.frame(Cluster = levels(plot_data$Cluster), Shape = c(19, 17, 15, 5, 1,12,10,13,20,14,16,18,21))

ggplot(plot_data, aes(x = PC1, y = PC2, color = Cluster,shape = Cluster)) +
  geom_point(size= 4) +
  labs(title = "Clustering of Sequences") +
  scale_shape_manual(values = shape_mapping$Shape) +
  theme_minimal()


# 创建散点图，每个点代表一个数据点，颜色表示聚类分配
ggplot(cluster_results, aes(x = Cluster)) +
  geom_bar() +
  labs(x = "Cluster", y = "Count") +
  ggtitle(paste("Clustering Results (Optimal K =", 7, ")"))+
  theme(axis.text = element_text(size = 16),  # 设置坐标轴标签字体大小
        axis.title = element_text(size = 18))  # 设置坐标轴标题字体大小


install.packages("ape")
library(ape)

# 计算序列之间的距离矩阵
distance_matrix <- dist(feature_matrix)
# 使用Neighbor-Joining方法构建树
nj_tree <- nj(distance_matrix)

# 可视化进化树
plot(nj_tree)


################



library(phangorn)
install.packages("phangorn")
# 创建数据框，包括序列和标识
fasta_file_path <- "SELEX-13-LFJ14896_L3_1stastic.fasta.sort.fasta"

# 读取FASTA文件内容
data <- readLines(fasta_file_path)

sequences <- data[seq(2, length(data), by = 2)]
identifiers <- data[seq(1, length(data), by = 2)]

df <- data.frame(Sequence = sequences, Identifier = identifiers, stringsAsFactors = FALSE)

# 进行特征提取，这里以计算碱基频率为例
# 你可以根据需要选择其他特征
A_freq <- sapply(df$Sequence, function(seq) sum(str_count(seq, "A")))
C_freq <- sapply(df$Sequence, function(seq) sum(str_count(seq, "C")))
G_freq <- sapply(df$Sequence, function(seq) sum(str_count(seq, "G")))
T_freq <- sapply(df$Sequence, function(seq) sum(str_count(seq, "T")))

df$A_frequency <- A_freq
df$C_frequency <- C_freq
df$G_frequency <- G_freq
df$T_frequency <- T_freq

df <- df[, -1]

# 创建一个距离矩阵
dist_matrix <- dist(df[, -1])

# 构建分类树
phylo_tree <- NJ(dist_matrix)  # 使用邻接法构建树，可以根据需要选择其他方法

# 可视化分类树
plot(phylo_tree, show.tip.label = TRUE)





dist_matrix <- dist(df[, -1])

cexRow <- 1.5  # 行数字的大小
cexCol <- 1.5  # 列数字的大小

# 创建热图
heatmap(as.matrix(dist_matrix), Rowv = NA, Colv = NA, col = colorRampPalette(c("white", "red"))(100), 
        scale = "none", main = "Sequence Similarity Heatmap", xlab = "Sequence", ylab = "Sequence",
        cexRow = cexRow, cexCol = cexCol)



