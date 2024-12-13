# Data reading
allData = as.matrix(read.table("k-means_clusters.txt", header = T, row.names = 1))

timePoints = c(-24, 0, 6, 12, 18, 24, 30, 42, 54, 96)
cluster_info <- read.table("clusters_in_waves.txt", 
                           header = T, sep = "\t", row.names = 1)

library(ggplot2)

###
# Cluster 1 
###

gene_list <- row.names(cluster_info)[which(cluster_info$Vawe == 1)]
sum(gene_list %in% row.names(allData)) == length(gene_list)
test = allData[gene_list,]
write.table(test, file = "WaveI.txt", quote = F, sep = "\t")

cluster2 = data.frame(rep(row.names(test), ncol(test)),
                 rep(timePoints, each = nrow(test)), 
                 matrix(test, dim(test)[1]*dim(test)[2], 1))
colnames(cluster2) = c("rownames", "Sample", "expression")
cluster2$TimePoint <- as.numeric(sub("Sample", "", cluster2$Sample))

ggplot(cluster2, aes(TimePoint, expression)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
  geom_line(aes(group = rownames), size = 0.5, alpha = 0.5, color = "grey") + 
#  geom_point(alpha = 0.3, aes(color = expression > 0)) + 
  geom_smooth(size = 2, se = FALSE, color = "black", method = loess) +
  scale_x_continuous(breaks = cluster2$TimePoint) +
  scale_y_continuous(limits = c(-5.5, 7.5)) +
#  scale_color_manual(values = c("brown", "red"), guide = FALSE) +
  labs(title = paste("Expression change in wave", "I"),
       x = "Time point",
       y = "Expression") +
  theme_classic()

###
# Cluster 2
###

gene_list <- row.names(cluster_info)[which(cluster_info$Vawe == 2)]
sum(gene_list %in% row.names(allData)) == length(gene_list)
test = allData[gene_list,]
write.table(test, file = "WaveII.txt", quote = F, sep = "\t")

cluster2 = data.frame(rep(row.names(test), ncol(test)),
                      rep(timePoints, each = nrow(test)), 
                      matrix(test, dim(test)[1]*dim(test)[2], 1))
colnames(cluster2) = c("rownames", "Sample", "expression")
cluster2$TimePoint <- as.numeric(sub("Sample", "", cluster2$Sample))

ggplot(cluster2, aes(TimePoint, expression)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
  geom_line(aes(group = rownames), size = 0.5, alpha = 0.5, color = "grey") + 
  #  geom_point(alpha = 0.3, aes(color = expression > 0)) + 
  geom_smooth(size = 2, se = FALSE, color = "red", method = loess) +
  scale_x_continuous(breaks = cluster2$TimePoint) +
  scale_y_continuous(limits = c(-5.5, 7.5)) +
  #  scale_color_manual(values = c("brown", "red"), guide = FALSE) +
  labs(title = paste("Expression change in wave", "II"),
       x = "Time point",
       y = "Expression") +
  theme_classic()

###
# Cluster 3
###

gene_list <- row.names(cluster_info)[which(cluster_info$Vawe == 3)]
sum(gene_list %in% row.names(allData)) == length(gene_list)
test = allData[gene_list,]
write.table(test, file = "WaveIII.txt", quote = F, sep = "\t")

cluster2 = data.frame(rep(row.names(test), ncol(test)),
                      rep(timePoints, each = nrow(test)), 
                      matrix(test, dim(test)[1]*dim(test)[2], 1))
colnames(cluster2) = c("rownames", "Sample", "expression")
cluster2$TimePoint <- as.numeric(sub("Sample", "", cluster2$Sample))

ggplot(cluster2, aes(TimePoint, expression)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
  geom_line(aes(group = rownames), size = 0.5, alpha = 0.5, color = "grey") + 
  #  geom_point(alpha = 0.3, aes(color = expression > 0)) + 
  geom_smooth(size = 2, se = FALSE, color = "orange", method = loess) +
  scale_x_continuous(breaks = cluster2$TimePoint) +
  scale_y_continuous(limits = c(-5.5, 7.5)) +
  #  scale_color_manual(values = c("brown", "red"), guide = FALSE) +
  labs(title = paste("Expression change in wave", "III"),
       x = "Time point",
       y = "Expression") +
  theme_classic()

###
# Cluster 4
###

gene_list <- row.names(cluster_info)[which(cluster_info$Vawe == 4)]
sum(gene_list %in% row.names(allData)) == length(gene_list)
test = allData[gene_list,]
write.table(test, file = "WaveIV.txt", quote = F, sep = "\t")

# Création de la dataframe pour ggplot2
cluster2 = data.frame(rep(row.names(test), ncol(test)),
                      rep(timePoints, each = nrow(test)), 
                      matrix(test, dim(test)[1]*dim(test)[2], 1))
colnames(cluster2) = c("rownames", "Sample", "expression")
cluster2$TimePoint <- as.numeric(sub("Sample", "", cluster2$Sample))

ggplot(cluster2, aes(TimePoint, expression)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
  geom_line(aes(group = rownames), size = 0.5, alpha = 0.5, color = "grey") + 
  #  geom_point(alpha = 0.3, aes(color = expression > 0)) + 
  geom_smooth(size = 2, se = FALSE, color = "green", method = loess) +
  scale_x_continuous(breaks = cluster2$TimePoint) +
  scale_y_continuous(limits = c(-5.5, 7.5)) +
  #  scale_color_manual(values = c("brown", "red"), guide = FALSE) +
  labs(title = paste("Expression change in wave", "IV"),
       x = "Time point",
       y = "Expression") +
  theme_classic()

###
# Cluster 5
###

gene_list <- row.names(cluster_info)[which(cluster_info$Vawe == 5)]
sum(gene_list %in% row.names(allData)) == length(gene_list)
test = allData[gene_list,]
write.table(test, file = "WaveV.txt", quote = F, sep = "\t")

cluster2 = data.frame(rep(row.names(test), ncol(test)),
                      rep(timePoints, each = nrow(test)), 
                      matrix(test, dim(test)[1]*dim(test)[2], 1))
colnames(cluster2) = c("rownames", "Sample", "expression")
cluster2$TimePoint <- as.numeric(sub("Sample", "", cluster2$Sample))

ggplot(cluster2, aes(TimePoint, expression)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
  geom_line(aes(group = rownames), size = 0.5, alpha = 0.5, color = "grey") + 
  #  geom_point(alpha = 0.3, aes(color = expression > 0)) + 
  geom_smooth(size = 2, se = FALSE, color = "blue", method = loess) +
  scale_x_continuous(breaks = cluster2$TimePoint) +
  scale_y_continuous(limits = c(-5.5, 7.5)) +
  #  scale_color_manual(values = c("brown", "red"), guide = FALSE) +
  labs(title = paste("Expression change in wave", "V"),
       x = "Time point",
       y = "Expression") +
  theme_classic()

