# Datafile is imported in R
data <- read.table("/expData_allGenes.txt", 
                   header = T, row.names = 1)

# Colnames are changed for easier use
colnames(data) <- unlist(strsplit(colnames(data), "_"))[seq(1, 40, by = 4)]

# The obtained table is verified
head(data)

# The library mixOmics is used
library(mixOmics)

# Experiments have to be in rows and genes in columns
X_data <- t(data)

# Statistics are calculated
result_pca <- pca(X_data, ncomp = nrow(X_data),
                  center = T, scale = T)

# Percentage of variability captured on each component
plot(result_pca, 
     main = "microarray data")

# Projection of time points in the 2D space defined with PC1 and PC2
plotIndiv(result_pca, 
          title = "Transcriptomics - Microarray dataset")

# Projection of time points in the 2D space defined with PC2 and PC3
plotIndiv(result_pca, comp = c(2,3),
          title = "Transcriptomics - Microarray dataset")

# Loadings on PC1, PC2 and PC3
plotLoadings(result_pca, comp = 1)
plotLoadings(result_pca, comp = 2)
plotLoadings(result_pca, comp = 3)

# A table with loadings for all genes is created
table_loadings <- cbind(selectVar(result_pca, comp = 1)$value[row.names(data),], 
                        selectVar(result_pca, comp = 2)$value[row.names(data),],
                        selectVar(result_pca, comp = 3)$value[row.names(data),])

# Rownames and colnames are changed
row.names(table_loadings) <- row.names(data)
colnames(table_loadings) <- c("PC1", "PC2", "PC3")

# Result is verified
head(table_loadings)

# Result writing
write.table(table_loadings, file = "PCA_loadings.txt", 
            quote = F, sep = "\t", row.names = T, col.names = T)

n_genes <- 50

# Variable selection
pca_out1 <- as.matrix(selectVar(result_pca, comp = 1)$value)[1:n_genes,]
write.table(pca_out1, file = "pca_out1.txt",
            quote = F, sep = "\t", col.names = F)

# Positive and negative contributions are separated
pca_out1_pos <- pca_out1[which(pca_out1 > 0)]
pca_out1_neg <- pca_out1[which(pca_out1 < 0)]

gene_list1_pos <- names(pca_out1_pos)
gene_list1_neg <- names(pca_out1_neg)

# Gene Lists are formatted for easy use in Pixel
print(paste(gene_list1_pos, collapse = ";"))
print(paste(gene_list1_neg, collapse = ";"))

# Same selection of variables, with PC2 this time
pca_out2 <- as.matrix(selectVar(result_pca, comp = 2)$value)[1:n_genes,]
write.table(pca_out2, file = "pca_out2.txt",
            quote = F, sep = "\t", col.names = F)

pca_out2_pos <- pca_out2[which(pca_out2 > 0)]
pca_out2_neg <- pca_out2[which(pca_out2 < 0)]

gene_list2_pos <- names(pca_out2_pos)
gene_list2_neg <- names(pca_out2_neg)

print(paste(gene_list2_pos, collapse = ";"))
print(paste(gene_list2_neg, collapse = ";"))

# Finally, selection with PC3
pca_out3 <- as.matrix(selectVar(result_pca, comp = 3)$value)[1:n_genes,]
write.table(pca_out3, file = "pca_out3.txt",
            quote = F, sep = "\t", col.names = F)

pca_out3_pos <- pca_out3[which(pca_out3 > 0)]
pca_out3_neg <- pca_out3[which(pca_out3 < 0)]

gene_list3_pos <- names(pca_out3_pos)
gene_list3_neg <- names(pca_out3_neg)

print(paste(gene_list3_pos, collapse = ";"))
print(paste(gene_list3_neg, collapse = ";"))

# Selection of genes based on PC1
boxplot(data[gene_list1_pos,], main = "PC1 - Positive contributions")
boxplot(data[gene_list1_neg,], main = "PC1 - Negative contributions")

# Selection of genes based on PC2
boxplot(data[gene_list2_pos,], main = "PC2 - Positive contributions")
boxplot(data[gene_list2_neg,], main = "PC2 - Negative contributions")

# Selection of genes based on PC3
boxplot(data[gene_list3_pos,], main = "PC3 - Positive contributions")
boxplot(data[gene_list3_neg,], main = "PC3 - Negative contributions")

