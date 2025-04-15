# Data reading
orthoData <- as.matrix(read.csv("Orthogroups.tsv", sep = "\t", row.names = 1))
DEgenes <- as.matrix(read.csv("Orthogroups_DE-genes.txt", 
                              sep = "\t", row.names = 1))

# Table to store orthologous genes for DE genes
orthoData2 <- NULL

for(i in 1:nrow(DEgenes)){
  # Get orthogroup number
  OrthoNum <- DEgenes[i, 3]
  DEgene   <- DEgenes[i,2]
  
  if(OrthoNum %in% row.names(orthoData)){
    # Get lists of genes
    orthoData2 <- rbind(orthoData2, c(DEgene, OrthoNum, orthoData[OrthoNum,]))
  }
  
# end of for()
}

# Writing of the results
write.table(orthoData2, file = "DE-genes_in_Orthogroups.txt",
            row.names = F, col.names = T, quote = F, sep = "\t")



