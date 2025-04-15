# data reading
files <- c("WaveI.txt",
           "WaveII.txt",
           "WaveIII.txt",
           "WaveIV.txt",
           "WaveV.txt")

# GO terms 
dataMetaGO = as.matrix(read.table("dataMETAGO_2.txt", header = T))

# Description of GO terms
GOinfo = as.matrix(read.csv("Level_GO.txt", header = T, sep = "\t"))
# Level to analyse GO terms;
L = 3
Col = which(colnames(dataMetaGO) == paste0("METAGO_L", L))

# New table to store the results
dataMetaGO2 = NULL

# Extract GO terms for each genes
for(g in 1:nrow(dataMetaGO)){
  
  print(paste(g, "/", nrow(dataMetaGO)))
  metaFun = unlist(strsplit(dataMetaGO[g,Col], "|", fixed = T))  
  
  if(length(metaFun) > 1){
    for(m in metaFun){
      newRow = dataMetaGO[g,]
      newRow[Col] = m
      
      dataMetaGO2 = rbind(dataMetaGO2, newRow)
    }
  }
  
  if(length(metaFun) == 1){
    dataMetaGO2 = rbind(dataMetaGO2, dataMetaGO[g,])
  }
}

# Analyze list of genes in each wave
for(f in files){
  
  expData  = as.matrix(read.table(paste0("./",f), header = T, row.names = 1))
  geneList = row.names(expData)
  subDataMetaGO2 = dataMetaGO2[dataMetaGO2[,"gene_id"] %in% geneList,] 

  resFun = table(subDataMetaGO2[,Col])
  
  # We use the following titorial :
  # http://pedagogix-tagc.univ-mrs.fr/courses/ASG1/practicals/go_statistics_td/go_statistics_td_2015.html

  allFunRes = NULL
  for(i in 1:length(resFun)){

    print(paste("Statistical evaluation for function", names(resFun)[i]))
  
    n <- sum(geneList %in% dataMetaGO2[,"gene_id"])  
    k <- resFun[i] 
    N <- length(unique(dataMetaGO2[, "gene_id"]))  
    M <- length(unique(dataMetaGO2[dataMetaGO2[,Col] == names(resFun)[i], "gene_id"])) 

    p.value <-  phyper(k-1, M, N-M, n, lower.tail = FALSE)

    # Final results
    funRes = c(names(resFun)[i], 
               paste0(GOinfo[which(GOinfo[, "GO_ID"] == names(resFun)[i]),], collapse = "|"),
               paste0(unique(subDataMetaGO2[subDataMetaGO2[,Col] == names(resFun)[i], "gene_id"]), collapse = "|"), 
               n, k, N, M, p.value)

    allFunRes = rbind(allFunRes, funRes)
  # end of for()
  }

  # P-values are finally adjusted (multiple testing)
  allFunRes2 = cbind(allFunRes, p.adjust(allFunRes[,ncol(allFunRes)], method = "BH", n = nrow(allFunRes)))

  colnames(allFunRes2) = c("GO_id", "GO_info", "List of genes in function", 
                           "n", "k", "N", "M", "P-val", "adj.P-val")

  write.table(allFunRes2, file = paste0(f,"_GO_analysis.res"), 
              sep = "\t", quote = F, row.names = F)
  
# End of for()
}
