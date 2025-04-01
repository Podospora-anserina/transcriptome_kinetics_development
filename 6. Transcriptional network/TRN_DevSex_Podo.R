# List of genes which are annotated as coding for transcription factors
TFs <- as.vector(as.matrix(read.table("list_TF_podo.txt")))
                 
# Gene expression profiles for all genes, and only for TF genes
expData <- as.matrix(read.table("data_all_genes.txt", 
                      header = T, sep = "\t"))
expDataTF <- expData[TFs,]

# Threshold to create the transcriptional network
T = 0.95

# List of edges
edgeList <- NULL

# Waves information
waves <- c("I", "II", "III", "IV", "V")
colorW <- c("grey", "red", "orange", "green", "blue")

# Search for putative targets for each TF
for(i in 1:length(TFs)){
  tf <- TFs[i]
  
  for(j in waves){
    expDataW <- as.matrix(read.table(paste0("2021-03-24_Wave", j, ".txt")))
    
    for(k in 1:nrow(expDataW)){
      corVal <- cor(expDataTF[tf,], expDataW[k,])
      if(corVal > T){
        edgeList <- rbind(edgeList, c(tf, row.names(expDataW)[k], corVal, j))
      }
    }
  }
}

edgeList <- edgeList[!(edgeList[,1] == edgeList[,2]),]
colnames(edgeList) <- c("TF", "target", "Cor", "Wave")
sort(table(edgeList[,1]), decreasing = TRUE)

# Get additional information for target genes
geneInfo <- read.csv("gene_infos.csv",
                     header = T, row.names = 1)

tfInfo     <- geneInfo[edgeList[,1], 1:2]
targetInfo <- geneInfo[edgeList[,2], 1:2]

allRes <- cbind(edgeList, tfInfo, targetInfo)
colnames(allRes)[5:8] <- c("tf_name", "tf_description",
                           "target_name", "target_description")

# Export the results
write.table(allRes, file = paste0("TRN_cor", T, ".txt"),
            sep = "\t", quote = F, row.names = F)

# This part of the script allows to create a graphical representation of the network
library(igraph)
TFnet = graph_from_edgelist(edgeList[,1:2], directed = F)
plot.igraph(TFnet, vertex.color = "red", vertex.size = 2, 
            vertex.label = NA,
            edge.color = "grey",
            edge.width = 0.5,
            main = "P. anserina transcriptional network",
            arrow.size = 0.1)
       #     layout=layout_on_sphere) 

nodeNames = get.vertex.attribute(TFnet, name = "name")

TFcolor = "white"
TargetColor = "pink"
vecColor = rep(TargetColor, length(nodeNames))
for(i in 1:length(nodeNames)){
  if(nodeNames[i] %in% edgeList[,2]){
    wave <- unique(edgeList[edgeList[,2] == nodeNames[i],4])
    if(wave == waves[1]){
      vecColor[i] <- colorW[1]
    }
    if(wave == waves[2]){
      vecColor[i] <- colorW[2]
    }
    if(wave == waves[3]){
      vecColor[i] <- colorW[3]
    }
    if(wave == waves[4]){
      vecColor[i] <- colorW[4]
    }
    if(wave == waves[5]){
      vecColor[i] <- colorW[5]
    }
  }
}
vecColor[nodeNames %in% edgeList[,1]] = TFcolor

plot.igraph(TFnet, vertex.color = vecColor, vertex.size = 5, 
            vertex.label = NA,
            arrow.size = 0.1, edge.color = "grey", 
            main = "P. anserina - Network")

TFSize = 10
TargetSize = 5
vecSize = rep(TargetSize, length(nodeNames))
vecSize[nodeNames %in% edgeList[,1]] = TFSize

plot.igraph(TFnet, vertex.color = vecColor, vertex.size = vecSize, vertex.label = NA,
            arrow.size = 0.1, edge.color = "grey", main = "P. anserina - Network")
