#####################################
######## SEARCH FOR DE GENES ########
#####################################

library(maanova)
load("Avalues.RData")

datafile="datafile.txt"
designfile="designfile.txt"

data=read.madata(datafile=datafile,designfile=designfile,arrayType="twoColor",spotflag=T,n.rep=4,avgreps=0,log.trans=T,
  probeid=1, metarow=3, metacol=4, row=5, col=6, intensity=7)
rd=read.table(file=datafile,header=T,sep="\t")

fit.fix=fitmaanova(data,formula=~Sample+Cinetique)

##### Differential analysis #######

fdr=function (p, method = c("stepup", "adaptive", "stepdown", "jsFDR"))
{
    method <- match.arg(method)
    if (method == "stepup" | method == "adaptive" | method ==
        "setdown") {
        m <- length(p)
        tmp <- sort(p, index.return = TRUE)
        sortp <- tmp$x
        idx <- tmp$ix
        if (method == "stepdown") {
            d <- m:1
            sortp <- (1 - (1 - sortp)^d) * d/m
            for (i in 1:(m - 1)) {
                if (sortp[i + 1] < sortp[i])
                  sortp[i + 1] <- sortp[i]
            }
        }
        else {
            if (method == "stepup")
                m0 <- m
            else if (method == "adaptive") {
                s <- sort(1 - sortp)/(1:m)
                m0raw <- m
                i <- m
                while (i > 1 && s[i] <= s[i - 1]) i <- i - 1
                if (i > 1)
                  m0raw <- 1/s[i - 1]
                else m0raw <- 1/s[1]
                m0 <- min(floor(1 + m0raw), m)
            }
            sortp <- sortp * m0/(1:m)
            for (i in (m - 1):1) {
                if (sortp[i] > sortp[i + 1])
                  sortp[i] <- sortp[i + 1]
            }
        }
        result <- NULL
        result[idx] <- sortp
    }
    else if (method == "jsFDR") {
        library(qvalue)
        result = qvalue(p)$qvalues
    }
    else stop("Need to specify FDR method (stepup, adaptive, stepdown, or jsFDR). To use jsFDR, one needs to install qvalue() package")
    result
}


calculeAmean <- function(Avalues,comp){
  library(limma)
  ind1=grep(comp[1],colnames(Avalues))
  ind2=grep(comp[2],colnames(Avalues))
  if(length(ind1)==0 | length(ind2)==0)
    stop("Une des conditions n'existe pas")
  A <- as.matrix(Avalues[,c(ind1,ind2)])
  Amean <- rowMeans(unwrapdups(A,ndups=4,spacing=1),na.rm=T)
  return(Amean)
}

fold <- function(x){
  return(sign(x)*(2^abs(x)))
}

# Design matrix to perform all pairwise
# comparisons between time points
C = matrix(0, nrow = 90, ncol = 10)

ligne = 1
for(i in 1:10){
  for(j in 1:10){
    if(i != j){
      C[ligne, i] = -1
      C[ligne, j] = 1
      ligne = ligne + 1
    }
  }
 
}

# Names of time point
SampNames = c("Tm24", "T0", "T6", "T12", "T18", "T24", "T30", "T42", "T54", "T96")

# Statistical analysis
test.fix=matest(data,fit.fix,Contrast=C,term="Sample",n.perm=300)

# Correction of pvalues
test.fix=adjPval(test.fix,method="stepup")

# Get results
result=data.frame(ID=unique(as.character(rd$ID)),Name=unique(as.character(rd$Name)))

# Table to combine all the results
diffRes = NULL

for(i in 1:10){
  for(j in 1:10){
    if(i != j){
      
      comp = c(SampNames[j],SampNames[i])
      LR=fit.fix$Sample[,j]-fit.fix$Sample[,i]
      FC=fold(LR)
      pval=test.fix$Fs$Pvalperm[,i]
      adjpval=test.fix$Fs$adjPvalperm[,i]
      Amean=calculeAmean(Avalues,comp)
      
      df.result=data.frame(result,logFoldChange=LR, FoldChange=FC, pvalue=pval,adjpvalue=adjpval,Amean=Amean)
      write.table(df.result,file=paste0(SampNames[j], "vs", SampNames[i], ".txt"),sep="\t",row.names=F)
      
      temp = cbind(LR, adjpval)
      colnames(temp) = c(paste0("LogFC.",SampNames[j], "vs", SampNames[i]),
                         paste0("adjPval.",SampNames[j], "vs", SampNames[i]))
      diffRes = cbind(diffRes, temp)
      
    }
  }
}

row.names(diffRes) = df.result[,"Name"]

# Results writing
write.table(diffRes, file = "diffRes_allComp.txt",
            sep = "\t", quote = F)