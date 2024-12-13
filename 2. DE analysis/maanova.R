#####################################
######## MAANOVA ####################
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


# Choose comparisons to be done
C=rbind(c(-1,1,0,0,0,0,0,0,0,0),       # comparison T0 vs Tm24
        c(0,-1,1,0,0,0,0,0,0,0),       # comparison T6 vs T0
        c(0,0,-1,1,0,0,0,0,0,0),       # comparison T12 vs T6
        c(0,0,0,-1,1,0,0,0,0,0),       # comparison T18 vs T12
        c(0,0,0,0,-1,1,0,0,0,0),       # comparison T24 vs T18
        c(0,0,0,0,0,-1,1,0,0,0),       # comparison T30 vs T24
        c(0,0,0,0,0,0,-1,1,0,0),       # comparison T42 vs T30
        c(0,0,0,0,0,0,0,-1,1,0),       # comparison T54 vs T42
        c(0,0,0,0,0,0,0,0,-1,1)        # comparison T96 vs T54
        )

test.fix=matest(data,fit.fix,Contrast=C,term="Sample",n.perm=300)

# Correction of p-values
test.fix=adjPval(test.fix,method="stepup")


result=data.frame(ID=unique(as.character(rd$ID)),Name=unique(as.character(rd$Name)))


# comparison T0 vs Tm24
comp=c("T0","Tm24")
LR=fit.fix$Sample[,2]-fit.fix$Sample[,1]
FC=fold(LR)
pval=test.fix$Fs$Pvalperm[,1]
adjpval=test.fix$Fs$adjPvalperm[,1]
Amean=calculeAmean(Avalues,comp)

df.result=data.frame(result,T0vsTm24=LR,FC.T0vsTm24=FC,pvalue.T0vsTm24=pval,adjpvalue.T0vsTm24=adjpval,Amean.T0vsTm24=Amean)
write.table(df.result,file="T0vsTm24.txt",sep="\t",row.names=F)


# comparison T6 vs T0
comp=c("T6","T0")
LR=fit.fix$Sample[,3]-fit.fix$Sample[,2]
FC=fold(LR)
pval=test.fix$Fs$Pvalperm[,2]
adjpval=test.fix$Fs$adjPvalperm[,2]
Amean=calculeAmean(Avalues,comp)

df.result=data.frame(result,T6vsT0=LR,FC.T6vsT0=FC,pvalue.T6vsT0=pval,adjpvalue.T6vsT0=adjpval,Amean.T6vsT0=Amean)
write.table(df.result,file="T6vsT0.txt",sep="\t",row.names=F)


# comparison T12 vs T6
comp=c("T12","T6")
LR=fit.fix$Sample[,4]-fit.fix$Sample[,3]
FC=fold(LR)
pval=test.fix$Fs$Pvalperm[,3]
adjpval=test.fix$Fs$adjPvalperm[,3]
Amean=calculeAmean(Avalues,comp)

df.result=data.frame(result,T12vsT6=LR,FC.T12vsT6=FC,pvalue.T12vsT6=pval,adjpvalue.T12vsT6=adjpval,Amean.T12vsT6=Amean)
write.table(df.result,file="T12vsT6.txt",sep="\t",row.names=F)


# comparison T18 vs T12
comp=c("T18","T12")
LR=fit.fix$Sample[,5]-fit.fix$Sample[,4]
FC=fold(LR)
pval=test.fix$Fs$Pvalperm[,4]
adjpval=test.fix$Fs$adjPvalperm[,4]
Amean=calculeAmean(Avalues,comp)

df.result=data.frame(result,T18vsT12=LR,FC.T18vsT12=FC,pvalue.T18vsT12=pval,adjpvalue.T18vsT12=adjpval,Amean.T18vsT12=Amean)
write.table(df.result,file="T18vsT12.txt",sep="\t",row.names=F)


# comparison T24 vs T18
comp=c("T24","T18")
LR=fit.fix$Sample[,6]-fit.fix$Sample[,5]
FC=fold(LR)
pval=test.fix$Fs$Pvalperm[,5]
adjpval=test.fix$Fs$adjPvalperm[,5]
Amean=calculeAmean(Avalues,comp)

df.result=data.frame(result,T24vsT18=LR,FC.T24vsT18=FC,pvalue.T24vsT18=pval,adjpvalue.T24vsT18=adjpval,Amean.T24vsT18=Amean)
write.table(df.result,file="T24vsT18.txt",sep="\t",row.names=F)


# comparison T30 vs T24
comp=c("T30","T24")
LR=fit.fix$Sample[,7]-fit.fix$Sample[,6]
FC=fold(LR)
pval=test.fix$Fs$Pvalperm[,6]
adjpval=test.fix$Fs$adjPvalperm[,6]
Amean=calculeAmean(Avalues,comp)

df.result=data.frame(result,T30vsT24=LR,FC.T30vsT24=FC,pvalue.T30vsT24=pval,adjpvalue.T30vsT24=adjpval,Amean.T30vsT24=Amean)
write.table(df.result,file="T30vsT24.txt",sep="\t",row.names=F)


# comparison T42 vs T30
comp=c("T42","T30")
LR=fit.fix$Sample[,8]-fit.fix$Sample[,7]
FC=fold(LR)
pval=test.fix$Fs$Pvalperm[,7]
adjpval=test.fix$Fs$adjPvalperm[,7]
Amean=calculeAmean(Avalues,comp)

df.result=data.frame(result,T42vsT30=LR,FC.T42vsT30=FC,pvalue.T42vsT30=pval,adjpvalue.T42vsT30=adjpval,Amean.T42vsT30=Amean)
write.table(df.result,file="T42vsT30.txt",sep="\t",row.names=F)


# comparison T54 vs T42
comp=c("T54","T42")
LR=fit.fix$Sample[,9]-fit.fix$Sample[,8]
FC=fold(LR)
pval=test.fix$Fs$Pvalperm[,8]
adjpval=test.fix$Fs$adjPvalperm[,8]
Amean=calculeAmean(Avalues,comp)

df.result=data.frame(result,T54vsT42=LR,FC.T54vsT42=FC,pvalue.T54vsT42=pval,adjpvalue.T54vsT42=adjpval,Amean.T54vsT42=Amean)
write.table(df.result,file="T54vsT42.txt",sep="\t",row.names=F)


# comparison T96 vs T54
comp=c("T96","T54")
LR=fit.fix$Sample[,10]-fit.fix$Sample[,9]
FC=fold(LR)
pval=test.fix$Fs$Pvalperm[,9]
adjpval=test.fix$Fs$adjPvalperm[,9]
Amean=calculeAmean(Avalues,comp)

df.result=data.frame(result,T96vsT54=LR,FC.T96vsT54=FC,pvalue.T96vsT54=pval,adjpvalue.T96vsT54=adjpval,Amean.T96vsT54=Amean)
write.table(df.result,file="T96vsT54.txt",sep="\t",row.names=F)
