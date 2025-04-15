# Data reading :
pfam_podo <- read.csv("data_Pfam_devSex.txt", 
                        header = T, sep = "\t")

res_fam = table(pfam_podo[, "accession"])

# http://pedagogix-tagc.univ-mrs.fr/courses/ASG1/practicals/go_statistics_td/go_statistics_td_2015.html

N <- nrow(pfam_podo) 
final_res <- NULL

for(w in c("I", "II", "III", "IV", "V")){
  
  pfam_podo_wave <- pfam_podo[pfam_podo[, "Wave"] == paste("Wave", w),]
  
  pfam_wave <- table(pfam_podo_wave[, "accession"])
  n <- nrow(pfam_podo_wave)
  
  for(i in 1:length(pfam_wave)){
    k <- pfam_wave[i] 
    M <- sum(pfam_podo[,"accession"] == names(pfam_wave[i]))
    p.value <-  phyper(k-1, M, N-M, n, lower.tail = FALSE)
    final_res <- rbind(final_res, 
                       c(w, names(pfam_wave[i]), k, M, p.value))
  }
}

colnames(final_res) <- c("Wave", "Pfam_access", "k", "M", "pval")

final_res2 = cbind(final_res, p.adjust(final_res[,"pval"], 
                                       method = "BH", 
                                       n = nrow(final_res)))

write.table(final_res2, 
            file = "analysis_Pfam.res", 
            sep = "\t", quote = F, row.names = F)
