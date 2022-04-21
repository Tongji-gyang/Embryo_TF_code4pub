setwd("~/work_space/4.ProjectET/analysis/figs7")
library(data.table)

T1kk <- data.frame(fread("all_Tfap2c_K9ac_K9me3.tab"),check.names = F)
T1kk <- T1kk[which(T1kk$morula_H3K9me3 > 3),]
T1kk <- T1kk[which(T1kk$`K9-morula` > 5),]
T1kk[,4:11] <- t(apply(T1kk[,4:11],1,order))
T1kk[,12:19] <- t(apply(T1kk[,12:19],1,order))  
T1kk <- T1kk[which(T1kk$morula_H3K9me3 > 3),]
T1kk <- T1kk[which(T1kk$`K9-morula` > 3),]

apply(T1kk[,4:11],1,order)
order(T1kk[,4:11])


colnames(T1kk) <- c("chr","start", "end", "m_H3K9me3", "I_H3K9me3", "T_H3K9me3", "m_H3K9ac", "I_H3K9ac", "T_H3K9ac")
T1kk <- T1kk[which(rowSums(T1kk[,4:6]) > 5),]
T1kk <- T1kk[which(rowSums(T1kk[,7:9]) > 10),]
T1kk_k9me3 <- T1kk[,1:6]
T1kk_k9ac <- T1kk[,c(1:3,7:9)]
T1kk_k9me3[,4:6] <- t(scale(t(T1kk_k9me3[,4:6]))) 
T1kk_k9ac[,4:6] <- t(scale(t(T1kk_k9ac[,4:6]))) 
T1kk <- merge(T1kk_k9me3,T1kk_k9ac, by=c("chr","start", "end"))
pheatmap(T1kk[,4:9],show_rownames = F,kmeans_k=15,cluster_cols = F)



T2kk <- data.frame(fread("morula_Tead4_K9ac_K9me3.tab"),check.names = F)
colnames(T1kk) <- c("chr","start", "end", "H3K9me3", "H3K9ac")
