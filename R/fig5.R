setwd('~/work_space/4.ProjectET/analysis/fig5')

# tmp1 <-  data.frame(fread("EPI_only_epi_signal.tab"),check.names = F)
# tmp1[,4:9] <- log2(tmp1[,4:9] + 1)
# tmp1$fc = tmp1$E65Exe_H3K9me3 - tmp1$E65Epi_H3K9me3
# tmp1 <- tmp1[tmp1$E65Epi_H3K27me3 < 4,]
# 
# tmp1 <-  data.frame(fread("EXE_only_epi_signal.tab"),check.names = F)
# tmp1[,4:9] <- log2(tmp1[,4:9] + 1)
# tmp1$fc = tmp1$E65Epi_H3K9me3 - tmp1$E65Exe_H3K9me3
# tmp1 <- tmp1[tmp1$E65Epi_H3K27me3 < 4,]


tmp1 <-  data.frame(fread("SB.tab"),check.names = F)
tmp1[,4:9] <- log2(tmp1[,4:9] + 1)
tmp1 <- tmp1[,4:9]
wilcox.test(tmp1$E65Epi_H3K27me3,tmp1$E65Exe_H3K27me3)
wilcox.test(tmp1$E65Epi_H3K4me3,tmp1$E65Exe_H3K4me3)
wilcox.test(tmp1$E65Epi_H3K9me3,tmp1$E65Exe_H3K9me3)
#tmp1$fc = tmp1$E65Epi_H3K9me3 - tmp1$E65Exe_H3K9me3
#tmp1 <- tmp1[tmp1$E65Epi_H3K27me3 < 4,]
tmp1 <- gather(tmp1,"name","value")
tmp1$sample <-  str_split_fixed(tmp1$name,"_",2)[,1]
tmp1$modi <-  str_split_fixed(tmp1$name,"_",2)[,2]

ggboxplot(tmp1,"sample","value",facet.by = "modi")
ggsave("E65_SB_epi.pdf")
