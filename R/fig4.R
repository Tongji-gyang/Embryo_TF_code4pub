library(data.table)
library(dplyr)
library(ggpubr)
library(pheatmap)
library(stringr)
library(tidyr)
setwd("~/work_space/4.ProjectET/analysis/fig4/1.epi_influence_binding")

######################################
########  Stage I: Stage-specific binding of TFAP2C
######################################

tt_signal <- data.frame(fread("all_stagePeaks_T1_signal.tab"),check.names = F)
tt_anno <- data.frame(fread("~/work_space/4.ProjectET/analysis/fig4/build/tf_peaks/3.zeroone/Peak.anno"),check.names = F)
tt_status <- data.frame(fread("~/work_space/4.ProjectET/analysis/fig4/build/tf_peaks/3.zeroone/merge_zeroone.matrix"),check.names = F)
colnames(tt_anno) <- c("chr","start","end","peak_index")
tt_status$peak_index <- tt_anno$peak_index
tt_status$status <- str_c(tt_status$`8Cell_zeroone.bed`,tt_status$EPI_zeroone.bed,tt_status$EXE_zeroone.bed)
tt_signal <- merge(tt_signal,tt_anno, by=c("chr","start","end"))

tt_status <- tt_status[order(tt_status$`8Cell_zeroone.bed`,tt_status$EPI_zeroone.bed,tt_status$EXE_zeroone.bed, decreasing = T), ]
rownames(tt_signal) <- tt_signal$peak_index
tt_signal <- tt_signal[tt_status$peak_index,]
tt_signal <- merge(tt_signal,tt_status,by=c("peak_index"))
rownames(tt_signal) <- tt_signal$peak_index
#tt_signal[,4:10] <- log2(tt_signal[,4:10] + 1)
tt_signal[,5:10] <- scale(log2(tt_signal[,5:10] + 1))
tt_status_select <- tt_status[which(tt_status$status %in% c("100","010","001")),]
#c(111,110,101,100,011,010,001
plt.tmp <- tt_signal[tt_status_select$peak_index,5:10]
plt.tmp <- plt.tmp[,c(1,5,4,2,6,3)]

anno_row <- data.frame(status=tt_status_select$status)
rownames(anno_row) <- tt_status_select$peak_index
p <- pheatmap(t(scale(t(plt.tmp))),cluster_rows = F,cluster_cols = F,show_rownames = F,
         color = colorRampPalette(colors = c("NavyBlue","black","Yellow1"))(100),
         annotation_row = anno_row)
ggsave("T1_binding_110_010_001.pdf",p$gtable,width = 4,height = 3)

# temp <- tt_signal[tt_status_select$peak_index,1:10]
fwrite(temp[,2:4], "T1binding_8c_epi_exe_100_010_001.bed",quote = F,sep = "\t",col.names = F)
fwrite(tt_signal[which(tt_signal$status == "100"),c("chr","start","end")], "T1binding_8c_epi_exe_100.bed",quote = F,sep = "\t",col.names = F)
fwrite(tt_signal[which(tt_signal$status == "010"),c("chr","start","end")], "T1binding_8c_epi_exe_010.bed",quote = F,sep = "\t",col.names = F)
fwrite(tt_signal[which(tt_signal$status == "001"),c("chr","start","end")], "T1binding_8c_epi_exe_001.bed",quote = F,sep = "\t",col.names = F)
fwrite(tt_signal[which(tt_signal$status == "111"),c("chr","start","end")], "T1binding_8c_epi_exe_111.bed",quote = F,sep = "\t",col.names = F)


# k9_signal <- data.frame(fread("T1_K9me3_signal.tab"),check.names = F)
# k9_signal <- merge(k9_signal, tt_anno, by=c("chr","start","end"))
# rownames(k9_signal) <- k9_signal$peak_index
# k9_signal[,4:9] <- scale(log2(k9_signal[,4:9] + 1))
# k9_signal <- k9_signal[tt_status_select$peak_index,]
# 
# k9_signal[,4:9] <- t(scale(t(k9_signal[,4:9])))
# k9_signal <- merge(k9_signal,tt_status,by="peak_index")
# 
# pheatmap(t(scale(t(k9_signal[,4:9]))),cluster_rows = F,cluster_cols = F,show_rownames = F,
#          color = colorRampPalette(colors = c("NavyBlue","black","Yellow1"))(100),
#          annotation_row = anno_row)

######################################
########  Stage II: GO of Stage-specific peaks
######################################
library(stringr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(ChIPseeker)
library(ggplot2)
txdb<- TxDb.Mmusculus.UCSC.mm10.knownGene

bed2genename <- function(bedfile){
  peakAnno <- annotatePeak(bedfile, tssRegion=c(-2000, 500),
                           TxDb=txdb, annoDb="org.Mm.eg.db")
  peakAnno <- data.frame(peakAnno)
  peakAnno <- peakAnno[which(peakAnno$annotation %in% c('Promoter (<=1kb)','Promoter (1-2kb)')),]
  genename <- sort(unique(peakAnno$SYMBOL))
  print(str_c("The gene list length is:",length(genename)))
  return(genename)
}



#gran_theme <- theme_classic() + theme(legend.title =element_blank(),legend.text = element_text(size = 12, face = "bold"),axis.title =element_text(size=12,face = "bold") ,axis.text=element_text(size=12))

setwd("~/work_space/4.ProjectET/analysis/fig4/1.epi_influence_binding")

test <- bed2genename("T1binding_8c_epi_exe_100.bed")
ego2 <- enrichGO(gene         = test,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
dotplot(ego2, showCategory=10) + ggtitle("GO for 8C-specific TFAP2C peaks") 
ggsave("Dotplot_GO_8C-specific.pdf",width=7 ,height = 4)

test <- bed2genename("T1binding_8c_epi_exe_010.bed")
ego2 <- enrichGO(gene         = test,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
dotplot(ego2, showCategory=10) + ggtitle("GO for EPI-specific TFAP2C peaks") 
ggsave("Dotplot_GO_EPI-specific.pdf",width=7 ,height = 4)

test <- bed2genename("T1binding_8c_epi_exe_001.bed")
ego2 <- enrichGO(gene         = test,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
dotplot(ego2, showCategory=10) + ggtitle("GO for EXE-specific TFAP2C peaks") 
ggsave("Dotplot_GO_EXE-specific.pdf",width=7 ,height = 4)

############## Stage-peaks ---> gene expr?
gexpr1 <-  data.frame(fread("~/work_space/Bio-Resource/embryos/RNA/Xiewei_2016_Nature_log2FPKM_RNA.csv"),check.names = F)
gexpr2 <- data.frame(fread("~/work_space/Public_data/q.E65_epiexe_xw_zy_RNA/Embryo_E65_RNA_yizhang_xw_collapse.csv"),check.names = F)
gexpr3 <- data.frame(fread("~/work_space/Public_data/q.E65_epiexe_xw_zy_RNA/GSE125318_XW_WT.RNA.FPKM.txt"),check.names = F)
gexpr3[,2:5] <- log2(gexpr3[,2:5] + 1)
colnames(gexpr3) <- c("Gene","E65Epi","E75Ect","E85Head","E95Head")
gexpr <- merge(gexpr1,gexpr2, by.x="Gene", by.y="name")
gexpr_bar <- gexpr[,c("Gene","8.cell","ICM","epi")]
gexpr_bar <- merge(gexpr_bar,gexpr3,by="Gene")
gexpr_bar <- gexpr_bar[c(1,2,3,5,6,7,8)]


####### gene expr 100 010 001
test <- bed2genename("T1binding_8c_epi_exe_100.bed")
temp <- gexpr[which(gexpr$Gene %in% test),2:4] # 974, merged 908
colnames(temp) <- c("8Cell","EPI","EXE")
wilcox.test(temp$`8Cell`,temp$EPI) #9.087e-10
wilcox.test(temp$`8Cell`,temp$EXE) #1.427e-12
temp <- gather(temp,"sample","expr")
p <- ggboxplot(temp, x = "sample", y = "expr",
               color = "black", fill="white",
               width = 0.7,notch = TRUE, outlier.shape = NA) + stat_boxplot(geom="errorbar",width=0.6) 

p <- ggpar(p,legend = "none",xlab = FALSE,ylim=c(0,12),x.text.angle = 45)
p
ggsave("RNA_100.pdf",dpi=300,width=3,height=4)

test <- bed2genename("T1binding_8c_epi_exe_001.bed")
temp <- gexpr[which(gexpr$Gene %in% test),2:4] # 341, merged 309
colnames(temp) <- c("8Cell","EPI","EXE")
wilcox.test(temp$`8Cell`,temp$EPI) #0.11
wilcox.test(temp$EXE,temp$EPI) #0.01
temp <- gather(temp,"sample","expr")
p <- ggboxplot(temp, x = "sample", y = "expr",
               color = "black", fill="white",
               width = 0.7,notch = TRUE, outlier.shape = NA) + stat_boxplot(geom="errorbar",width=0.6) 

p <- ggpar(p,legend = "none",xlab = FALSE,ylim=c(0,12),x.text.angle = 45)
p
ggsave("RNA_010.pdf",dpi=300,width=3,height=4)

test <- bed2genename("T1binding_8c_epi_exe_010.bed")
temp <- gexpr[which(gexpr$Gene %in% test),2:4] # 5611, merged 5186
colnames(temp) <- c("8Cell","EPI","EXE")
wilcox.test(temp$`8Cell`,temp$EXE) #< 2.2e-16
wilcox.test(temp$EPI,temp$EXE) #0.002
temp <- gather(temp,"sample","expr")
p <- ggboxplot(temp, x = "sample", y = "expr",
               color = "black", fill="white",
               width = 0.7,notch = TRUE, outlier.shape = NA) + stat_boxplot(geom="errorbar",width=0.6) 

p <- ggpar(p,legend = "none",xlab = FALSE,ylim=c(0,12),x.text.angle = 45)
p
ggsave("RNA_001.pdf",dpi=300,width=3,height=4)

####### EPI-specific, prime expr?

test <- bed2genename("T1binding_8c_epi_exe_010.bed")
temp <- gexpr3[which(gexpr3$Gene %in% test),] # 341, merged 291
colnames(temp) <- c("sample","E65Epi","E75Ect","E85Head","E95Head")
wilcox.test(temp$E65Epi,temp$E75Ect) #0.0005
wilcox.test(temp$E85Head,temp$E75Ect) #1.238e-05
wilcox.test(temp$E85Head,temp$E95Head) #0.16
temp <- gather(temp,"sample","expr")

p <- ggboxplot(temp, x = "sample", y = "expr",
               color = "black", fill="white",
               width = 0.7,notch = TRUE, outlier.shape = NA) + stat_boxplot(geom="errorbar",width=0.6)

p <- ggpar(p,legend = "none",xlab = FALSE,ylim=c(0,4),x.text.angle = 45)
p
ggsave("RNA_010_post_expression.pdf",dpi=300,width=3,height=4)

####### EPI-specific VS Super bivalency genes 

library(VennDiagram)
twosetfisher <- function(set1,set2,total_num,venn_name,venn_main=""){
  
  temp_int <- intersect(set1,set2)
  temp_merge <- unique(sort(c(set1,set2)))
  temp_mat <-
    matrix(c(length(temp_int),  length(setdiff(set1,temp_int)),length(setdiff(set2,temp_int)), total_num-length(temp_merge)),
           nrow = 2,
           dimnames = list(Guess = c("in", "out"),
                           Truth = c("in", "out")))
  print(fisher.test(temp_mat))
  
  pdf(venn_name)
  
  D1<-venn.diagram(list(set1=set1,set2=set2),filename=NULL,lwd=1,lty=1,col=c('red','blue'),
                   fill=c('red','blue'),cat.col=c('red','blue'),rotation.degree=180,main = venn_main)
  grid.draw(D1)
  dev.off()
  
}

sb_genes <- data.frame(fread("/home1/gyang/work_space/Bio-Resource/embryos/geneset/E65EPI_super_bivalency.txt",header = F),check.names = F)$V1
test <- bed2genename("T1binding_8c_epi_exe_010.bed")
twosetfisher(sb_genes,test,24225,"EPI_specific_SuperBivalency_Venn.pdf","(< 2.2e-16)")


######################################
########  Stage III: Verfication of the heterogeneity (Spatial–temporal expression)
######################################

test <- bed2genename("T1binding_8c_epi_exe_010.bed")
temp <- gexpr3[which(gexpr3$Gene %in% test),] # 341, merged 291
colnames(temp) <- c("sample","E65Epi","E75Ect","E85Head","E95Head")
temp$cf = temp$E75Ect-temp$E65Epi

gexpr_bar
plt.tmp <- data.frame(time=colnames(gexpr_bar)[-1],value=t(gexpr_bar[which(gexpr_bar$Gene == "Pax6"),][,2:7])[,1])
ggbarplot(plt.tmp, x = "time", y = "value",fill = "#F26A6B")
ggsave("Pax6_expr.pdf",dpi=300,width=3,height=4)
plt.tmp <- data.frame(time=colnames(gexpr_bar)[-1],value=t(gexpr_bar[which(gexpr_bar$Gene == "Ccnd2"),][,2:7])[,1])
ggbarplot(plt.tmp, x = "time", y = "value",fill = "#F26A6B")
ggsave("Ccnd2_expr.pdf",dpi=300,width=3,height=4)
plt.tmp <- data.frame(time=colnames(gexpr_bar)[-1],value=t(gexpr_bar[which(gexpr_bar$Gene == "Lhx5"),][,2:7])[,1])
ggbarplot(plt.tmp, x = "time", y = "value",fill = "#F26A6B")
ggsave("Lhx5_expr.pdf",dpi=300,width=3,height=4)
plt.tmp <- data.frame(time=colnames(gexpr_bar)[-1],value=t(gexpr_bar[which(gexpr_bar$Gene == "Hoxa2"),][,2:7])[,1])
ggbarplot(plt.tmp, x = "time", y = "value",fill = "#F26A6B")
ggsave("Hoxa2_expr.pdf",dpi=300,width=3,height=4)
plt.tmp <- data.frame(time=colnames(gexpr_bar)[-1],value=t(gexpr_bar[which(gexpr_bar$Gene == "Gbx2"),][,2:7])[,1])
ggbarplot(plt.tmp, x = "time", y = "value",fill = "#F26A6B")
ggsave("Gbx2_expr.pdf",dpi=300,width=3,height=4)
