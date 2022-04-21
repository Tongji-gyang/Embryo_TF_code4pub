library(data.table)
library(GGally)
library(pheatmap)



refresh_rna <- function(){
  rna <- read.csv("/home1/gyang/work_space/Public_data/1.wt_early_embryo_rna/wt_mean_expression.csv",row.names = 1)
  colnames(rna) <- c("MII","Zygote","2cell","4cell","8cell","Morula","Blastocyst","ICM","TE")
  rna <- rna[,c("MII","Zygote","2cell","4cell","8cell","Morula","ICM","TE")]
  rna$gene <- rownames(rna)
  
  rna_e2c <- read.csv("~/work_space/Public_data/9.Science_e2C_rna/wt_e2C_science_log2mean_expression.csv",col.names = c("gene","e2cell"))
  rna_e2c <- data.frame(tapply(rna_e2c$e2cell,rna_e2c$gene,mean))
  rna_e2c$gene <- rownames(rna_e2c)
  colnames(rna_e2c)[1] <- "e2cell"
  
  rna_m2c <- read.csv("~/work_space/Public_data/9.Science_e2C_rna/m2C/wt_m2C_science_log2mean_expression.csv",col.names = c("gene","m2cell"))
  rna_m2c <- data.frame(tapply(rna_m2c$m2cell,rna_m2c$gene,mean))
  rna_m2c$gene <- rownames(rna_m2c)
  colnames(rna_m2c)[1] <- "m2cell"
  
  rna_es <- data.frame(fread("~/work_space/1.Mouse_Acetylation/buffet/k9_analysis/3.peak/fig2_part_of_peak/b.density/es_fpkm_all_from_mliu.csv")[,c("geneid","J.S_M_10_P3","J.S_M_16_P3","J.S_F_17_P3","J.S_F_46_P3")])
  rna_es$mES <- log2(rowMeans(rna_es[,2:5]) + 1)
  rna_es <- rna_es[,c("geneid","mES")]
  
  rna_ts <- data.frame(fread("~/work_space/Public_data/c.gao_tsc_rna/TSC_rna_expr.txt")[,c("geneid","TSC")])
  
  rna <- merge(rna, rna_e2c, by.x = "gene", by.y = "gene")
  rna <- merge(rna, rna_m2c ,by.x = "gene", by.y = "gene")
  rna <- merge(rna, rna_es ,by.x = "gene", by.y = "geneid")
  rna <- merge(rna, rna_ts ,by.x = "gene", by.y = "geneid")
  colnames(rna)[4] <- "l2cell"
  
  rownames(rna) <- rna$gene
  return(rna)
}

########## 0. CR versus ChIP
setwd("~/work_space/4.ProjectET/analysis/figs1/0.CR_ChIP_TSC_compare/2.promoter")
promoter_T1 <- data.frame(fread("CR_ChIP_TSC_compare_promoter_2203.tab"))
promoter_T1[,4:6] <- log2(promoter_T1[,4:6] + 1)
cor(promoter_T1[,4:6])
TSC_ChIP_Sample1 TSC_ChIP_Sample2 TSC_Tfap2c_1k
# TSC_ChIP_Sample1        1.0000000        0.7997631     0.8038173
# TSC_ChIP_Sample2        0.7997631        1.0000000     0.8955168
# TSC_Tfap2c_1k           0.8038173        0.8955168     1.0000000

smoothScatter(promoter_T1$TSC_Tfap2c_1k,promoter_T1$TSC_ChIP_Sample2, xlim = c(0,8), ylim = c(0,8), xaxs="i", yaxs="i")

png("promoter_T1_CR_ChIP_cor_0.90.png",res = 300,width = 1500,height = 1500)
smoothScatter(promoter_T1$TSC_Tfap2c_1k,promoter_T1$TSC_ChIP_Sample2, xlim = c(0,8), ylim = c(0,8), xaxs="i", yaxs="i", xlab = "TSC_TFAP2C_CUT&RUN", ylab = "TSC_TFAP2C_ChIP-seq")
dev.off()

########## 1.accessibility_motif_selection
setwd("/home1/gyang/work_space/4.ProjectET/analysis/figs1/2.atac_motif")
all_dir <- list.files(pattern = "?peakAnalysis$")
all_dir <- c("early_2-cell_rep1","early_2-cell_rep2","2-cell_rep1","2-cell_rep2","4-cell_rep1","4-cell_rep2","8-cell_rep1","8-cell_rep2","ICM_rep1","ICM_rep2")

all_df <- data.frame()
for(i in all_dir){
  #i="early_2-cell_rep1" 
  km <- data.frame(fread(paste0("./",i,"_peakAnalysis/","knownResults.txt")))[,c(1,3)]
  km$P.value = as.numeric(gsub("1e-","",km$P.value))
  colnames(km) <- c('name','logpval') 
  km$time <- gsub("_peakAnalysis","",i)
  all_df <- rbind(all_df, km)
}

all_df_10plus <- all_df[which(all_df$logpval > 10),]
fwrite(all_df,"../Weixie_ATAC_motif_raw.csv")
fwrite(all_df_10plus,"../Weixie_ATAC_motif_raw_logpabove10.csv")

##################
setwd("/home1/gyang/work_space/4.ProjectET/analysis/figs1/2.dnase_motif")
all_dir <- list.files(pattern = "?peakAnalysis$")
all_dir <- c("one_c_peakAnalysis","two_c_peakAnalysis","four_c_peakAnalysis", "eight_c_peakAnalysis" ,"morula_peakAnalysis")

all_df <- data.frame()
for(i in all_dir){
  #i="early_2-cell_rep1" 
  km <- data.frame(fread(paste0("./",i,"/","knownResults.txt")))[,c(1,3)]
  km$P.value = as.numeric(gsub("1e-","",km$P.value))
  colnames(km) <- c('name','logpval') 
  km$time <- gsub("_peakAnalysis","",i)
  all_df <- rbind(all_df, km)
}

all_df_10plus <- all_df[which(all_df$logpval > 10),]
fwrite(all_df,"../Yizhang_dnase_motif_raw.csv")
fwrite(all_df_10plus,"../Yizhang_dnase_motif_raw_logpabove10.csv")

#################### Yizhang Dnase KLF5 morula
### Homer motif locus calc

setwd("~/work_space/Public_data/k.zhangyi_cell_2016_Dnase_mm9/tf_ana")
raw_bed <- fread("/home1/gyang/work_space/4.ProjectET/analysis/figs1/1.dnase_peak/morula_dhs.peak")
raw_bed$pos <- as.numeric(rownames(raw_bed))
rela_locus <- fread("~/work_space/Public_data/k.zhangyi_cell_2016_Dnase_mm9/tf_ana/morula_c_dhs.KLF5.peak")
bed <- merge(raw_bed,rela_locus, by.x = c("pos"), by.y = c("PositionID"))
bed$motif_fake_start <- bed$V2 + bed$Offset + 1 #fake means suitable for + strand
bed$motif_fake_end <- bed$motif_fake_start + nchar(bed$Sequence) -1
bed$start <- ifelse(bed$Strand == "+", bed$motif_fake_start, 2*bed$motif_fake_start - bed$motif_fake_end)
bed$end <- ifelse(bed$Strand == "+", bed$motif_fake_end, bed$motif_fake_start)
# bed <- merge(bed, promoter, by.x = c("V1","V2"), by.y = c("V1","V2"))[,c("V1","motif_start","motif_end","Sequence","V4")]
# colnames(bed) = c("chr","start","end","seq","gene")
write.table(bed[,c("V1","start","end")],"DNase_morula_KLF5_mm9.bed",sep="\t",quote = F,row.names = F,col.names = F)



embryo_gexpr <-  refresh_rna()
fwrite(embryo_gexpr,"../early_embryo_gexpr.csv")

####### 3.Tfap2c_TSC_cor_heatmap

setwd("~/work_space/4.ProjectET/analysis/figs1")
Tfap2c_TSC <- data.frame(fread("~/work_space/4.ProjectET/analysis/figs1/3.Tfap2c_TSC_cor_heatmap/Tfap2c_TSC_promoter.tab"))
ggcorr(Tfap2c_TSC[, 4:8], hjust = 1, label = TRUE, label_alpha = 1,label_round = 3,
       low = "white", mid = "white", high = "white")
ggsave("Tfap2c_TSC_cor_heatmap_triangle_raw.pdf")

####### 4.Tfap2c_Tead4_cor_heatmap
library(corrplot)
setwd("/home1/gyang/work_space/4.ProjectET/analysis/figs1")
Tfap2c <- data.frame(fread("/home1/gyang/work_space/4.ProjectET/analysis/figs1/4.Tfap2c_Tead4_embyo_cor/all-tfap2c_promoter_cor.txt"),check.names = F)
rownames(Tfap2c) <- Tfap2c$V1
Tfap2c<- Tfap2c[,2:17]
pheatmap(Tfap2c)
corrplot(as.matrix(Tfap2c), type = 'lower', order = 'original', tl.col = 'black', method = 'shade',
         cl.ratio = 0.2, tl.srt = 45,is.corr=F,col.lim = c(0.4,1),addCoef.col = 'black',col = colorRampPalette(colors = c("DodgerBlue1","PaleTurquoise","white","Wheat","RosyBrown1","IndianRed1"))(100),tl.cex=1, addrect = 5) 
#gsave("Tfap2c_samples_cor.pdf",width = 12,height = 12)
Tead4 <- data.frame(fread("/home1/gyang/work_space/4.ProjectET/analysis/figs1/4.Tfap2c_Tead4_embyo_cor/all-tead4_promoter_cor.txt"),check.names = F)
rownames(Tead4) <- Tead4$V1
Tead4<- Tead4[,2:11]
pheatmap(Tead4)
corrplot(as.matrix(Tead4), type = 'upper', order = 'original', tl.col = 'black', method = 'shade',
         cl.ratio = 0.2, tl.srt = 45,is.corr=F,col.lim = c(0.4,1),addCoef.col = 'black',col = colorRampPalette(colors = c("DodgerBlue1","PaleTurquoise","white","Wheat","RosyBrown1","IndianRed1"))(100),tl.cex=1, addrect = 5)
#gsave("Tead4_samples_cor.pdf",width = 12,height = 12)