library(data.table)
library(dplyr)
library(ggpubr)
library(pheatmap)
library(stringr)
library(tidyr)
setwd("~/work_space/4.ProjectET/analysis/fig1")

####################################################
########  Stage Zero: Expr: Tfap2c/Tead4
####################################################

gexpr1 <-  data.frame(fread("~/work_space/Bio-Resource/embryos/RNA/Xiewei_2016_Nature_log2FPKM_RNA.csv"),check.names = F)
gexpr2 <- data.frame(fread("~/work_space/Public_data/q.E65_epiexe_xw_zy_RNA/Embryo_E65_RNA_yizhang_xw_collapse.csv"),check.names = F)
gexpr4 <- data.frame(fread("~/work_space/Bio-Resource/embryos/RNA/Normal_embryo_noem2c_log2FPKM.csv"),check.names = F)[,c("gene","Morula","TE","TSC" )]
gexpr <- merge(gexpr1,gexpr2, by.x="Gene", by.y="name")
gexpr <- merge(gexpr,gexpr4, by.x="Gene", by.y="gene")
colnames(gexpr) <- c("Gene" ,"D10.oocyte","D14.oocyte","8w.oocyte" ,"MII.oocyte" ,"PN5.zygote","Early.2.cell","Late.2.cell" ,
                     "4.cell", "8.cell" ,"ICM" ,"ESC" , "EPI"  ,"EXE" ,"morula" ,"TE" ,"TSC" )
gexpr <- gexpr[,c("Gene" ,"D10.oocyte","D14.oocyte","8w.oocyte" ,"MII.oocyte" ,"PN5.zygote","Early.2.cell","Late.2.cell" ,
                  "4.cell", "8.cell" ,"morula" ,"ICM"  , "EPI" ,"ESC"  ,"TE" ,"EXE" ,"TSC" )]
tt_expr <- gexpr[gexpr$Gene %in% c("Tfap2c","Tead4"),]
tt_expr <- gather(tt_expr,"time","FPKM",-"Gene")
tt_expr$Gene <- factor(tt_expr$Gene,levels=c("Tfap2c","Tead4"))

p <- ggbarplot(tt_expr, "time", "FPKM",
          fill = "Gene", color = "black", palette = c("#F9C0C0","#F37F72"),
          position = position_dodge(0.6),width=0.5)
ggpar(p,x.text.angle = 45)
ggsave("TT_expr.pdf",width = 8,height = 4)

####################################################
########  Stage I: Dynamic Change: Tfap2c/Tead4
####################################################

######### Total number

t1_num <- data.frame(fread("Tfap2c_rpkm15_peaknum.stat"),check.names = F)
t1_num$Time <- gsub("_Tfap2c","",t1_num$Time)
t1_num$time <- 1:8

pdf("Tfap2c_peaknum.pdf",width = 5,height = 5)
twoord.plot(lx = t1_num$time, ly = t1_num$Total, rx = t1_num$time,
            ry = t1_num$pro_ratio, type=c('bar','line'),
            lcol = 'blue', rcol = 'red',
            ylab = 'Peak number', rylab = 'Ratio of proximal peaks', main = 'Tfap2c peak dynamics',
            xtickpos=t1_num$time, xticklab=t1_num$Time,lylim=c(0,30000),rylim=c(0,100))
dev.off()

t2_num <- data.frame(fread("Tead4_rpkm15_peaknum.stat"),check.names = F)
t2_num$Time <- gsub("_Tead4","",t2_num$Time)
t2_num$time <- 1:6

pdf("Tead4_peaknum.pdf",width = 8,height = 8)
twoord.plot(lx = t2_num$time, ly = t2_num$Total, rx = t2_num$time,
            ry = t2_num$pro_ratio, type=c('bar','line'),
            lcol = 'blue', rcol = 'red',
            ylab = 'Peak number', rylab = 'Ratio of proximal peaks', main = 'Tead4 peak dynamics',
            xtickpos=t2_num$ime, xticklab=t2_num$Time,lylim=c(0,18000),rylim=c(0,100))
dev.off()

######### Gain and lost

setwd("~/work_space/4.ProjectET/analysis/fig1/1.dynamic_change/gain_loss")
t1_gl <- data.frame(fread("./tfap2c/3.zeroone/merge_zeroone.matrix"),check.names = F)
colnames(t1_gl) <- gsub("_zeroone.bed","",colnames(t1_gl))
#"8Cell"  "Morula" "ICM"    "EPI"    "TE"     "EXE"
table(t1_gl$Morula-t1_gl$`8Cell`) 
# -1     0     1 
# 4703 38353 10521 
table(t1_gl$ICM-t1_gl$Morula)
# -1     0     1 
# 9586 37710  6281 
table(t1_gl$TE-t1_gl$Morula)
# -1     0     1 
# 8823 34751 10003 
table(t1_gl$EPI-t1_gl$ICM)
# -1     0     1 
# 16463 36112  1002
table(t1_gl$EXE-t1_gl$TE)
# -1     0     1 
# 9786 33153 10638 

t2_gl <- data.frame(fread("./tead4/3.zeroone/merge_zeroone.matrix"),check.names = F)
colnames(t2_gl) <- gsub("_zeroone.bed","",colnames(t2_gl))
#colnames(t2_gl) "8Cell"  "Morula" "TE"     "EXE" 
table(t2_gl$Morula-t2_gl$`8Cell`) 
# -1     0     1 
# 7850 14562  5484 
table(t2_gl$TE-t2_gl$Morula)
# -1     0     1 
# 10561 16124  1211 
table(t2_gl$EXE-t2_gl$TE)
# -1     0     1 
# 1997 23128  2771 

gl_sum <- data.frame(fread("Gain_loss.csv"),check.names = F)
colnames(gl_sum)[1] <- "time"
#   time  Gain  Loss
# 1 8M_T1 10521  4703
# 2 MI_T1  6281  9586
# 3 MT_T1 10003  8823
# 4 IE_T1  1002 16463
# 5 TE_T1 10638  9786
# 6 8M_T2  5484  7850
# 7 MT_T2  1211 10561
# 8 TE_T2  2771  1997
gl_sum <- gather(gl_sum,"class","num",-"time") 
gl_sum$time <- factor(gl_sum$time,levels = c("8M_T1", "MI_T1" ,"MT_T1" ,"IE_T1" ,"TE_T1", "8M_T2", "MT_T2" ,"TE_T2"))
library(ggpubr)
p <- ggbarplot(gl_sum, "time", "num",
               fill = "class", color = "class", palette =  c("IndianRed1","SteelBlue"),
               lab.pos = "in",width=0.5) 
ggpar(p,x.text.angle=30)
ggsave("../../TT_gainloss_barplot.pdf",width =5,height=4)

####################################################
########  Stage II: TT co-binding and expr
####################################################

############# 1. TT co-binding dynamics
setwd("~/work_space/4.ProjectET/analysis/fig1/2.TF_binding_gexpr")
tt_pro <- data.frame(fread("TT_pro_summary.txt"),check.names = F)
colnames(tt_pro)[1] = "stage"
tmp=tt_pro$stage
tt_pro <-  gather(tt_pro,"class","num",-"stage")
tt_pro$stage <- factor(tt_pro$stage,levels = tmp)
tt_pro$class <- factor(tt_pro$class,levels = c("Tfap2c_only","TT_both","Tead4_only"))
p <- ggbarplot(tt_pro, "stage", "num",
               fill = "class", color = "white", palette =  c("#E585B6","#40E0D0","#E0AD69"),
               lab.pos = "in",width=0.5,alpha=0.8) 
ggpar(p,x.text.angle=45)
ggsave("../TT_co-binding_barplot.pdf",width =4,height=5)


############# 2. TT co-binding ~ gexpr
gexpr <- data.frame(fread("~/work_space/4.ProjectET/analysis/fig1/build/0.whole_transcriptome/Normal_embryo_log2FPKM.csv"),check.names = F)
gexpr <- gexpr[,c("gene","8cell")]
t1_only <- data.frame(fread("~/work_space/4.ProjectET/analysis/fig1/2.TF_binding_gexpr/T1_only_8Cell.promoter",header = F),check.names = F)
t1_only$class <- "Tfap2c_only"
t2_only <- data.frame(fread("~/work_space/4.ProjectET/analysis/fig1/2.TF_binding_gexpr/T2_only_8Cell.promoter",header = F),check.names = F)
t2_only$class <- "Tead4_only"
tt_both <- data.frame(fread("~/work_space/4.ProjectET/analysis/fig1/2.TF_binding_gexpr/T1T2_8Cell.promoter",header = F),check.names = F)
tt_both$class <- "TT_both"
tt_details <- rbind(t1_only,t2_only,tt_both)



gexpr_tt <- merge(gexpr,tt_details,by.x="gene",by.y="V1")
p <- ggboxplot(gexpr_tt, "class", "8cell", color = "black", 
               width=0.5,alpha=0.8) 
ggpar(p,x.text.angle=45)
ggsave("../TT_co-binding_8c_gexpr_8c_barplot.pdf",width =3,height=3)
# wilcox.test(gexpr_tt[which(gexpr_tt$class == "Tead4_only"),"8cell"],gexpr_tt[which(gexpr_tt$class == "Tfap2c_only"),"8cell"])
# 
# Wilcoxon rank sum test with continuity correction
# 
# data:  gexpr_tt[which(gexpr_tt$class == "Tead4_only"), "8cell"] and gexpr_tt[which(gexpr_tt$class == "Tfap2c_only"), "8cell"]
# W = 1326706, p-value = 0.0003429
# alternative hypothesis: true location shift is not equal to 0

gexpr <- data.frame(fread("~/work_space/4.ProjectET/analysis/fig1/build/0.whole_transcriptome/Normal_embryo_log2FPKM.csv"),check.names = F)
gexpr <- gexpr[,c("gene","Morula")]
t1_only <- data.frame(fread("~/work_space/4.ProjectET/analysis/fig1/2.TF_binding_gexpr/T1_only_8Cell.promoter",header = F),check.names = F)
t1_only$class <- "Tfap2c_only"
t2_only <- data.frame(fread("~/work_space/4.ProjectET/analysis/fig1/2.TF_binding_gexpr/T2_only_8Cell.promoter",header = F),check.names = F)
t2_only$class <- "Tead4_only"
tt_both <- data.frame(fread("~/work_space/4.ProjectET/analysis/fig1/2.TF_binding_gexpr/T1T2_8Cell.promoter",header = F),check.names = F)
tt_both$class <- "TT_both"
tt_details <- rbind(t1_only,t2_only,tt_both)



gexpr_tt <- merge(gexpr,tt_details,by.x="gene",by.y="V1")
wilcox.test(gexpr_tt[which(gexpr_tt$class == "Tead4_only"),"Morula"],gexpr_tt[which(gexpr_tt$class == "Tfap2c_only"),"Morula"])
# Wilcoxon rank sum test with continuity correction
# 
# data:  gexpr_tt[which(gexpr_tt$class == "Tead4_only"), "Morula"] and gexpr_tt[which(gexpr_tt$class == "Tfap2c_only"), "Morula"]
# W = 1332858, p-value = 0.0001456
# alternative hypothesis: true location shift is not equal to 0
p <- ggboxplot(gexpr_tt, "class", "Morula", color = "black", 
               width=0.5,alpha=0.8) 
ggpar(p,x.text.angle=45)
ggsave("../TT_co-binding_8c_gexpr_morula_barplot.pdf",width =3,height=3)



####################################################
########  Stage III: Dynamic Change: Tfap2c/Tead4
####################################################
library(tidyr)
library(ggplot2)
library(ggsci)
setwd("~/work_space/4.ProjectET/analysis/fig1/build/2.Tfap2c_peak/rpkm15/tmp/PeakAnno")
anno <- read.table("all_log2.anno",check.names = F)
anno <- anno[,c("8Cell_Tfap2c","Morula_Tfap2c","ICM_Tfap2c","EPI_E65_Tfap2c","ESC_Tfap2c","TE_Tfap2c","EXE_E65_Tfap2c","TSC_Tfap2c")]
colnames(anno) <- c("8-Cell","Morula","ICM","EPI","ESC","TE","EXE","TSC")
anno$class <- rownames(anno)
anno <- anno[c("Promoter","5UTR","Exon","Intron","3UTR","Intergenic","LINE","SINE","Low_complexity","Satellite"),]
anno2 <- gather(anno,timepoint,value,-class)
anno2$class <- factor(anno2$class,level=c("Promoter","5UTR","Exon","Intron","3UTR","Intergenic","LINE","SINE","Low_complexity","Satellite"))
anno2$timepoint <- factor(anno2$timepoint,levels=c("8-Cell","Morula","ICM","EPI","ESC","TE","EXE","TSC"))


ggplot(anno2,aes(x=timepoint,y=value)) +
  geom_bar(stat="identity",position=position_dodge(),aes(fill=class))  +
  scale_fill_npg() +
  theme_classic()
ggsave("~/work_space/4.ProjectET/analysis/fig1/Homer_Tfap2c_anno_bar.pdf",width = 7,height = 3)


setwd("~/work_space/4.ProjectET/analysis/fig1/build/2.Tead4_peak/rpkm15/tmp/PeakAnno")
anno <- read.table("all_log2.anno",check.names = F)
anno <- anno[,c("8Cell_Tead4","Morula_Tead4","TE_Tead4","EXE_E65_Tead4","TSC_Tead4")]
colnames(anno) <- c("8-Cell","Morula","TE","EXE","TSC")
anno$class <- rownames(anno)
anno <- anno[c("Promoter","5UTR","Exon","Intron","3UTR","Intergenic","LINE","SINE","Low_complexity","Satellite"),]
anno2 <- gather(anno,timepoint,value,-class)
anno2$class <- factor(anno2$class,level=c("Promoter","5UTR","Exon","Intron","3UTR","Intergenic","LINE","SINE","Low_complexity","Satellite"))
anno2$timepoint <- factor(anno2$timepoint,levels=c("8-Cell","Morula","TE","EXE","TSC"))


ggplot(anno2,aes(x=timepoint,y=value)) +
  geom_bar(stat="identity",position=position_dodge(),aes(fill=class))  +
  scale_fill_npg() +
  theme_classic()
ggsave("~/work_space/4.ProjectET/analysis/fig1/Homer_Tead4_anno_bar.pdf",width = 7,height = 3)


