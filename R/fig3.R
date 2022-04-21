library(data.table)
library(dplyr)
library(ggpubr)
library(pheatmap)
library(stringr)
library(tidyr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(ChIPseeker)
txdb<- TxDb.Mmusculus.UCSC.mm10.knownGene

setwd("~/work_space/4.ProjectET/analysis/fig3")

############
##### 1. PE TF: bubble plot
############
setwd("~/work_space/4.ProjectET/analysis/fig3/3.motif_bubble")
motif_gexpr <- data.frame(fread("motif_gexpr.txt"),check.names = F)
motif_gexpr$x1 = "a"

enh <- motif_gexpr[1:11,1:5]
enh$Gene <- factor(enh$Gene,levels = rev(enh$Gene))

p <- ggplot(enh,aes(x=x1, y=Gene)) +geom_point(aes(size=gexpr ,colour=pvalue)) 
p +  scale_color_gradient(low = "DarkSlateBlue", high = "Yellow3") + theme_bw()
ggsave("../Enhancer_cofactor_bubble.pdf",width = 3.5,height = 6)

pro <- motif_gexpr[12:21,1:5]
pro$Gene <- factor(pro$Gene,levels = rev(pro$Gene))

p <- ggplot(pro,aes(x=x1, y=Gene)) +geom_point(aes(size=gexpr ,colour=pvalue)) 
p +  scale_color_gradient(low = "DarkSlateBlue", high = "Yellow3") + theme_bw()
ggsave("../Promoter_cofactor_bubble.pdf",width = 3.5,height = 6)


############ selected_TF_expr
selected_motif_gexpr <- data.frame(fread("selected_motif_prolonged_gexpr.txt"),check.names = F)
colnames(selected_motif_gexpr) <- c("gene","D10.oocyte" ,"D14.oocyte" ,  "8w.oocyte" ,  "MII.oocyte"  , "PN5.zygote" ,  "Early.2.cell", "Late.2.cell" , "4.cell", "8.cell")
rownames(selected_motif_gexpr) <- selected_motif_gexpr$gene
selected_motif_gexpr <- selected_motif_gexpr[,-1]
p <- pheatmap(selected_motif_gexpr,cluster_rows = T,cluster_cols = F,color = colorRampPalette(colors = c("NavyBlue","black","Yellow1"))(100))
ggsave("../selected_TF_expr_heatmap.pdf",p$gtable,dpi=300,width=5,height=4)
# selected_motif_gexpr <- gather(selected_motif_gexpr,"stage","fpkm",-"gene")
# 
# p <- ggline(selected_motif_gexpr, "stage", "fpkm",
#             color = "gene", palette =  "paired")
# p<- p+ scale_color_npg()
# ggpar(p,x.text.angle = 90)
# ggsave("genomecov_aar_arr.pdf",dpi=300,width=3,height=4)


############
##### 2. PE-TT occupancy
############
setwd("~/work_space/4.ProjectET/analysis/fig3/4.PEloop_TToccu")
TToccu <- data.frame(fread("stat_on_TToccu_PE.txt"),check.names = F)
TToccu <- gather(TToccu,"type","num",-"class")
TToccu$type <- factor(TToccu$type,levels = c( "co_occu","only_Occu_by_Tfap2c","only_Occu_by_Tead4","none"))
ggbarplot(TToccu, "class", "num",
          fill = "type", color = "black", palette = c("#B22222","#FA8072","#FFC0CB","#FFFAFA"),
          label = TRUE, lab.col = "black", lab.pos = "in",width=0.5)
ggsave("../PEloop_TToccu_barplot.pdf",width = 3,height = 5)

############
##### 3. Science_TTkd_RNA
############
TTkd_raw <- data.frame(fread("~/work_space/Public_data/o.science_TTkd_8C_RNA/T1T2_KD_normat.csv"),check.names = F)
TTkd_raw <- TTkd_raw[which(rowSums(TTkd_raw[,2:l3]) > 8),]
TTkd_raw$ctrl <- rowMeans(TTkd_raw[,2:5])
TTkd_raw$T1kd <- rowMeans(TTkd_raw[,6:9])
TTkd_raw$T2kd <- rowMeans(TTkd_raw[,10:l3])
TTkd <- TTkd_raw[,c("name","ctrl","T1kd","T2kd")]
TTkd$t1x <- TTkd$T1kd - TTkd$ctrl
TTkd$t2y <- TTkd$T2kd - TTkd$ctrl

########九宫格版

all_fpkm2=TTkd
colnames(all_fpkm2)[1]="gene"
all_fpkm2$color = ifelse(abs(all_fpkm2$t1x) < 0.585, "grey", ifelse(abs(all_fpkm2$t2y) < 0.585, "grey", "black"))
# dux_target <- data.frame(fread("~/work_space/1.Mouse_Acetylation/buffet/resource_set/dux_affected_gene.txt",col.names = "V1"),check.names = F)$"V1"
# plt.temp <- all_fpkm2[which(all_fpkm2$gene %in% dux_target),]
# plt.temp <- spread(plt.temp, class, fpkm)

#plt.temp$fc <- plt.temp$dux - plt.temp$nt

#ggboxplot(plt.temp,"class","fpkm")

key_influ <- c("Dnmt3b","Nanog","Tfap2c","Tead4")
all_fpkm2$highlight <- ifelse(all_fpkm2$gene %in% key_influ, "yes", "no")
#all_fpkm2 <- spread(all_fpkm2,"class","fpkm")
all_fpkm2$lab = all_fpkm2$gene
all_fpkm2$lab <- ifelse(all_fpkm2$lab %in% key_influ,all_fpkm2$lab,NA)

p <-
  ggscatter(all_fpkm2, x = "t1x", y = "t2y",alpha = 1,color = "color", palette = c("Firebrick4","grey"),label = "lab",repel = F, legend="none" ) +
  geom_hline(linetype = "dashed",yintercept = 0.585,color = "grey") +
  geom_hline(linetype = "dashed",yintercept = -0.585,color = "grey") +
  geom_vline(linetype = "dashed",xintercept = 0.585,color = "grey") +
  geom_vline(linetype = "dashed",xintercept = -0.585,color = "grey") +
  xlim(-6,6) + ylim(-4,4)+coord_cartesian(expand = F)
p + border()
ggsave("../TTkd_RNA_scatter.pdf",height = 5,width = 5)


temp_gene <- c(all_fpkm2[which(all_fpkm2$t1x < -0.585 & all_fpkm2$t2y > 0.585),"gene"],
               all_fpkm2[which(all_fpkm2$t1x < -0.585 & all_fpkm2$t2y < 0.585 & all_fpkm2$t2y > -0.585),"gene"],
               all_fpkm2[which(all_fpkm2$t1x < -0.585 & all_fpkm2$t2y < -0.585),"gene"],
               all_fpkm2[which(all_fpkm2$t2y < -0.585 & all_fpkm2$t1x < 0.585 & all_fpkm2$t1x > -0.585),"gene"],
               all_fpkm2[which(all_fpkm2$t2y < -0.585 & all_fpkm2$t1x > 0.585 ),"gene"]
) 
temp_group <- c(rep("G1",length(all_fpkm2[which(all_fpkm2$t1x < -0.585 & all_fpkm2$t2y > 0.585),"gene"])),
                rep("G2",length(all_fpkm2[which(all_fpkm2$t1x < -0.585 & all_fpkm2$t2y < 0.585 & all_fpkm2$t2y > -0.585),"gene"])),
                rep("G3",length(all_fpkm2[which(all_fpkm2$t1x < -0.585 & all_fpkm2$t2y < -0.585),"gene"])),
                rep("G4",length(all_fpkm2[which(all_fpkm2$t2y < -0.585 & all_fpkm2$t1x < 0.585 & all_fpkm2$t1x > -0.585),"gene"])),
                rep("G5",length(all_fpkm2[which(all_fpkm2$t2y < -0.585 & all_fpkm2$t1x > 0.585 ),"gene"]))
) 
temp = data.frame(gene=temp_gene,group=temp_group)
fwrite(temp,"../TTkd_RNA_scatter_group.csv")


selected_gene <- data.frame(fread("~/work_space/4.ProjectET/analysis/fig3/TTkd_selected_Tenary_plot.txt"),check.names = F)
gene_sort <- selected_gene$name
selected_gene <- merge(selected_gene,TTkd,by="name")
df2 <- selected_gene[,c("name","group","ctrl","T1kd","T2kd")]
df2 <- gather(df2,"sample","fpkm",-"name",-"group")
df2$sample <- factor(df2$sample,levels = c("T1kd","ctrl","T2kd"))
df2$name <- factor(df2$name,levels = gene_sort)
df2 <- df2[order(df2$group,decreasing = T),]

ggdotchart(df2, x = "name", y = "fpkm",rotate=F,
           color = "sample", size = 3,
           add = "segment",
           add.params = list(color = "lightgray", size = 1.5),
           position = position_dodge(0.3),
           palette = "npg",
           ggtheme = theme_pubclean(),
           sorting="none"
)
ggsave("../TTkd_RNA_selected_dotchart.pdf",height = 5,width = 10)

# ############# Tenary plot
# 
# selected_gene <- data.frame(fread("~/work_space/4.ProjectET/analysis/fig3/TTkd_all_selected_Tenary_plot.txt"),check.names = F)
# selected_gene <- merge(selected_gene,TTkd,by="name")
# library("ggtern")
# p1<-ggtern(data=selected_gene,aes(x=ctrl,y=T1kd,z=T2kd))+
#   geom_point(aes(size=ctrl,color=group),alpha=0.8)


############
##### 4. Klf5 peak dynamics
############

# 40036 8Cell_Klf5_ds_rpkm15filter.bed
# 21268 8Cell_siTfap2c_Klf5_ds_rpkm15filter.bed
# 19591 Morula_Klf5_ds_rpkm15filter.bed

klf5_dynamics <- data.frame(sample=c("8Cell","8Cell_Tfap2cKD","Morula"),peaknum=c(40036,21268,19591))
ggbarplot(klf5_dynamics,"sample","peaknum",fill = "#12A64A",width = 0.5,alpha=0.8)
ggsave("../PEloop_TToccu_barplot.pdf",height = 4,width = 3)

############
##### 5. RAR/ALDH1 family expression
############
gexpr1 <-  data.frame(fread("~/work_space/Bio-Resource/embryos/RNA/Xiewei_2016_Nature_log2FPKM_RNA.csv"),check.names = F)
gexpr2 <- data.frame(fread("~/work_space/Public_data/q.E65_epiexe_xw_zy_RNA/Embryo_E65_RNA_yizhang_xw_collapse.csv"),check.names = F)
gexpr4 <- data.frame(fread("~/work_space/Bio-Resource/embryos/RNA/Normal_embryo_noem2c_log2FPKM.csv"),check.names = F)[,c("gene","Morula","TE","TSC" )]
gexpr <- merge(gexpr1,gexpr2, by.x="Gene", by.y="name")
gexpr <- merge(gexpr,gexpr4, by.x="Gene", by.y="gene")
colnames(gexpr) <- c("Gene" ,"D10.oocyte","D14.oocyte","8w.oocyte" ,"MII.oocyte" ,"PN5.zygote","Early.2.cell","Late.2.cell" ,
                     "4.cell", "8.cell" ,"ICM" ,"ESC" , "EPI"  ,"EXE" ,"morula" ,"TE" ,"TSC" )
gexpr <- gexpr[,c("Gene" ,"D10.oocyte","D14.oocyte","8w.oocyte" ,"MII.oocyte" ,"PN5.zygote","Early.2.cell","Late.2.cell" ,
                  "4.cell", "8.cell" ,"morula"  )]
tt_expr <- gexpr[gexpr$Gene %in% c("Rara","Rarb","Rarg"),]
tt_expr <- gather(tt_expr,"time","FPKM",-"Gene")
tt_expr$Gene <- factor(tt_expr$Gene,levels=c("Rara","Rarb","Rarg"))

p <- ggbarplot(tt_expr, "time", "FPKM",
               fill = "Gene", color = "black", palette = c("#E64B35","#4DBBD5","#00A087"),
               position = position_dodge(0.6),width=0.5)
ggpar(p,x.text.angle = 45)
ggsave("../Rar_expr.pdf",width = 6,height = 4)

tt_expr <- gexpr[gexpr$Gene %in% c("Aldh1a1","Aldh1a2","Aldh1a3"),]
tt_expr <- gather(tt_expr,"time","FPKM",-"Gene")
tt_expr$Gene <- factor(tt_expr$Gene,levels=c("Aldh1a1","Aldh1a2","Aldh1a3"))
p <- ggbarplot(tt_expr, "time", "FPKM",
               fill = "Gene", color = "black", palette = c("#E64B35","#4DBBD5","#00A087"),
               position = position_dodge(0.6),width=0.5)
ggpar(p,x.text.angle = 45)
ggsave("../Aldh1_expr.pdf",width = 6,height = 4)

############
##### 6. rarg ~ sine
############
setwd("~/work_space/4.ProjectET/analysis/fig3/6.Rarg_motif_repeats_demethy/2.rarg_sine")
rarg_stat <-  data.frame(fread("rarg_repeats_stat.txt"),check.names = F)
rarg_stat_tmp <- rarg_stat[rarg_stat$type == "class",]
ggpie(rarg_stat_tmp, "num", label = "anno",palette=pal_npg("nrc", alpha = 1)(10)[c(1,3,9)],
      fill = "anno", color = "white")
ggsave("../../rarg_repeatsanno_class_pie.pdf",width = 2)

rarg_stat_tmp <- rarg_stat[rarg_stat$type == "family",]
ggpie(rarg_stat_tmp, "num", label = "anno",palette=pal_npg("nrc", alpha = 1)(10)[c(2,5,7,6)],
      fill = "anno", color = "white")
ggsave("../../rarg_repeatsanno_family_pie.pdf",width = 3)

rarg_stat_tmp <- rarg_stat[rarg_stat$type == "name",]
ggpie(rarg_stat_tmp, "num", label = "anno",palette=pal_npg("nrc", alpha = 1)(10)[c(1,3,2,6,9,10)],
      fill = "anno", color = "white")
ggsave("../../rarg_repeatsanno_name_pie.pdf",width = 4)


############
##### 7. Rarg-SINE three levels - demethylation
############

####### Highlight quick demethylation in level 1
setwd("~/work_space/4.ProjectET/analysis/fig3/6.Rarg_motif_repeats_demethy/3.rarg_3levels_demethy") 
l1 <- data.frame(fread("l1"))
l1$zygote <- (l1$GSM1386019_oocyte + l1$GSM1386020_sperm)/2 
l1$tag <- "l1"
l2 <- data.frame(fread("l2"))
l2$zygote <- (l2$GSM1386019_oocyte + l2$GSM1386020_sperm)/2 
l2$tag <- "l2"
l3 <- data.frame(fread("l3"))
l3$zygote <- (l3$GSM1386019_oocyte + l3$GSM1386020_sperm)/2
l3$tag <- "l3"
l123 <- rbind(l1,l2) 
l123 <- rbind(l123, l3) 
colnames(l123) <-c("chr", "start" ,"end","oocyte","sperm","2cell", "4cell","zygote", "tag") 
#1123 <- 1123[, C"chr" , "start" , "end" , "2cell", "4cell","zygote", "tag")]

l123 <- gather(l123, "time","ratio",-"chr",-"start",-"end",-"tag")
l123$time <- factor(l123$time, levels = c("oocyte","sperm", "zygote","2cell", "4cell"))

p <- ggboxplot(l123, "time","ratio", facet.by = "tag")
ggpar(p,x.text.angle = 45)

ggsave("SINE_Rarg_3Levels.pdf",width = 7, height = 3)

####### Select quick demethylation in level 2
l2_quick <- l2
l2_quick <- l2_quick[which(l2_quick$GSM1386022_4cell - l2_quick$zygote < -0.8),]
fwrite(l2_quick[,c(1,2,3)] ,file = "L2_quick_demethylation.bed", quote = F,sep="\t",col.names = F)
peakAnno <- annotatePeak("L2_quick_demethylation.bed", tssRegion=c(-2000, 500),
                         TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnno_df <- data.frame(peakAnno)
peakAnno_df <- peakAnno_df[which(abs(peakAnno_df$distanceToTSS) < 20000),]
dql <- data.frame(fread("~/work_space/Bio-Resource/embryos/RNA/Deng_Science_allstage_collapse.txt"),check.names = F)
peakAnno_df <- merge(peakAnno_df,dql,by.x="SYMBOL",by.y="gene")
peakAnno_df <- merge(l2_quick,peakAnno_df,by.y=c("seqnames","end"),by.x=c("chr","end"))

rna_df1 <- peakAnno_df[,c(25:31)]
rna_df1 <- rna_df1[which(rowSums(rna_df1) > 5),]
#rna_df <- t(scale(t(rna_df)))
boxplot(rna_df1,las=2,ylim=c(0,10))

dql <- dql[-which(dql$gene %in% peakAnno_df$SYMBOL),]
rna_df2 <- dql[,c(2:8)]
rna_df2 <- rna_df2[which(rowSums(rna_df2) > 5),]
#rna_df <- t(scale(t(rna_df)))
boxplot(rna_df2,las=2,ylim=c(0,10))

wilcox.test(rna_df1$early2cell,rna_df2$early2cell) #2.5e-4
wilcox.test(rna_df1$late2cell,rna_df2$late2cell) #< 2.2e-16
wilcox.test(rna_df1$`4cell`,rna_df2$`4cell`)#< 2.2e-16

tmp0 <- data.frame(num=c(rna_df1$early2cell,rna_df2$early2cell),class=c(rep("quick",nrow(rna_df1)),rep("others",nrow(rna_df2))))
tmp1 <- data.frame(num=c(rna_df1$late2cell,rna_df2$late2cell),class=c(rep("quick",nrow(rna_df1)),rep("others",nrow(rna_df2))))
tmp2 <- data.frame(num=c(rna_df1$`4cell`,rna_df2$`4cell`),class=c(rep("quick",nrow(rna_df1)),rep("others",nrow(rna_df2))))
tmp <- rbind(tmp0,tmp1,tmp2)
tmp$stage <- c(rep("m2C",nrow(tmp0)),rep("l2C",nrow(tmp1)),rep("4C",nrow(tmp2)))
tmp$stage <- factor(tmp$stage,levels = c("m2C","l2C","4C"))
ggboxplot(tmp,"class","num",facet.by = "stage")
ggsave("../../l2_quick_demethy_gexpr.pdf",width = 5,height = 5)


############
##### 8. T1KD-epi
############
setwd("~/work_space/4.ProjectET/analysis/fig3/7.t1kd_epi")
epi_peaks_dy <- data.frame(fread("epi_peaks_change.stat"))
colnames(epi_peaks_dy) <- c("num","sample")
epi_peaks_dy$sample <- factor(epi_peaks_dy$sample,levels = c("8Cell_Ctrl_H3K4me3","8Cell_siTfap2c_H3K4me3","8Cell_Ctrl_H3K27ac_rep2","8Cell_siTfap2c_H3K27ac_rep1"))
p <- ggbarplot(epi_peaks_dy, "sample", "num", color = "black",fill="grey60",width=0.5,alpha=0.8)
ggpar(p,x.text.angle=45)
ggsave("../T1KD_epi_peaknum.pdf",width = 3, height =5)

epi_bw <- data.frame(fread("build/T1kd_epi.tab"),check.names = F)
epi_bw$ratio27 <- epi_bw$`8Cell_Ctrl_H3K27ac_rep2` - epi_bw$`8Cell_siTfap2c_H3K27ac_rep1`
epi_bw$ratio4 <- epi_bw$`8Cell_Ctrl_H3K4me3_ds` - epi_bw$`8Cell_siTfap2c_H3K4me3_ds`
enhancer_anno <- data.frame(fread("~/work_space/4.ProjectET/analysis/fig3/build/loop/8Cloop_raw.bedpe"))[,c(1,2,3,7)]
epi_bw <- merge(epi_bw,enhancer_anno,by.x=c("chr","start","end"),by.y=c("V1","V2","V3"))

epi_bw_thin <- gather(epi_bw[,4:7])
epi_bw_thin$key <- factor(epi_bw_thin$key,levels = c("8Cell_Ctrl_H3K4me3_ds","8Cell_siTfap2c_H3K4me3_ds","8Cell_Ctrl_H3K27ac_rep2","8Cell_siTfap2c_H3K27ac_rep1"))
p<- ggboxplot(epi_bw_thin,"key","value")
ggpar(p,legend = "none",xlab = FALSE,ylim=c(0,90),x.text.angle=45)
ggsave("../T1KD_downDEG_enhancer_epi_boxplot.pdf",width = 3,height = 4)

epi_bw2 <- data.frame(fread("build/T1kd_T2be_epi.tab"),check.names = F)
epi_bw2_thin <- gather(epi_bw[,4:7])
epi_bw2_thin$key <- factor(epi_bw2_thin$key,levels = c("8Cell_Ctrl_H3K4me3_ds","8Cell_siTfap2c_H3K4me3_ds","8Cell_Ctrl_H3K27ac_rep2","8Cell_siTfap2c_H3K27ac_rep1"))
p<- ggboxplot(epi_bw2_thin,"key","value")
ggpar(p,legend = "none",xlab = FALSE,ylim=c(0,90),x.text.angle=45)
ggsave("../T1KD_T2binding_enhancer_epi_boxplot.pdf",width = 3,height = 4)
