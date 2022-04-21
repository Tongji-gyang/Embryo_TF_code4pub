setwd("~/work_space/4.ProjectET/analysis/figs3")
library(data.table)
library(stringr)

############
##### 1. Output allele-biased peaks
############
  
allele_df <- data.frame(fread("~/work_space/4.ProjectET/analysis/figs3/project_et_t2_allele.counts"),check.names = F)
allele_df <- allele_df[which(allele_df$chr != "chrX" ),]
allele_df <- allele_df[which(allele_df$chr != "chrY" ),]
allele_df <- allele_df[which(allele_df$chr != "chrM" ),]
#2463178 Kb

for (i in 1:6) {
  name=gsub(".c57","",colnames(allele_df)[2*i+2])
  name=paste0(name,".ASraw.txt")
  df <- allele_df[,c(1,2,3,2*i+2,2*i+3)]
  df$total <- df[,4] + df[,5]
  df$lfc <- abs(log2(df[,4]+1) - log2(df[,5]+1))
  df <- df[which(df$total >= 100),]
  fwrite(df[,1:7],name,quote = F,sep = "\t",col.names = F)
}

library(stringr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(ChIPseeker)
txdb<- TxDb.Mmusculus.UCSC.mm10.knownGene

setwd("~/work_space/4.ProjectET/analysis/figs3/AS_stat")
peak_list <- list.files(pattern="peak.bed")

out_supp <- data.frame()
for (i in peak_list){
  name = gsub(".peak.bed","",i)
peakAnno <- data.frame(annotatePeak(i, tssRegion=c(-2500, 500),
                         TxDb=txdb, annoDb="org.Mm.eg.db"))
peakAnno <- peakAnno[,c("seqnames","start","end","SYMBOL", "annotation")]
peakAnno$origin = name
out_supp <- rbind(out_supp,peakAnno)
}
fwrite(out_supp,"../Allele_biased_peaks.txt",quote = F,sep = "\t",col.names = F)


############
##### 2. ABP versus SNP？
############
setwd("~/work_space/4.ProjectET/analysis/figs3/tf_snp/1.ABP_loci/check_snp_pm")

motif_snp_total <- data.frame(fread("tf_motif_snp_stat_tib.txt"),check.names = F)
motif_snp_total$xaxis <- paste0(motif_snp_total$tf,"_",motif_snp_total$allelic_type)
motif_snp_total$xaxis <- factor(motif_snp_total$xaxis,levels = c("Tfap2c_maternal","Tfap2c_paternal","Tead4_maternal","Tead4_paternal"))

p <- ggbarplot(motif_snp_total, "xaxis", "num",
               fill = "snp_type", color = "snp_type", palette = c("IndianRed1","SteelBlue"),width = 0.5)
p <- ggpar(p,legend = "none",x.text.angle = 45)
ggsave("~/work_space/4.ProjectET/analysis/figs3/snp_on_motif_summarybar.pdf",width = 3,height = 5.5)

motif_snp_logo <- data.frame(fread("TF_SNP_position_stat.txt"),check.names = F)
temp <- motif_snp_logo[motif_snp_logo$tf == "Tfap2c",]
# p <- ggline(temp, "position", "snp_percent",
#             color = "allelic_type", palette =  c("IndianRed1","SteelBlue"))
# ggpar(p,width=6,height=4)
# ggsave("~/work_space/4.ProjectET/analysis/figs3/snp_on_motif_base_line_tfap2c.pdf",width = 5,height = 3)

temp$allelic_type <- factor(temp$allelic_type,levels = c("maternal","paternal"))
p <- ggbarplot(temp, "position", "snp_percent",
               fill = "allelic_type", color = "white", palette = c("#282A73","#96CB6B"),alpha=0.8,width = 0.5,position = position_dodge(0.5))
p <- ggpar(p,legend = "none")
ggsave("~/work_space/4.ProjectET/analysis/figs3/snp_on_motif_base_bar_tfap2c.pdf",width = 5,height = 3)

motif_snp_logo <- data.frame(fread("TF_SNP_position_stat.txt"),check.names = F)
temp <- motif_snp_logo[motif_snp_logo$tf == "Tead4",]
# p <- ggline(temp, "position", "snp_percent",
#             color = "allelic_type", palette =  c("IndianRed1","SteelBlue"))
# ggpar(p,width=6,height=4)
# ggsave("~/work_space/4.ProjectET/analysis/figs3/snp_on_motif_base_line_tead4.pdf",width = 5,height = 3)
temp$allelic_type <- factor(temp$allelic_type,levels = c("maternal","paternal"))
p <- ggbarplot(temp, "position", "snp_percent",
               fill = "allelic_type", color = "white", palette = c("#282A73","#96CB6B"),alpha=0.8,width = 0.5,position = position_dodge(0.5))
p <- ggpar(p,legend = "none")
ggsave("~/work_space/4.ProjectET/analysis/figs3/snp_on_motif_base_bar_tead4.pdf",width = 5,height = 3)

############
#### 3. Epi at sight
############
library(ggpubr)
setwd("~/work_space/4.ProjectET/analysis/figs3/epi_at_sight")
fl=list.files(pattern = "epi_atsight.tab")
for(i in fl){
  print(i)
  df <- data.frame(fread(i))
  df <- df[,4:ncol(df)]
  df <- gather(df,"epi","value")
  p <- ggboxplot(df,"epi","value")
  p <- ggpar(p,legend = "none",xlab = FALSE,x.text.angle=45)
  ggsave(gsub(".tab",".pdf",i),plot=p,width = 4,height = 3,useDingbats=F)
  #boxplot(df[,4:ncol(df)],las=2,ylim = c(-3,3))
  
}


asg <- data.frame(fread("/home1/gyang/work_space/Bio-Resource/embryos/RNA/all_PWK_snp.csv"),check.names = F)[,c(1,10:21)]
fwrite(data.frame(asg[asg$`8cell-c57`-asg$`8cell-pwk` > 1.58,"name"]),
       "/home1/gyang/work_space/4.ProjectET/analysis/figs3/AS_stat/ASG_RNA_snp/8Cell_mat.csv")
fwrite(data.frame(asg[asg$`8cell-c57`-asg$`8cell-pwk` < -1.58,"name"]),
       "/home1/gyang/work_space/4.ProjectET/analysis/figs3/AS_stat/ASG_RNA_snp/8Cell_pat.csv")
fwrite(data.frame(asg[asg$`BP_Morula.c57`-asg$`BP_Morula.pwk` > 1.58,"name"]),
       "/home1/gyang/work_space/4.ProjectET/analysis/figs3/AS_stat/ASG_RNA_snp/Morula_mat.csv")
fwrite(data.frame(asg[asg$`BP_Morula.c57`-asg$`BP_Morula.pwk` < -1.58,"name"]),
       "/home1/gyang/work_space/4.ProjectET/analysis/figs3/AS_stat/ASG_RNA_snp/Morula_pat.csv")
fwrite(data.frame(asg[asg$`ICM-c57`-asg$`ICM-pwk` > 1.58,"name"]),
       "/home1/gyang/work_space/4.ProjectET/analysis/figs3/AS_stat/ASG_RNA_snp/ICM_mat.csv")
fwrite(data.frame(asg[asg$`ICM-c57`-asg$`ICM-pwk` < -1.58,"name"]),
       "/home1/gyang/work_space/4.ProjectET/analysis/figs3/AS_stat/ASG_RNA_snp/ICM_pat.csv")
fwrite(data.frame(asg[asg$`E65EPI-c57dba`-asg$`E65EPI-pwk` > 1.58,"name"]),
       "/home1/gyang/work_space/4.ProjectET/analysis/figs3/AS_stat/ASG_RNA_snp/EPI_E65_mat.csv")
fwrite(data.frame(asg[asg$`E65EPI-c57dba`-asg$`E65EPI-pwk` < -1.58,"name"]),
       "/home1/gyang/work_space/4.ProjectET/analysis/figs3/AS_stat/ASG_RNA_snp/EPI_E65_pat.csv")
fwrite(data.frame(asg[asg$`E65EXE-c57dba`-asg$`E65EXE-pwk` > 1.58,"name"]),
       "/home1/gyang/work_space/4.ProjectET/analysis/figs3/AS_stat/ASG_RNA_snp/EXE_E65_mat.csv")
fwrite(data.frame(asg[asg$`E65EXE-c57dba`-asg$`E65EXE-pwk` < -1.58,"name"]),
       "/home1/gyang/work_space/4.ProjectET/analysis/figs3/AS_stat/ASG_RNA_snp/EXE_E65_pat.csv")
fwrite(data.frame(asg[asg$`TE_c57dba`-asg$`TE_pwk` > 1.58,"name"]),
       "/home1/gyang/work_space/4.ProjectET/analysis/figs3/AS_stat/ASG_RNA_snp/TE_mat.csv")
fwrite(data.frame(asg[asg$`TE_c57dba`-asg$`TE_pwk`< -1.58,"name"]),
       "/home1/gyang/work_space/4.ProjectET/analysis/figs3/AS_stat/ASG_RNA_snp/TE_pat.csv")

# tf <- c('Tfap2c',"Tead4")
# all <- c('mat',"pat")
# for(i in tf) {
#   for(j in all) {
#     print(i)
#     print(j)
#     fl=list.files(pattern = paste0(i,"_",j))
#     df <- data.frame(fread(fl[1]))
#     for(k in fl[2:length(fl)]) {
#       df <- merge(df,data.frame(fread(k)))
#     }
#     
#   }
# }



######
###  Cross validation failed
######
# dql_16c <- data.frame(fread("~/work_space/4.ProjectET/analysis/figs3/Dql_16C_allele_merged.csv"))
# m_mat <- data.frame(fread("~/work_space/4.ProjectET/analysis/figs3/AS_stat/ASG_RNA_snp/Morula_mat.csv"))[,1]
# p_mat <- data.frame(fread("~/work_space/4.ProjectET/analysis/figs3/AS_stat/ASG_RNA_snp/Morula_pat.csv"))[,1]
# boxplot(log2(dql_16c[which(dql_16c$Gene %in% m_mat),c("CAST_hits","C57_hits")]+1))

all_abp_gene <- data.frame(fread("~/work_space/4.ProjectET/analysis/figs3/AS_stat/ASP_10k/all.rna.bed"))
all_modi <- data.frame(fread("~/work_space/4.ProjectET/analysis/figs3/epi_at_sight/all.modi.tab"))
all_modi <- merge(all_modi,all_abp_gene,by=c("V1","V2","V3"))
colnames(all_modi) <- c("chr","start","end","m_K27me3","p_k27me3","m_K9me3","p_k9me3","class","related_gene")  
setwd("~/work_space/4.ProjectET/analysis/figs3/epi_at_sight")
all_modi <- all_modi[order(all_modi$class),]
fwrite(all_modi,"ABP_allele_epi_genebias.csv")
tmp <- c("1700012B09Rik", "1700102H20Rik", "4933406M09Rik", "AA388235", "Adamts14", "Adamts9", "Adora1", "Ankrd33b", "Apbb1", "Arcn1", "Ascl2", "BC048562", "Ccdc36", "Ccnd3", "Chpf", "Chrna1", "Copg2", "Cyb5r1", "Cyth3", "D930015E06Rik", "Dcaf4", "Dhcr7", "Dnaaf1", "Dok5", "E030030I06Rik", "Eif4enif1", "Eif4enif1", "Fam89a", "Fancg", "Fut4", "Gab1", "Gm15910", "Gna15", "H19", "H19", "H2-K2", "H2-K2", "Hrc", "Hsdl1", "Htatip2", "Ift46", "Ildr1", "Itgb7", "Kl", "Mbtps1", "Meg3", "Mest", "Mest", "Mir335", "Mir335", "Mir675", "Mir675", "Msrb3", "Nadsyn1", "Nup210l", "Obsl1", "Peg3", "Peg3", "Peg3", "Phlda2", "Pigo", "Pisd-ps1", "Pisd-ps1", "Pisd-ps1", "Pisd-ps1", "Ppp2r2c", "Pxt1", "Raet1a", "Rarg", "Rasa3", "Rims2", "Rpl14-ps1", "Rsph1", "Rsph1", "Rsph1", "Rsph1", "S1pr4", "Sema5a", "Sfi1", "Sfi1", "Sfi1", "Sfi1", "Sfi1", "Sfi1", "Sfi1", "Sfi1", "Sfi1", "Sfi1", "Sfi1", "Sfi1", "Sfi1", "Sfi1", "Sfmbt2", "Sfmbt2", "Sh3bp5", "Sh3bp5", "Sh3bp5", "Slc22a18", "Smoc1", "Smpd1", "Snx30", "Stoml2", "Tfb1m", "Tmem198", "Trim62", "Trpm4", "Ubash3a", "Usp29", "Usp29", "Usp29", "Vegfa", "Vps52")
tmp <- asg[which(asg$name %in% tmp),]


####### Vegfa Ppp2r2c expr barplot

library(ggpubr)
library(stringr)
library(tidyr)
colnames(asg) <- gsub("Morula.","Morula-",colnames(asg))
colnames(asg) <- gsub("dba","",colnames(asg))
colnames(asg) <- gsub("BP_","",colnames(asg))
colnames(asg) <- gsub("E65","",colnames(asg))
colnames(asg) <- gsub("c57","Maternal",colnames(asg))
colnames(asg) <- gsub("pwk","Paternal",colnames(asg))
colnames(asg) <- gsub("_","-",colnames(asg))

tmp <- asg[asg$name %in% c("Gab1","Meg3","Ankrd33b"),]
tmp <- gather(tmp,"class","value",-"name")
tmp$stage <- str_split_fixed(tmp$class,"-",2)[,1]
tmp$allele <- str_split_fixed(tmp$class,"-",2)[,2]
tmp$value <- 2^(tmp$value) -1
ggbarplot(tmp, "stage", "value",
          fill = "allele", color = "allele", palette = c("#2E2F73","#E14B36"),
          label = F, lab.col = "white", lab.pos = "in",facet.by = "name")
ggsave("Gab1_Meg3_Ankrd33b_allele_gexpr.pdf",width = 6.5,height = 4)

tmp <- asg[asg$name %in% c("Ppp2r2c","Pxt1","E030030I06Rik"),]
tmp <- gather(tmp,"class","value",-"name")
tmp$stage <- str_split_fixed(tmp$class,"-",2)[,1]
tmp$allele <- str_split_fixed(tmp$class,"-",2)[,2]
tmp$value <- 2^(tmp$value) -1
ggbarplot(tmp, "stage", "value",
          fill = "allele", color = "allele", palette = c("#2E2F73","#E14B36"),
          label = F, lab.col = "white", lab.pos = "in",facet.by = "name")
ggsave("Ppp2r2c_Pxt1_E030030I06Rik_allele_gexpr.pdf",width = 6.5,height = 4)

tmp <- asg[asg$name == "Ppp2r2c",]
tmp <- gather(tmp,"class","value",-"name")
tmp$stage <- str_split_fixed(tmp$class,"-",2)[,1]
tmp$allele <- str_split_fixed(tmp$class,"-",2)[,2]
tmp$value <- 2^(tmp$value) -1
ggbarplot(tmp, "stage", "value",
          fill = "allele", color = "allele", palette = c("#2E2F73","#E14B36"),
          label = F, lab.col = "white", lab.pos = "in")
ggsave("Ppp2r2c_allele_gexpr.pdf",width = 2.5,height = 4)

tmp <- asg[asg$name == "Vegfa",]
tmp <- gather(tmp,"class","value",-"name")
tmp$stage <- str_split_fixed(tmp$class,"-",2)[,1]
tmp$allele <- str_split_fixed(tmp$class,"-",2)[,2]
tmp$value <- 2^(tmp$value) -1
ggbarplot(tmp, "stage", "value",
          fill = "allele", color = "allele", palette = c("#2E2F73","#E14B36"),
          label = F, lab.col = "white", lab.pos = "in")
ggsave("Vegfa_allele_gexpr.pdf",width = 2.5,height = 4)

tmp <- asg[asg$name == "Gab1",]
tmp <- gather(tmp,"class","value",-"name")
tmp$stage <- str_split_fixed(tmp$class,"-",2)[,1]
tmp$allele <- str_split_fixed(tmp$class,"-",2)[,2]
tmp$value <- 2^(tmp$value) -1
ggbarplot(tmp, "stage", "value",
          fill = "allele", color = "allele", palette = c("#2E2F73","#E14B36"),
          label = F, lab.col = "white", lab.pos = "in")
ggsave("Gab1_allele_gexpr.pdf",width = 2.5,height = 4)

tmp <- asg[asg$name == "Peg3",]
tmp <- gather(tmp,"class","value",-"name")
tmp$stage <- str_split_fixed(tmp$class,"-",2)[,1]
tmp$allele <- str_split_fixed(tmp$class,"-",2)[,2]
tmp$value <- 2^(tmp$value) -1
ggbarplot(tmp, "stage", "value",
          fill = "allele", color = "allele", palette = c("#2E2F73","#E14B36"),
          label = F, lab.col = "white", lab.pos = "in")
ggsave("Peg3_allele_gexpr.pdf",width = 2.5,height = 4)

tmp <- asg[asg$name == "Pxt1",]
tmp <- gather(tmp,"class","value",-"name")
tmp$stage <- str_split_fixed(tmp$class,"-",2)[,1]
tmp$allele <- str_split_fixed(tmp$class,"-",2)[,2]
tmp$value <- 2^(tmp$value) -1
ggbarplot(tmp, "stage", "value",
          fill = "allele", color = "allele", palette = c("#2E2F73","#E14B36"),
          label = F, lab.col = "white", lab.pos = "in")
ggsave("Pxt1_allele_gexpr.pdf",width = 2.5,height = 4)

tmp <- asg[asg$name == "Ankrd33b",]
tmp <- gather(tmp,"class","value",-"name")
tmp$stage <- str_split_fixed(tmp$class,"-",2)[,1]
tmp$allele <- str_split_fixed(tmp$class,"-",2)[,2]
tmp$value <- 2^(tmp$value) -1
ggbarplot(tmp, "stage", "value",
          fill = "allele", color = "allele", palette = c("#2E2F73","#E14B36"),
          label = F, lab.col = "white", lab.pos = "in")
ggsave("Ankrd33b_allele_gexpr.pdf",width = 2.5,height = 4)

tmp <- asg[asg$name == "E030030I06Rik",]
tmp <- gather(tmp,"class","value",-"name")
tmp$stage <- str_split_fixed(tmp$class,"-",2)[,1]
tmp$allele <- str_split_fixed(tmp$class,"-",2)[,2]
tmp$value <- 2^(tmp$value) -1
ggbarplot(tmp, "stage", "value",
          fill = "allele", color = "allele", palette = c("#2E2F73","#E14B36"),
          label = F, lab.col = "white", lab.pos = "in")
ggsave("E030030I06Rik_allele_gexpr.pdf",width = 2.5,height = 4)

tmp <- asg[asg$name == "Rsph1",]
tmp <- gather(tmp,"class","value",-"name")
tmp$stage <- str_split_fixed(tmp$class,"-",2)[,1]
tmp$allele <- str_split_fixed(tmp$class,"-",2)[,2]
tmp$value <- 2^(tmp$value) -1
ggbarplot(tmp, "stage", "value",
          fill = "allele", color = "allele", palette = c("#2E2F73","#E14B36"),
          label = F, lab.col = "white", lab.pos = "in")
ggsave("Rsph1_allele_gexpr.pdf",width = 2.5,height = 4)

E65_m_sele <- data.frame(fread("/home1/gyang/work_space/4.ProjectET/analysis/figs3/E65_m",header=FALSE))$V1
tmp <- asg[asg$name %in% E65_m_sele,c("name","EXE-Maternal","EXE-Paternal")]
tmp <- tmp[tmp$`EXE-Maternal` + tmp$`EXE-Paternal` > 5,]
tmp[,2:3] <- 2^(tmp[,2:3]) + 1
tmp$rel <- tmp$`EXE-Maternal`/tmp$`EXE-Paternal`
tmp$fake <- 1 
ind <- tmp$name[order(tmp$rel,decreasing = T)]
tmp <- tmp[,c("name","rel","fake")]
tmp <- gather(tmp,"class","value",-"name")
tmp$name <- factor(tmp$name, levels = ind)
tmp$class <- factor(tmp$class, levels = c('fake','rel'))

library(ggplot2)
library(ggbreak)
ggplot(tmp,aes(name,value,fill=class)) + 
  geom_bar(stat="identity",position=position_dodge(width = 0.7),width = 0.5,) +
  scale_y_break(c(6, 20),scales = 0.6) + scale_fill_manual(values=c("#2E2F73","#E14B36")) + theme_bw()  +
  theme(panel.grid.major = element_blank(), #主网格线
          panel.grid.minor = element_blank(), #次网格线
          #panel.border = element_blank(), #边框
          #axis.title = element_blank(),  #轴标题
          #axis.text = element_blank(), # 文本
          #axis.ticks = element_blank()
        ) +  geom_hline(linetype = "dashed",yintercept = 2,color = "grey") +
  theme(axis.text.x = element_text(angle = 45)) 
ggsave("mat_allele_gexpr_abrplot.pdf",width = 8,height = 4)


E65_p_sele <- data.frame(fread("/home1/gyang/work_space/4.ProjectET/analysis/figs3/E65_p",header=FALSE))$V1
tmp <- asg[asg$name %in% E65_p_sele,c("name","EXE-Maternal","EXE-Paternal")]
tmp <- tmp[tmp$`EXE-Maternal` + tmp$`EXE-Paternal` > 5,]
tmp[,2:3] <- 2^(tmp[,2:3]) + 1
tmp$rel <- tmp$`EXE-Paternal`/tmp$`EXE-Maternal`
tmp$fake <- 1 
ind <- tmp$name[order(tmp$rel,decreasing = T)]
tmp <- tmp[,c("name","rel","fake")]
tmp <- gather(tmp,"class","value",-"name")
tmp$name <- factor(tmp$name, levels = ind)
tmp$class <- factor(tmp$class, levels = c('fake','rel'))

ggplot(tmp,aes(name,value,fill=class)) + 
  geom_bar(stat="identity",position=position_dodge(width = 0.7),width = 0.5,) +
  scale_y_break(c(5, 90),scales = 0.3) + scale_fill_manual(values=c("#2E2F73","#E14B36")) + theme_bw()  +
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        #panel.border = element_blank(), #边框
        #axis.title = element_blank(),  #轴标题
        #axis.text = element_blank(), # 文本
        #axis.ticks = element_blank()
        ) +  geom_hline(linetype = "dashed",yintercept = 2,color = "grey") +
  theme(axis.text.x = element_text(angle = 45)) 
ggsave("pat_allele_gexpr_barplot.pdf",width = 8,height = 4)

abp_2groups <- data.frame(fread("~/work_space/4.ProjectET/analysis/figs3/epi_at_sight/ABP_Morula_Dnase.tab",header=T))
colnames(abp_2groups) <- c("V1","V2","V3","AG","PG")
abp_2groups <- merge(abp_2groups,abp_2groups_anno, by=c("V1","V2","V3"))
abp_2groups_anno <-  data.frame(fread("~/work_space/4.ProjectET/analysis/figs3/epi_at_sight/ABP.anno",header=F))
abp_2groups <- abp_2groups[,4:7]
abp_2groups$abs <- abs(abp_2groups$AG - abp_2groups$PG)
abp_2groups <- abp_2groups[,3:5]
colnames(abp_2groups) <- c("class","prop","value")
#abp_2groups <- gather(abp_2groups,"class","value",-"Prop")
wilcox.test(abp_2groups[abp_2groups$class=="tranditional","value"],abp_2groups[abp_2groups$class=="denovo","value"])
p <- ggboxplot(abp_2groups, x = "class", y = "value",
          color = "black",
          width = 0.5, outlier.shape = NA)
ggsave("~/work_space/4.ProjectET/analysis/figs3/epi_at_sight/ABP_2group_morula_dnase.pdf")

#################
###  SNP analysis
#################

tf_snp <- data.frame(fread("~/work_space/4.ProjectET/analysis/figs3/tf_pwk_motif_snp_info_final_subset.txt"))
tf_snp <- tf_snp[which(tf_snp$exp_sum > 2),]
tmp <- tf_snp[which(!(tf_snp$annotated_imprints == T)),c(13,34:37)]
#tmp$var <- apply(tmp[,2:5],1,var)
#tmp$min <- apply(tmp[,2:5],1,min)
#tmp[,2:5] <- (tmp[,2:5] - tmp$min)/tmp$var
tmp[,2:5] <- t(scale(t(tmp[,2:5])))

 tmp$pp <- tmp$bp_e95_pla_pwk + tmp$pb_e95_pla_pwk
 tmp$bb <- tmp$bp_e95_pla_b6 + tmp$pb_e95_pla_b6
 tmp$class <- "snp"

tf_snp_2 <- data.frame(fread("~/work_space/1.Mouse_Acetylation/buffet/resource_set/mouse_imprinting_genes_v2.txt",header = F))$V1
tf_snp_2 <- gene_expression[which(gene_expression$name %in% tf_snp_2),]
tf_snp_2 <- tf_snp_2[rowSums(tf_snp_2[,1:4]) > 2,]
tmp2 <- tf_snp_2[,c(5,1,2,3,4)]

# tmp2$var <- apply(tmp2[,2:5],1,var)
# tmp2$min <- apply(tmp2[,2:5],1,min)
# tmp2[,2:5] <- (tmp2[,2:5] - tmp2$min)/tmp2$var

tmp2[,2:5] <- t(scale(t(tmp2[,2:5])))
tmp2$bb <- tmp2$BP_E95pla.B6 + tmp2$PB_E95pla.B6
tmp2$pp <- tmp2$BP_E95pla.PWK + tmp2$PB_E95pla.PWK
tmp2$class <- "imprinting"

tmp_all <- rbind(tmp[,6:8],tmp2[,6:8])


#tmp_all$name2 <- str_c(tmp_all$class,tmp_all$name)

tf_nosnp <- data.frame(fread("~/work_space/4.ProjectET/analysis/figs3/epi_at_sight/de_novo.gene",header = F))$V1
tf_nosnp <- gene_expression[gene_expression$name %in% tf_nosnp,]
tf_nosnp$exp_sum <- rowSums(tf_nosnp[,1:4])
tf_nosnp <- tf_nosnp[which(tf_nosnp$exp_sum > 2),]
#tmp$var <- apply(tmp[,2:5],1,var)
#tmp$min <- apply(tmp[,2:5],1,min)
#tmp[,2:5] <- (tmp[,2:5] - tmp$min)/tmp$var
tmp3 <- tf_nosnp[,c(5,1:4)]
tmp3[,2:5] <- t(scale(t(tmp3[,2:5])))

tmp3$pp <- tmp3$BP_E95pla.PWK + tmp3$PB_E95pla.PWK
tmp3$bb <- tmp3$BP_E95pla.B6+ tmp3$PB_E95pla.B6
tmp3$class <- "denovo"
unknown_gene <- tmp3[abs(tmp3$pp) < 1,"name"]

tmp3 <- tmp3[,6:8]
tmp_all <- rbind(tmp_all ,tmp3)

tmp_all <- gather(tmp_all,"name","value",-"class")
#ggscatter(tmp_all,"bb","pp",color = "name")

tmp_all <- tmp_all[!(tmp_all$class == "denovo"), ]
p <- ggdensity(tmp_all,"value",color = "class", rug=T)
ggpar(p,xlim = c(-1.75,1.75))
ggsave("/home1/gyang/work_space/4.ProjectET/analysis/figs3/Imprinting_SNP_backcross.pdf",width = 3.5,height = 3)

library(clusterProfiler)
library(ChIPseeker)
library(data.table)
library(ggpubr)
library(org.Mm.eg.db)
library(stringr)
library(tidyr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene


# peakAnno <- annotatePeak("~/work_space/4.ProjectET/analysis/figs3/tf_snp/3.ABP10k_ABW_motif_genomeanno/ABP10k_ABW_motif.bed", tssRegion=c(-2000, 500),
#                          TxDb=txdb, annoDb="org.Mm.eg.db")
# plotPeakProf2(peak = readPeakFile("~/work_space/4.ProjectET/analysis/figs3/tf_snp/3.ABP10k_ABW_motif_genomeanno/ABP10k_ABW_motif.bed"), upstream = rel(0.2), downstream = rel(0.2),
#               conf = 0.95, by = "gene", type = "body", nbin = 800,
#               TxDb = txdb,ignore_strand = F)
# plotAnnoBar(peakAnno)

df <- data.frame(
  group = c("Intron", "Others"),
  value = c(73, 47))
ggpie(df, "value", label = "group",fill = "group")
ggsave("/home1/gyang/work_space/4.ProjectET/analysis/figs3/AAM_genome_anno.pdf",width = 2.5,height = 3)

df <- data.frame(
  group = c("Intron", "Others"),
  value = c(752, 1885-752))
ggpie(df, "value", label = "group",fill = "group")
ggsave("/home1/gyang/work_space/4.ProjectET/analysis/figs3/EQTL_HG38_TFAP2C_anno.pdf",width = 2.5,height = 3)

df <- data.frame(
  group = c("Imprinting", "Others"),
  value = c(35, 526))
ggpie(df, "value", label = "group",fill = "group")
ggsave("/home1/gyang/work_space/4.ProjectET/analysis/figs3/Imprinting_ABPgene_composition.pdf",width = 2.5,height = 3)

df <- data.frame(
  group = c("Imprinting", "Others"),
  value = c(35, 526))
ggpie(df, "value", label = "group",fill = "group")
ggsave("/home1/gyang/work_space/4.ProjectET/analysis/figs3/Imprinting_ABPgene_composition.pdf",width = 2.5,height = 3)

df <- data.frame(
  species = c("Rat", "Pig","Zebrafish"),
  expr = c(20, 11,3),
  ctrl = c(34, 15,0))
df[,2:3] = df[,2:3]/37
df <- gather(df,group,value,-species)
p <- ggbarplot(df, "species", "value",
          fill = "group", color = "black", palette = c("blue","grey60"),width = 0.5,
          position = position_dodge(0.6))
ggpar(p,ylim = c(0,1))
ggsave("/home1/gyang/work_space/4.ProjectET/analysis/figs3/Evolution_SNP_enhancer_liftover_box.pdf",width = 2.5,height = 3)
######## Backcrossing bias definition

gevo <- data.frame(fread("/home1/gyang/work_space/Bio-Resource/evolution/MouseEvo_20groups.csv"),check.names = F)[,1:2]
colnames(gevo) <- c('gene',"group")
ids <- bitr(gevo[,"gene"], fromType="ENSEMBL", toType=c("SYMBOL"), OrgDb="org.Mm.eg.db")
ids <- merge(ids,gevo,by.x="ENSEMBL",by.y="gene")
fwrite(back_cross,"/home1/gyang/work_space/Bio-Resource/evolution/MouseEvo_20groups_transid.csv")


back_cross <- gene_expression
back_cross$t1 <- back_cross$BP_E95pla.PWK - back_cross$BP_E95pla.B6
back_cross$t2 <- back_cross$PB_E95pla.PWK - back_cross$PB_E95pla.B6
back_cross$color = ifelse(abs(back_cross$t1) > 1.5 & abs(back_cross$t2) > 1.5, ifelse(back_cross$t1 * back_cross$t2 > 0, "red", "blue"),'grey')
back_cross$color <- factor(back_cross$color, levels = c("red","blue","grey"))

# a <- ids[ids$SYMBOL %in% back_cross[back_cross$color == "red","name"],c("SYMBOL","group")]
# a$class="strain"
# b <- ids[ids$SYMBOL %in% back_cross[back_cross$color == "blue","name"],c("SYMBOL","group")]
# b$class="none"
# ab <- rbind(a,b)
# c <- data.frame(num=0:11, time=c(454.6,361.2,324.5,220.2,176.1,104.7,97.4,91,44.2,29.6,18.8,15.1))
# ab <- merge(ab,c,by.x="branch",by.y="num")
# ggecdf(ab, x = "time",color = "class")

 ggo <- groupGO(gene     = a$SYMBOL,
                OrgDb    = org.Mm.eg.db,
                ont      = "BP",
                keyType  = 'SYMBOL',
                level    = 3)
 barplot(ggo, showCategory=10) 
 
ego2 <- enrichGO(gene         = a$SYMBOL,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
dotplot(ego2, showCategory=10) + ggtitle("GO for Down-regulated genes") 

ego2 <- enrichGO(gene         = b$SYMBOL,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
dotplot(ego2, showCategory=10) + ggtitle("GO for Down-regulated genes") 

fwrite(back_cross,"/home1/gyang/work_space/4.ProjectET/analysis/figs3/Backcross_RNA_Scattergroup_anno.csv")
write.csv(back_cross[back_cross$color == "red","name"],"/home1/gyang/work_space/4.ProjectET/analysis/figs3/Backcross_RNA_Scattergroup_anno_redname.csv",quote = F,col.names = F,row.names = F)
write.csv(back_cross[back_cross$color == "grey","name"],"/home1/gyang/work_space/4.ProjectET/analysis/figs3/Backcross_RNA_Scattergroup_anno_greyname.csv",quote = F,col.names = F,row.names = F)
write.csv(back_cross[back_cross$color == "blue","name"],"/home1/gyang/work_space/4.ProjectET/analysis/figs3/Backcross_RNA_Scattergroup_anno_bluename.csv",quote = F,col.names = F,row.names = F)

t1 <- data.frame(fread("~/work_space/4.ProjectET/analysis/figs3/red_genebody.bed.enhancer.bed.pwk.snpcount",header = F))$V4
t2 <- data.frame(fread("~/work_space/4.ProjectET/analysis/figs3/blue_genebody.bed.enhancer.bed.pwk.snpcount",header = F))$V4
t3 <- data.frame(fread("~/work_space/4.ProjectET/analysis/figs3/grey_genebody.bed.enhancer.bed.pwk.snpcount",header = F))$V4

wilcox.test(t1,t2) #5.387e-05
wilcox.test(t1,t3) # 0.02725

plt.tmp <- data.frame(value=c(t1,t2,t3),tag=c(rep("red",length(t1)),rep("blue",length(t2)),rep("grey",length(t3))))
p <- ggboxplot(plt.tmp,"tag","value",width = 0.4)
ggpar(p,ylim=c(0,25))
ggsave("/home1/gyang/work_space/4.ProjectET/analysis/figs3/backcross_3groups_enhancer_PWksnp_boxplot.pdf",width = 2.5,height = 3.5)

t1 <- data.frame(fread("~/work_space/4.ProjectET/analysis/figs3/red_genebody.bed.enhancer.bed.CAST.snpcount",header = F))$V4
t2 <- data.frame(fread("~/work_space/4.ProjectET/analysis/figs3/blue_genebody.bed.enhancer.bed.CAST.snpcount",header = F))$V4
t3 <- data.frame(fread("~/work_space/4.ProjectET/analysis/figs3/grey_genebody.bed.enhancer.bed.CAST.snpcount",header = F))$V4

wilcox.test(t1,t2) #7.558e-06
wilcox.test(t1,t3) #0.3642

plt.tmp <- data.frame(value=c(t1,t2,t3),tag=c(rep("red",length(t1)),rep("blue",length(t2)),rep("grey",length(t3))))
p <- ggboxplot(plt.tmp,"tag","value",width = 0.4)
ggpar(p,ylim=c(0,20))
ggsave("/home1/gyang/work_space/4.ProjectET/analysis/figs3/backcross_3groups_enhancer_CASTsnp_boxplot.pdf",width = 2.5,height = 3.5)
#all_fpkm2$highlight <- ifelse(all_fpkm2$name %in% key_influ, "yes", "no")
#all_fpkm2 <- spread(all_fpkm2,"class","fpkm")
p <- ggscatter(back_cross, x = "t1", y = "t2",alpha = 1,color = "color",palette = c("#F16A6B","#4682B4","grey60"),
               font.label= c(12, "black")) +
  geom_hline(linetype = "dashed",yintercept = 1.5,color = "grey") +
  geom_hline(linetype = "dashed",yintercept = -1.5,color = "grey") +
  geom_vline(linetype = "dashed",xintercept = 1.5,color = "grey") +
  geom_vline(linetype = "dashed",xintercept = -1.5,color = "grey") +
  xlim(-8.5,8.5) + ylim(-8.5,8.5)+coord_cartesian(expand = F)
p <- ggpar(p,legend = "none")
ggsave("/home1/gyang/work_space/4.ProjectET/analysis/figs3/Backcross_RNA_Scattergroup.pdf",width = 4,height = 4)

############# Unknown remaining: De novo
asg <- data.frame(fread("/home1/gyang/work_space/Bio-Resource/embryos/RNA/all_PWK_snp.csv"),check.names = F)[,c(1,10:21)]
unknown_gene_tro <- data.frame(fread("~/work_space/4.ProjectET/analysis/figs3/epi_at_sight/de_novo_tropho.gene",header = F))$V1
unknown_gene_tro <- asg[asg$name %in% unknown_gene_tro,] 
unknown_gene_tro <- unknown_gene_tro[unknown_gene_tro$name %in% unknown_gene,]
unknown_gene_tro <- merge(unknown_gene_tro, gene_expression, by="name")
plt.tmp <- data.frame(name=unknown_gene_tro[,"name"])
str_all<-  c()
for(i in 1:7){
  index=2*i+1
  str=gsub(".pwk","",colnames(unknown_gene_tro)[index])
  str_all <- c(str_all,str)
  plt.tmp[,i+1] <- unknown_gene_tro[,index-1]- unknown_gene_tro[,index]
}
colnames(plt.tmp)[2:8] <- str_all
plt.tmp[,2:8] <- abs(plt.tmp[,2:8])
boxplot(plt.tmp[,2:8],las=2)

