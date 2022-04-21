setwd("~/work_space/4.ProjectET/analysis/figs3")

library(dplyr)
library(ggpubr)
library(stringr)
library(tidyr)

wt_expr <- data.frame(fread("~/work_space/Bio-Resource/embryos/RNA/Normal_embryo_noem2c_log2FPKM.csv"),check.names = F)
wt_expr_xw <- data.frame(fread("~/work_space/Bio-Resource/embryos/RNA/Xiewei_2016_Nature_log2FPKM_RNA.csv",check.names=FALSE),check.names = F)
wt_expr <- wt_expr[,c("gene","Morula")]
wt_expr_all <- merge(wt_expr, wt_expr_xw, by.x = "gene", by.y = "Gene")
wt_expr_all <- wt_expr_all[,c("gene","MII.oocyte","PN5.zygote","Early.2.cell","Late.2.cell","4.cell","8.cell" ,"ICM" , "mESC" )]
rownames(wt_expr_all) <- wt_expr_all$gene
wt_klf_family <- wt_expr_all[str_which(wt_expr_all$gene, "Klf"), 1:ncol(wt_expr_all)]

#wt_klf_family[,2:ncol(wt_klf_family)] <- t(scale(t(wt_klf_family[,2:ncol(wt_klf_family)])))


wt_klf_family <- gather(wt_klf_family,"stage", "rpkm", -"gene")
wt_klf_family$color <- ifelse(wt_klf_family$gene == "Klf5", "colored", "grey")
wt_klf_family$stage <-  factor(wt_klf_family$stage ,levels=c("MII.oocyte","PN5.zygote","Early.2.cell","Late.2.cell","4.cell","8.cell" ,"ICM" , "mESC"))
p <- ggline(wt_klf_family, "stage", "rpkm",  color = "gene", palette =  c("grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","IndianRed1","grey","grey","grey","grey"))

# ggplot(wt_klf_family) +
#   geom_line(aes(stage,rpkm,color=color,group=gene),)

ggpar(p,x.text.angle = 90)
ggsave("Klf_family_xieweiRNA_log2FPKM.pdf",dpi=300,width=4,height=7)
