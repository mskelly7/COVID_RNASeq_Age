# BRAVE Kids RNA Sequencing Analysis
# Aditya Mohan (MD/PhD candidate)  / Matthew Kelly, MD, MPH 
# Figure 2
# Last update: July 26, 2023

remove(list=ls())
setwd("______________________") 
set.seed(2341)
version

if(any(grepl("package:plyr", search()))) detach("package:plyr") else message("plyr not loaded")
library(readr)
library(dplyr)
library(readxl)
library(writexl)
library(ggplot2)
library(tidyr) 
library(ggpubr)
library(plyr) 
library(gridExtra) 
library(cowplot) 
library(forcats) 
library(ggrepel)
options(ggrepel.max.overlaps = Inf)
library(org.Hs.eg.db)
library(tibble)
library(msigdbr)
library(reshape2)
library(data.table)

# Specify expressions for figures
label_log2FC <- expression(bold(log[2] ~ fold ~ change))
label_adjp <- expression(bold(-log[10] ~ p[adj]))

# ANALYSES ADJUSTING FOR CIBERSORT CELL POPULATIONS

# Upload DESeq2 output files
dds_np_Pos_Neg <- read.csv("2_COVID_Pos_vs_Neg/genes_np_Pos_Neg.csv")
dds_pax_Pos_Neg <- read.csv("2_COVID_Pos_vs_Neg/genes_pax_Pos_Neg.csv")
# Upload FGSEA output files
fgsea_np_Pos_Neg <- data.frame(read_excel("2_COVID_Pos_vs_Neg/modules_np_Pos_Neg.xlsx"))
fgsea_pax_Pos_Neg <- data.frame(read_excel("2_COVID_Pos_vs_Neg/modules_pax_Pos_Neg.xlsx"))

# Volcano plots of differentially expressed genes

dds_np_Pos_Neg <- dds_np_Pos_Neg[order(dds_np_Pos_Neg$padj),]
top_genes <- bind_rows(dds_np_Pos_Neg %>% filter(Expression == 'Up-regulated') %>% arrange(desc(abs(log2FoldChange))) %>% head(15),
                       dds_np_Pos_Neg %>% filter(Expression == 'Down-regulated') %>% arrange(desc(abs(log2FoldChange))) %>% head(5))
volcano_np_Pos_Neg <- ggplot(dds_np_Pos_Neg, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(color = Expression), size = 1, alpha=0.8) + theme_bw() +
  xlab(label_log2FC) + ylab(label_adjp) + 
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_hline(yintercept=1.3, linetype="dashed", color = "black", linewidth=0.5) +
  geom_vline(xintercept = -0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  geom_vline(xintercept = 0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  scale_color_manual(name="Legend", labels = c("Up-regulated"="Up-regulated", "Down-regulated"="Down-regulated", "Unchanged"="Unchanged"), 
                     values = c("Up-regulated"="firebrick3", "Down-regulated"="blue", "Unchanged"="gray50")) +
  theme(legend.position = "none", axis.title.x = element_text(size=13, face="bold", margin = ggplot2::margin(t = 5, r = 10, b = 0, l = 0)),
        axis.text.x=element_text(size=11), axis.text.y=element_text(size=11), plot.title=element_blank(),
        axis.title.y = element_text(size=13, face="bold", margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0))) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange<0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(NA,-1)) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange>0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(2,NA)) +
  scale_x_continuous(limits=c(-6,6), breaks=c(-6,-4,-2,0,2,4,6)) +
  scale_y_continuous(limits=c(0,50), breaks=c(0,10,20,30,40,50))
volcano_np_Pos_Neg

dds_pax_Pos_Neg <- dds_pax_Pos_Neg[order(dds_pax_Pos_Neg$padj),]
top_genes <- bind_rows(dds_pax_Pos_Neg %>% filter(Expression == 'Up-regulated') %>% arrange(desc(abs(log2FoldChange))) %>% head(15),
                       dds_pax_Pos_Neg %>% filter(Expression == 'Down-regulated') %>% arrange(desc(abs(log2FoldChange))) %>% head(5))
volcano_pax_Pos_Neg <- ggplot(dds_pax_Pos_Neg, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(color = Expression), size = 1, alpha=0.8) + theme_bw() +
  xlab(label_log2FC) + ylab(label_adjp) + 
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_hline(yintercept=1.3, linetype="dashed", color = "black", linewidth=0.5) +
  geom_vline(xintercept = -0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  geom_vline(xintercept = 0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  scale_color_manual(name="Legend", labels = c("Up-regulated"="Up-regulated", "Down-regulated"="Down-regulated", "Unchanged"="Unchanged"), 
                     values = c("Up-regulated"="firebrick3", "Down-regulated"="blue", "Unchanged"="gray50")) +
  theme(legend.position = "none", axis.title.x = element_text(size=13, face="bold", margin = ggplot2::margin(t = 5, r = 10, b = 0, l = 0)),
        axis.text.x=element_text(size=11), axis.text.y=element_text(size=11), plot.title=element_blank(),
        axis.title.y = element_text(size=13, face="bold", margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0))) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange<0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(NA,-1)) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange>0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(2,NA)) +
  scale_x_continuous(limits=c(-6,6), breaks=c(-6,-4,-2,0,2,4,6)) +
  scale_y_continuous(limits=c(0,50), breaks=c(0,10,20,30,40,50))
volcano_pax_Pos_Neg

# Plots of differentially expressed modules

top_pathways <- bind_rows(fgsea_np_Pos_Neg %>% filter(Expression == 'Up-regulated') %>% arrange(desc(abs(NES))) %>% head(10),
  fgsea_np_Pos_Neg %>% filter(Expression == 'Down-regulated') %>% arrange(desc(abs(NES))) %>% head(10))
module_np_Pos_Neg <- ggplot(fgsea_np_Pos_Neg, aes(x=NES, y=-log10(padj))) + theme_bw() +
  geom_point(aes(color = Expression), size = (fgsea_np_Pos_Neg$GeneRatio*4)^2, alpha=0.4)+
  xlab("Normalized enrichment score") + ylab(label_adjp) +  
  scale_color_manual(values = c( "gray50", "firebrick")) +
  guides(colour = guide_legend(override.aes = list(size=1.5)))+
  geom_hline(yintercept=1.3, linetype="dashed", color = "black", linewidth=0.5) +
  geom_vline(xintercept = -1, linetype="dotdash", color = "black", linewidth=0.5) +
  geom_vline(xintercept = 1, linetype="dotdash", color = "black", linewidth=0.5) +
  theme(axis.text.x=element_text(size=11), axis.text.y=element_text(size=11), 
        axis.title.x = element_text(size=13, face="bold", margin = ggplot2::margin(t = 5, r = 10, b = 0, l = 0)), plot.title = element_blank(),
        axis.title.y = element_text(size=13, face="bold", margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)), legend.position = "none") + 
  geom_label_repel(data = top_pathways, mapping = aes(NES, -log10(padj), label=pathway), size = 4.0, xlim=c(NA,3)) +
  scale_x_continuous(limits=c(-2,4), breaks=c(-2,-1,0,1,2,3,4)) +
  scale_y_continuous(limits=c(0,50), breaks=c(0,10,20,30,40,50))
module_np_Pos_Neg

top_pathways <- bind_rows(fgsea_pax_Pos_Neg %>% filter(Expression == 'Up-regulated') %>% arrange(desc(abs(NES))) %>% head(10),
                          fgsea_pax_Pos_Neg %>% filter(Expression == 'Down-regulated') %>% arrange(desc(abs(NES))) %>% head(10))
module_pax_Pos_Neg <- ggplot(fgsea_pax_Pos_Neg, aes(x=NES, y=-log10(padj))) + theme_bw() +
  geom_point(aes(color = Expression), size = (fgsea_pax_Pos_Neg$GeneRatio*4)^2, alpha=0.4)+
  xlab("Normalized enrichment score") + ylab(label_adjp) + ggtitle("Peripheral blood") +
  scale_color_manual(values = c( "gray50", "firebrick")) +
  guides(colour = guide_legend(override.aes = list(size=1.5)))+
  geom_hline(yintercept=1.3, linetype="dashed", color = "black", linewidth=0.5) +
  geom_vline(xintercept = -1, linetype="dotdash", color = "black", linewidth=0.5) +
  geom_vline(xintercept = 1, linetype="dotdash", color = "black", linewidth=0.5) +
  theme(axis.text.x=element_text(size=11), axis.text.y=element_text(size=11), 
        axis.title.x = element_text(size=13, face="bold", margin = ggplot2::margin(t = 5, r = 10, b = 0, l = 0)), plot.title = element_blank(),
        axis.title.y = element_text(size=13, face="bold", margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)), legend.position = "none") +   
  geom_label_repel(data = top_pathways, mapping = aes(NES, -log10(padj), label=pathway), size = 4.0, xlim=c(NA,3)) +
  scale_x_continuous(limits=c(-2,4), breaks=c(-2,-1,0,1,2,3,4)) +
  scale_y_continuous(limits=c(0,50), breaks=c(0,10,20,30,40,50))
module_pax_Pos_Neg

title_np <- ggdraw() + draw_label("Upper respiratory tract", size=14, fontface='bold')
plot_np <- plot_grid(volcano_np_Pos_Neg, module_np_Pos_Neg, labels=NULL, nrow=1, align="h", rel_widths=c(1,1))
title_pax <- ggdraw() + draw_label("Peripheral blood", size=14, fontface='bold')
plot_pax <- plot_grid(volcano_pax_Pos_Neg, module_pax_Pos_Neg, labels=NULL, nrow=1, align="h", rel_widths=c(1,1))

png(file="Figures/Figure_2.png", width = 16, height = 10, units = 'in', res = 1200)
plot_grid(title_np, plot_np, title_pax, plot_pax, labels=c("a","","b",""), rel_heights=c(0.1,1,0.1,1), nrow=4, align="v") 
dev.off()