# BRAVE Kids RNA Sequencing Analysis
# Aditya Mohan (MD/PhD candidate)  / Matthew Kelly, MD, MPH 
# Figure S3 - volcano plots for COVID+ vs. COVID- by age category (Upper respiratory)
# Last update: Sept. 18, 2023

remove(list=ls())
setwd("____________________________") 
set.seed(1234)
version

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
library(tibble)
library(msigdbr)
library(reshape2)
library(patchwork)
library(rstatix)
library(data.table)
library(stringr)

# Specify expressions for figures
label_log2FC <- expression(bold(log[2] ~ fold ~ change))
label_adjp <- expression(bold(-log[10] ~ p[adj]))

# Upload DESeq2 output files
dds_np_0to5_Pos_Neg <- read.csv("Statistical_Analyses/3_COVID_Pos_by_Age/genes_np_0to5_Pos_Neg.csv")
dds_np_6to13_Pos_Neg <- read.csv("Statistical_Analyses/3_COVID_Pos_by_Age/genes_np_6to13_Pos_Neg.csv")
dds_np_14to20_Pos_Neg <- read.csv("Statistical_Analyses/3_COVID_Pos_by_Age/genes_np_14to20_Pos_Neg.csv")

top <- 10

top_genes <- bind_rows(dds_np_0to5_Pos_Neg %>% filter(Expression == 'Up-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top),
                       dds_np_0to5_Pos_Neg %>% filter(Expression == 'Down-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top))
volcano_np_0to5_Pos_Neg <- ggplot(dds_np_0to5_Pos_Neg, aes(x=log2FoldChange, y=-log10(padj))) + theme_bw() +
  geom_point(aes(color = Expression), size = 1, alpha=0.8) + xlab(label_log2FC) + ylab(label_adjp) + ggtitle("0-5 years") +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_hline(yintercept=1.3, linetype="dashed", color = "black", linewidth=0.5) +
  geom_vline(xintercept = -0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  geom_vline(xintercept = 0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain"), legend.position = "none") + 
  scale_color_manual(name="Legend", labels = c("Up-regulated"="Up-regulated", "Down-regulated"="Down-regulated", "Unchanged"="Unchanged"), 
                     values = c("Up-regulated"="firebrick3", "Down-regulated"="blue", "Unchanged"="gray50")) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange<0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(NA,-4)) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange>0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(3,NA)) +
  scale_x_continuous(limits=c(-10,10), breaks=c(-10,-8,-6,-4,-2,0,2,4,6,8,10)) +
  scale_y_continuous(limits=c(0,6), breaks=c(0,1,2,3,4,5,6))
volcano_np_0to5_Pos_Neg

top_genes <- bind_rows(dds_np_6to13_Pos_Neg %>% filter(Expression == 'Up-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top),
                       dds_np_6to13_Pos_Neg %>% filter(Expression == 'Down-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top))
volcano_np_6to13_Pos_Neg <- ggplot(dds_np_6to13_Pos_Neg, aes(x=log2FoldChange, y=-log10(padj))) + theme_bw() +
  geom_point(aes(color = Expression), size = 1, alpha=0.8) + xlab(label_log2FC) + ylab("") + ggtitle("6-13 years") +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_hline(yintercept=1.3, linetype="dashed", color = "black", linewidth=0.5) +
  geom_vline(xintercept = -0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  geom_vline(xintercept = 0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain"), legend.position = "none") + 
  scale_color_manual(name="Legend", labels = c("Up-regulated"="Up-regulated", "Down-regulated"="Down-regulated", "Unchanged"="Unchanged"), 
                     values = c("Up-regulated"="firebrick3", "Down-regulated"="blue", "Unchanged"="gray50")) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange<0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(NA,-3)) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange>0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(3,NA)) +
  scale_x_continuous(limits=c(-10,10), breaks=c(-10,-8,-6,-4,-2,0,2,4,6,8,10)) +
  scale_y_continuous(limits=c(0,6), breaks=c(0,1,2,3,4,5,6))
volcano_np_6to13_Pos_Neg

top_genes <- bind_rows(dds_np_14to20_Pos_Neg %>% filter(Expression == 'Up-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top),
                       dds_np_14to20_Pos_Neg %>% filter(Expression == 'Down-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top))
volcano_np_14to20_Pos_Neg <- ggplot(dds_np_14to20_Pos_Neg, aes(x=log2FoldChange, y=-log10(padj))) + theme_bw() +
  geom_point(aes(color = Expression), size = 1, alpha=0.8) + xlab(label_log2FC) + ylab("") + ggtitle("14-20 years") +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_hline(yintercept=1.3, linetype="dashed", color = "black", linewidth=0.5) +
  geom_vline(xintercept = -0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  geom_vline(xintercept = 0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain"), legend.position = "none") + 
  scale_color_manual(name="Legend", labels = c("Up-regulated"="Up-regulated", "Down-regulated"="Down-regulated", "Unchanged"="Unchanged"), 
                     values = c("Up-regulated"="firebrick3", "Down-regulated"="blue", "Unchanged"="gray50")) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange<0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(NA,-3)) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange>0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(3,NA)) +
  scale_x_continuous(limits=c(-10,10), breaks=c(-10,-8,-6,-4,-2,0,2,4,6,8,10)) +
  scale_y_continuous(limits=c(0,6), breaks=c(0,1,2,3,4,5,6))
volcano_np_14to20_Pos_Neg

title_np <- ggdraw() + draw_label("Figure S3. Differential expression of individual genes in the upper respiratory tract associated with SARS-CoV-2 infection among children and adolescents", 
                                  size=16, fontface='bold')
volcano_np <- plot_grid(volcano_np_0to5_Pos_Neg, volcano_np_6to13_Pos_Neg, volcano_np_14to20_Pos_Neg, rel_widths=c(1,0.99,0.99), labels=NULL, nrow=1, align="h")

set.seed(1234)
png(file="Statistical_Analyses/Figures/Figure_S3.png", width = 20, height = 5, units = 'in', res = 1200)
plot_grid(title_np, volcano_np, rel_heights=c(0.15,1), nrow=2, align="v") 
dev.off()