# BRAVE Kids RNA Sequencing Analysis
# Aditya Mohan (MD/PhD candidate)  / Matthew Kelly, MD, MPH 
# Figure S3 - comparison of upper respiratory gene expression among SARS-CoV-2-infected vs. uninfected children and adolescents by age
# Last update: June 8, 2024

remove(list=ls())
setwd("____________________________") 
set.seed(123456)
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
dds_np_0to5_pos_neg_nocibersort <- read.csv("Statistical_Analyses/3_COVID_Pos_by_Age/genes_np_0to5_pos_neg_nocibersort.csv")
dds_np_6to13_pos_neg_nocibersort <- read.csv("Statistical_Analyses/3_COVID_Pos_by_Age/genes_np_6to13_pos_neg_nocibersort.csv")
dds_np_14to20_pos_neg_nocibersort <- read.csv("Statistical_Analyses/3_COVID_Pos_by_Age/genes_np_14to20_pos_neg_nocibersort.csv")

top <- 10

top_genes <- bind_rows(dds_np_0to5_pos_neg_nocibersort %>% filter(Expression == 'Up-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top),
                       dds_np_0to5_pos_neg_nocibersort %>% filter(Expression == 'Down-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top))
volcano_np_0to5_Pos_Neg <- ggplot(dds_np_0to5_pos_neg_nocibersort, aes(x=log2FoldChange, y=-log10(padj))) + theme_bw() +
  geom_point(aes(color = Expression), size = 1, alpha=0.8) + xlab(label_log2FC) + ylab(label_adjp) + ggtitle("0-5 years") +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_hline(yintercept=1.3, linetype="dashed", color = "black", linewidth=0.5) +
  geom_vline(xintercept = -0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  geom_vline(xintercept = 0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain"), legend.position = "none") + 
  scale_color_manual(name="Legend", labels = c("Up-regulated"="Up-regulated", "Down-regulated"="Down-regulated", "Unchanged"="Unchanged"), 
                     values = c("Up-regulated"="firebrick3", "Down-regulated"="blue", "Unchanged"="gray50")) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange<0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(NA,-2)) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange>0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(3,NA)) +
  scale_x_continuous(limits=c(-10,10), breaks=c(-10,-8,-6,-4,-2,0,2,4,6,8,10)) +
  scale_y_continuous(limits=c(0,10.5), breaks=c(0,2,4,6,8,10,12))
volcano_np_0to5_Pos_Neg

top_genes <- bind_rows(dds_np_6to13_pos_neg_nocibersort %>% filter(Expression == 'Up-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top),
                       dds_np_6to13_pos_neg_nocibersort %>% filter(Expression == 'Down-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top))
volcano_np_6to13_Pos_Neg <- ggplot(dds_np_6to13_pos_neg_nocibersort, aes(x=log2FoldChange, y=-log10(padj))) + theme_bw() +
  geom_point(aes(color = Expression), size = 1, alpha=0.8) + xlab(label_log2FC) + ylab("") + ggtitle("6-13 years") +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_hline(yintercept=1.3, linetype="dashed", color = "black", linewidth=0.5) +
  geom_vline(xintercept = -0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  geom_vline(xintercept = 0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain"), legend.position = "none") + 
  scale_color_manual(name="Legend", labels = c("Up-regulated"="Up-regulated", "Down-regulated"="Down-regulated", "Unchanged"="Unchanged"), 
                     values = c("Up-regulated"="firebrick3", "Down-regulated"="blue", "Unchanged"="gray50")) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange<0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(NA,-2)) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange>0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(3,NA)) +
  scale_x_continuous(limits=c(-10,10), breaks=c(-10,-8,-6,-4,-2,0,2,4,6,8,10)) +
  scale_y_continuous(limits=c(0,10.5), breaks=c(0,2,4,6,8,10,12))
volcano_np_6to13_Pos_Neg

top_genes <- bind_rows(dds_np_14to20_pos_neg_nocibersort %>% filter(Expression == 'Up-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top),
                       dds_np_14to20_pos_neg_nocibersort %>% filter(Expression == 'Down-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top))
volcano_np_14to20_Pos_Neg <- ggplot(dds_np_14to20_pos_neg_nocibersort, aes(x=log2FoldChange, y=-log10(padj))) + theme_bw() +
  geom_point(aes(color = Expression), size = 1, alpha=0.8) + xlab(label_log2FC) + ylab("") + ggtitle("14-20 years") +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_hline(yintercept=1.3, linetype="dashed", color = "black", linewidth=0.5) +
  geom_vline(xintercept = -0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  geom_vline(xintercept = 0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain"), legend.position = "none") + 
  scale_color_manual(name="Legend", labels = c("Up-regulated"="Up-regulated", "Down-regulated"="Down-regulated", "Unchanged"="Unchanged"), 
                     values = c("Up-regulated"="firebrick3", "Down-regulated"="blue", "Unchanged"="gray50")) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange<0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(NA,-2)) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange>0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(3,NA)) +
  scale_x_continuous(limits=c(-10,10), breaks=c(-10,-8,-6,-4,-2,0,2,4,6,8,10)) +
  scale_y_continuous(limits=c(0,10.5), breaks=c(0,2,4,6,8,10,12))
volcano_np_14to20_Pos_Neg

volcano_np <- plot_grid(volcano_np_0to5_Pos_Neg, volcano_np_6to13_Pos_Neg, volcano_np_14to20_Pos_Neg, rel_widths=c(1,0.99,0.99), labels=NULL, nrow=1, align="h")

set.seed(1234)
png(file="Statistical_Analyses/Figures/Figure_S3.png", width = 21, height = 4.5, units = 'in', res = 1200)
plot_grid(volcano_np) 
dev.off()