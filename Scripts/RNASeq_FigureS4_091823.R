# BRAVE Kids RNA Sequencing Analysis
# Aditya Mohan (MD/PhD candidate)  / Matthew Kelly, MD, MPH 
# Figure S4 - volcano plots for COVID+ vs. COVID- by age category (Peripheral blood)
# Last update: Sept. 18, 2023

remove(list=ls())
setwd("G:/My Drive/Research/BRAVE_Kids/RNA_Sequencing/Statistical_Analyses/") 
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
dds_pax_0to5_Pos_Neg <- read.csv("3_COVID_Pos_by_Age/genes_pax_0to5_Pos_Neg.csv")
dds_pax_6to13_Pos_Neg <- read.csv("3_COVID_Pos_by_Age/genes_pax_6to13_Pos_Neg.csv")
dds_pax_14to20_Pos_Neg <- read.csv("3_COVID_Pos_by_Age/genes_pax_14to20_Pos_Neg.csv")
dds_pax_adult_Pos_Neg <- read.csv("3_COVID_Pos_by_Age/genes_pax_adult_Pos_Neg.csv")

top <- 10

top_genes <- bind_rows(dds_pax_0to5_Pos_Neg %>% filter(Expression == 'Up-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top),
                       dds_pax_0to5_Pos_Neg %>% filter(Expression == 'Down-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top))
volcano_pax_0to5_Pos_Neg <- ggplot(dds_pax_0to5_Pos_Neg, aes(x=log2FoldChange, y=-log10(padj))) + theme_bw() +
  geom_point(aes(color = Expression), size = 1, alpha=0.8) + xlab("") + ylab(label_adjp) + ggtitle("0-5 years") +
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
  scale_y_continuous(limits=c(0,20), breaks=c(0,5,10,15,20))
volcano_pax_0to5_Pos_Neg

top_genes <- bind_rows(dds_pax_6to13_Pos_Neg %>% filter(Expression == 'Up-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top),
                       dds_pax_6to13_Pos_Neg %>% filter(Expression == 'Down-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top))
volcano_pax_6to13_Pos_Neg <- ggplot(dds_pax_6to13_Pos_Neg, aes(x=log2FoldChange, y=-log10(padj))) + theme_bw() +
  geom_point(aes(color = Expression), size = 1, alpha=0.8) + xlab("") + ylab("") + ggtitle("6-13 years") +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_hline(yintercept=1.3, linetype="dashed", color = "black", linewidth=0.5) +
  geom_vline(xintercept = -0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  geom_vline(xintercept = 0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain"), legend.position = "none") + 
  scale_color_manual(name="Legend", labels = c("Up-regulated"="Up-regulated", "Down-regulated"="Down-regulated", "Unchanged"="Unchanged"), 
                     values = c("Up-regulated"="firebrick3", "Down-regulated"="blue", "Unchanged"="gray50")) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange<0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(NA,-3)) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange>0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(2,NA)) +
  scale_x_continuous(limits=c(-10,10), breaks=c(-10,-8,-6,-4,-2,0,2,4,6,8,10)) +
  scale_y_continuous(limits=c(0,20), breaks=c(0,5,10,15,20))
volcano_pax_6to13_Pos_Neg

top_genes <- bind_rows(dds_pax_14to20_Pos_Neg %>% filter(Expression == 'Up-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top),
                       dds_pax_14to20_Pos_Neg %>% filter(Expression == 'Down-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top))
volcano_pax_14to20_Pos_Neg <- ggplot(dds_pax_14to20_Pos_Neg, aes(x=log2FoldChange, y=-log10(padj))) + theme_bw() +
  geom_point(aes(color = Expression), size = 1, alpha=0.8) + xlab(label_log2FC) + ylab(label_adjp) + ggtitle("14-20 years") +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_hline(yintercept=1.3, linetype="dashed", color = "black", linewidth=0.5) +
  geom_vline(xintercept = -0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  geom_vline(xintercept = 0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain"), legend.position = "none") + 
  scale_color_manual(name="Legend", labels = c("Up-regulated"="Up-regulated", "Down-regulated"="Down-regulated", "Unchanged"="Unchanged"), 
                     values = c("Up-regulated"="firebrick3", "Down-regulated"="blue", "Unchanged"="gray50")) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange<0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(NA,-4)) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange>0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(3,NA)) +
  scale_x_continuous(limits=c(-30,30), breaks=c(-30,-20,-10,0,10,20,30)) +
  scale_y_continuous(limits=c(0,20), breaks=c(0,5,10,15,20))
volcano_pax_14to20_Pos_Neg

top_genes <- bind_rows(dds_pax_adult_Pos_Neg %>% filter(Expression == 'Up-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top),
                       dds_pax_adult_Pos_Neg %>% filter(Expression == 'Down-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top))
volcano_pax_adult_Pos_Neg <- ggplot(dds_pax_adult_Pos_Neg, aes(x=log2FoldChange, y=-log10(padj))) + theme_bw() +
  geom_point(aes(color = Expression), size = 1, alpha=0.8) + xlab(label_log2FC) + ylab("") + ggtitle("Adult") +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_hline(yintercept=1.3, linetype="dashed", color = "black", linewidth=0.5) +
  geom_vline(xintercept = -0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  geom_vline(xintercept = 0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain"), legend.position = "none") + 
  scale_color_manual(name="Legend", labels = c("Up-regulated"="Up-regulated", "Down-regulated"="Down-regulated", "Unchanged"="Unchanged"), 
                     values = c("Up-regulated"="firebrick3", "Down-regulated"="blue", "Unchanged"="gray50")) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange<0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(NA,-2.5)) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange>0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(3,NA)) +
  scale_x_continuous(limits=c(-10,10), breaks=c(-10,-8,-6,-4,-2,0,2,4,6,8,10)) +
  scale_y_continuous(limits=c(0,20), breaks=c(0,5,10,15,20))
volcano_pax_adult_Pos_Neg

title_pax <- ggdraw() + draw_label("Figure S4. Differential expression of individual genes in peripheral blood associated with SARS-CoV-2 infection among children, adolescents, and adults", 
                                   size=14, fontface='bold')
volcano_pax_1 <- plot_grid(volcano_pax_0to5_Pos_Neg, volcano_pax_6to13_Pos_Neg, labels=NULL, rel_widths=c(1,0.99), nrow=1, align="h")
volcano_pax_2 <- plot_grid(volcano_pax_14to20_Pos_Neg, volcano_pax_adult_Pos_Neg, labels=NULL, rel_widths=c(1,0.99), nrow=1, align="h")

png(file="Figures/Figure_S4.png", width = 15, height = 10.5, units = 'in', res = 1200)
plot_grid(title_pax, volcano_pax_1, volcano_pax_2, rel_heights=c(0.1,0.99,1), nrow=3, align="vh") 
dev.off()