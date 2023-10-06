# BRAVE Kids RNA Sequencing Analysis
# Aditya Mohan (MD/PhD candidate)  / Matthew Kelly, MD, MPH 
# Figure S1 - volcano plots among COVID-negative by age
# Last update: Sept. 18, 2023

remove(list=ls())
setwd("G:/My Drive/Research/BRAVE_Kids/RNA_Sequencing/") 
set.seed(2222)
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
dds_np_Neg_0to5_6to13 <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_np_Neg_0to5_6to13.csv")
dds_np_Neg_0to5_14to20 <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_np_Neg_0to5_14to20.csv")
dds_np_Neg_6to13_14to20 <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_np_Neg_6to13_14to20.csv")
dds_pax_Neg_0to5_6to13 <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_pax_Neg_0to5_6to13.csv")
dds_pax_Neg_0to5_14to20 <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_pax_Neg_0to5_14to20.csv")
dds_pax_Neg_6to13_14to20 <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_pax_Neg_6to13_14to20.csv")
dds_pax_Neg_0to5_Adult <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_pax_Neg_0to5_Adult.csv")
dds_pax_Neg_6to13_Adult <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_pax_Neg_6to13_Adult.csv")
dds_pax_Neg_14to20_Adult <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_pax_Neg_14to20_Adult.csv")

top <- 10

top_genes <- bind_rows(dds_np_Neg_0to5_6to13 %>% filter(Expression == 'Up-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top),
                       dds_np_Neg_0to5_6to13 %>% filter(Expression == 'Down-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top))
volcano_np_Neg_0to5_6to13 <- ggplot(dds_np_Neg_0to5_6to13, aes(x=log2FoldChange, y=-log10(padj))) + theme_bw() +
  geom_point(aes(color = Expression), size = 1, alpha=0.8) + xlab(label_log2FC) + ylab(label_adjp) + ggtitle("0-5 vs. 6-13 years") +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_hline(yintercept=1.3, linetype="dashed", color = "black", linewidth=0.5) +
  geom_vline(xintercept = -0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  geom_vline(xintercept = 0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain"), legend.position = "none") + 
  scale_color_manual(name="Legend", labels = c("Up-regulated"="Up-regulated", "Down-regulated"="Down-regulated", "Unchanged"="Unchanged"), 
                     values = c("Up-regulated"="firebrick3", "Down-regulated"="blue", "Unchanged"="gray50")) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange<0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(NA,-4)) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange>0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(3,NA)) +
  scale_x_continuous(limits=c(-15,15), breaks=c(-15,-10,-5,0,5,10,15)) +
  scale_y_continuous(limits=c(0,8), breaks=c(0,2,4,6,8))
volcano_np_Neg_0to5_6to13

top_genes <- bind_rows(dds_np_Neg_0to5_14to20 %>% filter(Expression == 'Up-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top),
                       dds_np_Neg_0to5_14to20 %>% filter(Expression == 'Down-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top))
volcano_np_Neg_0to5_14to20 <- ggplot(dds_np_Neg_0to5_14to20, aes(x=log2FoldChange, y=-log10(padj))) + theme_bw() +
  geom_point(aes(color = Expression), size = 1, alpha=0.8) + xlab(label_log2FC) + ylab("") + ggtitle("0-5 vs. 14-20 years") +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_hline(yintercept=1.3, linetype="dashed", color = "black", linewidth=0.5) +
  geom_vline(xintercept = -0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  geom_vline(xintercept = 0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain"), legend.position = "none") + 
  scale_color_manual(name="Legend", labels = c("Up-regulated"="Up-regulated", "Down-regulated"="Down-regulated", "Unchanged"="Unchanged"), 
                     values = c("Up-regulated"="firebrick3", "Down-regulated"="blue", "Unchanged"="gray50")) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange<0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(NA,-4)) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange>0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(2,NA)) +
  scale_x_continuous(limits=c(-15,15), breaks=c(-15,-10,-5,0,5,10,15)) +
  scale_y_continuous(limits=c(0,8), breaks=c(0,2,4,6,8))
volcano_np_Neg_0to5_14to20

top_genes <- bind_rows(dds_np_Neg_6to13_14to20 %>% filter(Expression == 'Up-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top),
                       dds_np_Neg_6to13_14to20 %>% filter(Expression == 'Down-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top))
volcano_np_Neg_6to13_14to20 <- ggplot(dds_np_Neg_6to13_14to20, aes(x=log2FoldChange, y=-log10(padj))) + theme_bw() +
  geom_point(aes(color = Expression), size = 1, alpha=0.8) + xlab(label_log2FC) + ylab("") + ggtitle("6-13 vs. 14-20 years") +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_hline(yintercept=1.3, linetype="dashed", color = "black", linewidth=0.5) +
  geom_vline(xintercept = -0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  geom_vline(xintercept = 0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain"), legend.position = "none") + 
  scale_color_manual(name="Legend", labels = c("Up-regulated"="Up-regulated", "Down-regulated"="Down-regulated", "Unchanged"="Unchanged"), 
                     values = c("Up-regulated"="firebrick3", "Down-regulated"="blue", "Unchanged"="gray50")) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange<0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(NA,-4)) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange>0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(3,NA)) +
  scale_x_continuous(limits=c(-15,15), breaks=c(-15,-10,-5,0,5,10,15)) +
  scale_y_continuous(limits=c(0,8), breaks=c(0,2,4,6,8))
volcano_np_Neg_6to13_14to20

top_genes <- bind_rows(dds_pax_Neg_0to5_6to13 %>% filter(Expression == 'Up-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top),
                       dds_pax_Neg_0to5_6to13 %>% filter(Expression == 'Down-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top))
volcano_pax_Neg_0to5_6to13 <- ggplot(dds_pax_Neg_0to5_6to13, aes(x=log2FoldChange, y=-log10(padj))) + theme_bw() +
  geom_point(aes(color = Expression), size = 1, alpha=0.8) + xlab("") + ylab(label_adjp) + ggtitle("0-5 vs. 6-13 years") +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_hline(yintercept=1.3, linetype="dashed", color = "black", linewidth=0.5) +
  geom_vline(xintercept = -0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  geom_vline(xintercept = 0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain"), legend.position = "none") + 
  scale_color_manual(name="Legend", labels = c("Up-regulated"="Up-regulated", "Down-regulated"="Down-regulated", "Unchanged"="Unchanged"), 
                     values = c("Up-regulated"="firebrick3", "Down-regulated"="blue", "Unchanged"="gray50")) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange<0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(NA,-4)) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange>0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(3,NA)) +
  scale_x_continuous(limits=c(-15,15), breaks=c(-15,-10,-5,0,5,10,15)) +
  scale_y_continuous(limits=c(0,8), breaks=c(0,2,4,6,8))
volcano_pax_Neg_0to5_6to13

top_genes <- bind_rows(dds_pax_Neg_0to5_14to20 %>% filter(Expression == 'Up-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top),
                       dds_pax_Neg_0to5_14to20 %>% filter(Expression == 'Down-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top))
volcano_pax_Neg_0to5_14to20 <- ggplot(dds_pax_Neg_0to5_14to20, aes(x=log2FoldChange, y=-log10(padj))) + theme_bw() +
  geom_point(aes(color = Expression), size = 1, alpha=0.8) + xlab("") + ylab("") + ggtitle("0-5 vs. 14-20 years") +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_hline(yintercept=1.3, linetype="dashed", color = "black", linewidth=0.5) +
  geom_vline(xintercept = -0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  geom_vline(xintercept = 0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain"), legend.position = "none") + 
  scale_color_manual(name="Legend", labels = c("Up-regulated"="Up-regulated", "Down-regulated"="Down-regulated", "Unchanged"="Unchanged"), 
                     values = c("Up-regulated"="firebrick3", "Down-regulated"="blue", "Unchanged"="gray50")) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange<0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(NA,-4)) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange>0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(3,NA)) +
  scale_x_continuous(limits=c(-15,15), breaks=c(-15,-10,-5,0,5,10,15)) +
  scale_y_continuous(limits=c(0,8), breaks=c(0,2,4,6,8))
volcano_pax_Neg_0to5_14to20

top_genes <- bind_rows(dds_pax_Neg_6to13_14to20 %>% filter(Expression == 'Up-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top),
                       dds_pax_Neg_6to13_14to20 %>% filter(Expression == 'Down-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top))
volcano_pax_Neg_6to13_14to20 <- ggplot(dds_pax_Neg_6to13_14to20, aes(x=log2FoldChange, y=-log10(padj))) + theme_bw() +
  geom_point(aes(color = Expression), size = 1, alpha=0.8) + xlab("") + ylab("") + ggtitle("6-13 vs. 14-20 years") +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_hline(yintercept=1.3, linetype="dashed", color = "black", linewidth=0.5) +
  geom_vline(xintercept = -0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  geom_vline(xintercept = 0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain"), legend.position = "none") + 
  scale_color_manual(name="Legend", labels = c("Up-regulated"="Up-regulated", "Down-regulated"="Down-regulated", "Unchanged"="Unchanged"), 
                     values = c("Up-regulated"="firebrick3", "Down-regulated"="blue", "Unchanged"="gray50")) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange<0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(NA,-4)) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange>0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(3,NA)) +
  scale_x_continuous(limits=c(-15,15), breaks=c(-15,-10,-5,0,5,10,15)) +
  scale_y_continuous(limits=c(0,8), breaks=c(0,2,4,6,8))
volcano_pax_Neg_6to13_14to20

top_genes <- bind_rows(dds_pax_Neg_0to5_Adult %>% filter(Expression == 'Up-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top),
                       dds_pax_Neg_0to5_Adult %>% filter(Expression == 'Down-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top))
volcano_pax_Neg_0to5_Adult <- ggplot(dds_pax_Neg_0to5_Adult, aes(x=log2FoldChange, y=-log10(padj))) + theme_bw() +
  geom_point(aes(color = Expression), size = 1, alpha=0.8) + xlab(label_log2FC) + ylab(label_adjp) + ggtitle("0-5 years vs. Adult") +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_hline(yintercept=1.3, linetype="dashed", color = "black", linewidth=0.5) +
  geom_vline(xintercept = -0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  geom_vline(xintercept = 0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain"), legend.position = "none") + 
  scale_color_manual(name="Legend", labels = c("Up-regulated"="Up-regulated", "Down-regulated"="Down-regulated", "Unchanged"="Unchanged"), 
                     values = c("Up-regulated"="firebrick3", "Down-regulated"="blue", "Unchanged"="gray50")) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange<0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(NA,-4)) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange>0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(3,NA)) +
  scale_x_continuous(limits=c(-15,15), breaks=c(-15,-10,-5,0,5,10,15)) +
  scale_y_continuous(limits=c(0,8), breaks=c(0,2,4,6,8))
volcano_pax_Neg_0to5_Adult

top_genes <- bind_rows(dds_pax_Neg_6to13_Adult %>% filter(Expression == 'Up-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top),
                       dds_pax_Neg_6to13_Adult %>% filter(Expression == 'Down-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top))
volcano_pax_Neg_6to13_Adult <- ggplot(dds_pax_Neg_6to13_Adult, aes(x=log2FoldChange, y=-log10(padj))) + theme_bw() +
  geom_point(aes(color = Expression), size = 1, alpha=0.8) + xlab(label_log2FC) + ylab("") + ggtitle("6-13 years vs. Adult") +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_hline(yintercept=1.3, linetype="dashed", color = "black", linewidth=0.5) +
  geom_vline(xintercept = -0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  geom_vline(xintercept = 0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain"), legend.position = "none") + 
  scale_color_manual(name="Legend", labels = c("Up-regulated"="Up-regulated", "Down-regulated"="Down-regulated", "Unchanged"="Unchanged"), 
                     values = c("Up-regulated"="firebrick3", "Down-regulated"="blue", "Unchanged"="gray50")) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange<0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(NA,-4)) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange>0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(3,NA)) +
  scale_x_continuous(limits=c(-15,15), breaks=c(-15,-10,-5,0,5,10,15)) +
  scale_y_continuous(limits=c(0,8), breaks=c(0,2,4,6,8))
volcano_pax_Neg_6to13_Adult

top_genes <- bind_rows(dds_pax_Neg_14to20_Adult %>% filter(Expression == 'Up-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top),
                       dds_pax_Neg_14to20_Adult %>% filter(Expression == 'Down-regulated') %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(top))
volcano_pax_Neg_14to20_Adult <- ggplot(dds_pax_Neg_14to20_Adult, aes(x=log2FoldChange, y=-log10(padj))) + theme_bw() +
  geom_point(aes(color = Expression), size = 1, alpha=0.8) + xlab(label_log2FC) + ylab("") + ggtitle("14-20 years vs. Adult") +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_hline(yintercept=1.3, linetype="dashed", color = "black", linewidth=0.5) +
  geom_vline(xintercept = -0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  geom_vline(xintercept = 0.5, linetype="dotdash", color = "black", linewidth=0.5) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain"), legend.position = "none") + 
  scale_color_manual(name="Legend", labels = c("Up-regulated"="Up-regulated", "Down-regulated"="Down-regulated", "Unchanged"="Unchanged"), 
                     values = c("Up-regulated"="firebrick3", "Down-regulated"="blue", "Unchanged"="gray50")) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange<0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(NA,-4)) +
  geom_label_repel(data = top_genes[top_genes$log2FoldChange>0,], mapping = aes(log2FoldChange, -log10(padj), label = Gene), size = 4.0, xlim=c(3,NA)) +
  scale_x_continuous(limits=c(-15,15), breaks=c(-15,-10,-5,0,5,10,15)) +
  scale_y_continuous(limits=c(0,8), breaks=c(0,2,4,6,8))
volcano_pax_Neg_14to20_Adult

title_np <- ggdraw() + draw_label("     Upper respiratory", size=16, fontface='bold')
volcano_np <- plot_grid(volcano_np_Neg_0to5_6to13, volcano_np_Neg_0to5_14to20, volcano_np_Neg_6to13_14to20, labels=NULL, rel_widths=c(1,0.99,0.99), nrow=1, align="h")
title_pax <- ggdraw() + draw_label("     Peripheral blood", size=16, fontface='bold')
volcano_pax_peds <- plot_grid(volcano_pax_Neg_0to5_6to13, volcano_pax_Neg_0to5_14to20, volcano_pax_Neg_6to13_14to20, labels=NULL, rel_widths=c(1,0.99,0.99), nrow=1, align="h")
volcano_pax_adult <- plot_grid(volcano_pax_Neg_0to5_Adult, volcano_pax_Neg_6to13_Adult, volcano_pax_Neg_14to20_Adult, labels=NULL, rel_widths=c(1,0.99,0.99), nrow=1, align="h")

# NO TITLE
#set.seed(1234)
#png(file="Statistical_Analyses/Figures/Figure_S1.png", width = 18, height = 12, units = 'in', res = 1200)
#plot_grid(title_np, volcano_np, title_pax, volcano_pax_peds, volcano_pax_adult, rel_heights=c(0.15,1,0.15,1,1), nrow=5, align="v", labels=c("a","","b","","")) 
#dev.off()

# ADDITION OF TITLE

figS1_plot <- plot_grid(title_np, volcano_np, title_pax, volcano_pax_peds, volcano_pax_adult, rel_heights=c(0.15,1,0.15,1,1), nrow=5, align="v", labels=c("a","","b","","")) 
figS1_title <- ggdraw() + draw_label("Figure S1. Differential host gene expression in upper respiratory and peripheral blood samples associated with SARS-CoV-2 infection among children, adolescents, and adults", 
                                    size=14, fontface='bold')

set.seed(1234)
png(file="Statistical_Analyses/Figures/Figure_S1.png", width = 16.5, height = 12, units = 'in', res = 1200)
plot_grid(figS1_title, NULL, figS1_plot, labels=NULL, ncol=1, rel_heights=c(0.03,0.02,0.97), align="v") 
dev.off()