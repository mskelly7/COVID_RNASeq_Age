# BRAVE Kids RNA Sequencing Analysis
# Aditya Mohan (MD/PhD candidate) / Matthew Kelly, MD, MPH 
# Figure 4 - heatmap of gene module expression by symptom presence in SARS-CoV-2-infected individuals
# Analyses of gene module expression in peripheral blood samples adjusted for imputed cell proportions
# Last update: June 13, 2024

remove(list=ls())
setwd("____________________________") 
set.seed(1234)
version

library(readr)
library(circlize)
library(ComplexHeatmap)
library(data.table)
library(Hmisc)
library(corrplot)
library(readxl)
library(writexl)
library(readr)
library(dplyr)
library(readxl)
library(writexl)
library(ggplot2)
library(gridExtra) 
library(cowplot) 
library(forcats) 
library(tibble)
library(reshape2)
library(phyloseq)
library(ggcorrplot)
library(plyr)
library(matrixStats)
library(stringr)

vars_to_keep <- c("Variable","pathway","NES","padj")
innate_pathways <- c("Innate Immune Cell Activation","Interferon Response","Type I Interferon Signaling","Type II Interferon Signaling","Type III Interferon Signaling", 
                     "TNF Signaling","Chemokine Signaling","TLR Signaling","NLR Signaling","RNA Sensing","Phagocytosis","Myeloid Activation",
                     "Myeloid Inflammation","Inflammasomes","NK Activity","Complement System","Coagulation","Leukotriene and Prostaglandin Inflammation")
adaptive_pathways <- c("Adaptive Immune Response","MHC Class I Antigen Presentation","MHC Class II Antigen Presentation","Mononuclear Cell Migration","Lymphocyte Trafficking",
                       "BCR Signaling","NF-kappaB Signaling","TCR Signaling","JAK-STAT Signaling","Immune Memory","Immune Exhaustion","Treg Differentiation")

# NP Heatmap

# Load the required datasets
fgsea_np_fever_nocibersort <- read_xlsx("Statistical_Analyses/4_Symptoms/modules_np_fever_nocibersort.xlsx")
fgsea_np_fever_nocibersort$Variable <- "Fever"
fgsea_np_fever_nocibersort <- fgsea_np_fever_nocibersort[,vars_to_keep]
fgsea_np_cough_nocibersort <- read_xlsx("Statistical_Analyses/4_Symptoms/modules_np_cough_nocibersort.xlsx")
fgsea_np_cough_nocibersort$Variable <- "Cough"
fgsea_np_cough_nocibersort <- fgsea_np_cough_nocibersort[,vars_to_keep]
fgsea_np_rhinorrhea_nocibersort <- read_xlsx("Statistical_Analyses/4_Symptoms/modules_np_rhinorrhea_nocibersort.xlsx")
fgsea_np_rhinorrhea_nocibersort$Variable <- "Rhinorrhea"
fgsea_np_rhinorrhea_nocibersort <- fgsea_np_rhinorrhea_nocibersort[,vars_to_keep]
fgsea_np_congestion_nocibersort <- read_xlsx("Statistical_Analyses/4_Symptoms/modules_np_congestion_nocibersort.xlsx")
fgsea_np_congestion_nocibersort$Variable <- "Congestion"
fgsea_np_congestion_nocibersort <- fgsea_np_congestion_nocibersort[,vars_to_keep]
fgsea_np_headache_nocibersort <- read_xlsx("Statistical_Analyses/4_Symptoms/modules_np_headache_nocibersort.xlsx")
fgsea_np_headache_nocibersort$Variable <- "Headache"
fgsea_np_headache_nocibersort <- fgsea_np_headache_nocibersort[,vars_to_keep]
fgsea_np_abd_pain_nocibersort <- read_xlsx("Statistical_Analyses/4_Symptoms/modules_np_abd_pain_nocibersort.xlsx")
fgsea_np_abd_pain_nocibersort$Variable <- "Abdominal pain"
fgsea_np_abd_pain_nocibersort <- fgsea_np_abd_pain_nocibersort[,vars_to_keep]
fgsea_np_anosmia_nocibersort <- read_xlsx("Statistical_Analyses/4_Symptoms/modules_np_anosmia_nocibersort.xlsx")
fgsea_np_anosmia_nocibersort$Variable <- "Loss of smell"
fgsea_np_anosmia_nocibersort <- fgsea_np_anosmia_nocibersort[,vars_to_keep]
fgsea_np_dysgeusia_nocibersort <- read_xlsx("Statistical_Analyses/4_Symptoms/modules_np_dysgeusia_nocibersort.xlsx")
fgsea_np_dysgeusia_nocibersort$Variable <- "Loss of taste"
fgsea_np_dysgeusia_nocibersort <- fgsea_np_dysgeusia_nocibersort[,vars_to_keep]
fgsea_np_myalgias_nocibersort <- read_xlsx("Statistical_Analyses/4_Symptoms/modules_np_myalgias_nocibersort.xlsx")
fgsea_np_myalgias_nocibersort$Variable <- "Myalgias"
fgsea_np_myalgias_nocibersort <- fgsea_np_myalgias_nocibersort[,vars_to_keep]
fgsea_np <- rbind(fgsea_np_fever_nocibersort, fgsea_np_cough_nocibersort, fgsea_np_rhinorrhea_nocibersort, fgsea_np_congestion_nocibersort, fgsea_np_headache_nocibersort, fgsea_np_abd_pain_nocibersort, 
                  fgsea_np_anosmia_nocibersort, fgsea_np_dysgeusia_nocibersort, fgsea_np_myalgias_nocibersort)
fgsea_np$NES_if_sig <- fgsea_np$NES
fgsea_np$NES_if_sig[fgsea_np$padj>=0.05] <- NA
fgsea_np$Variable <- factor(fgsea_np$Variable, levels=c("Loss of taste","Loss of smell","Rhinorrhea","Congestion","Myalgias","Abdominal pain","Headache","Cough","Fever"))

fgsea_np_innate <- subset(fgsea_np, pathway %in% innate_pathways)
fgsea_np_innate$pathway <- str_to_sentence(fgsea_np_innate$pathway)
fgsea_np_innate$pathway[fgsea_np_innate$pathway=="Type i interferon signaling"] <- "Type I interferon signaling"
fgsea_np_innate$pathway[fgsea_np_innate$pathway=="Type ii interferon signaling"] <- "Type II interferon signaling"
fgsea_np_innate$pathway[fgsea_np_innate$pathway=="Type iii interferon signaling"] <- "Type III interferon signaling"
fgsea_np_innate$pathway[fgsea_np_innate$pathway=="Tnf signaling"] <- "TNF signaling"
fgsea_np_innate$pathway[fgsea_np_innate$pathway=="Tlr signaling"] <- "Toll-like receptor signaling"
fgsea_np_innate$pathway[fgsea_np_innate$pathway=="Nlr signaling"] <- "Nod-like receptor signaling"
fgsea_np_innate$pathway[fgsea_np_innate$pathway=="Tlr signaling"] <- "Toll-like receptor signaling"
fgsea_np_innate$pathway[fgsea_np_innate$pathway=="Rna sensing"] <- "RNA sensing"
fgsea_np_innate$pathway[fgsea_np_innate$pathway=="Nk activity"] <- "NK cell activity"
fgsea_np_innate$pathway[fgsea_np_innate$pathway=="Leukotriene and prostaglandin inflammation"] <- "Leukotriene/prostaglandin inflammation"
fgsea_np_innate$pathway <- factor(fgsea_np_innate$pathway, levels=c("Innate immune cell activation","Interferon response","Type I interferon signaling","Type II interferon signaling",
                                                                      "Type III interferon signaling","TNF signaling","Toll-like receptor signaling","Nod-like receptor signaling",
                                                                      "Chemokine signaling","RNA sensing","Phagocytosis","Myeloid inflammation","Myeloid activation","Inflammasomes", "NK cell activity",
                                                                      "Complement system","Coagulation","Leukotriene/prostaglandin inflammation"))
heatmap_np_innate <- fgsea_np_innate %>% ggplot(aes(x = pathway, y = Variable, fill = NES, label=round(NES_if_sig,2))) + geom_tile() + 
  labs(x = "Innate immunity", y = "Upper respiratory", fill = "NES") + scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-4,4)) + geom_text() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", fill = NA, linewidth=0.5), 
        axis.title.y = element_text(size=12, face="bold", margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=10), axis.text.x = element_blank()) 

fgsea_np_adaptive <- subset(fgsea_np, pathway %in% adaptive_pathways)
fgsea_np_adaptive$pathway <- str_to_sentence(fgsea_np_adaptive$pathway)
fgsea_np_adaptive$pathway[fgsea_np_adaptive$pathway=="Mhc class i antigen presentation"] <- "MHC class I presentation"
fgsea_np_adaptive$pathway[fgsea_np_adaptive$pathway=="Mhc class ii antigen presentation"] <- "MHC class II presentation"
fgsea_np_adaptive$pathway[fgsea_np_adaptive$pathway=="Bcr signaling"] <- "BCR signaling"
fgsea_np_adaptive$pathway[fgsea_np_adaptive$pathway=="Nf-kappab signaling"] <- "NF-kappaB signaling"
fgsea_np_adaptive$pathway[fgsea_np_adaptive$pathway=="Tcr signaling"] <- "TCR signaling"
fgsea_np_adaptive$pathway[fgsea_np_adaptive$pathway=="Jak-stat signaling"] <- "JAK-STAT signaling"
fgsea_np_adaptive$pathway <- factor(fgsea_np_adaptive$pathway, levels=c("Adaptive immune response","MHC class I presentation","MHC class II presentation","TCR signaling","BCR signaling",
                                                                        "NF-kappaB signaling","JAK-STAT signaling","Mononuclear cell migration","Lymphocyte trafficking",
                                                                        "Immune memory","Immune exhaustion","Treg differentiation"))
heatmap_np_adaptive <- fgsea_np_adaptive %>% ggplot(aes(x = pathway, y = Variable, fill = NES, label=round(NES_if_sig,2))) + geom_tile() + 
  labs(x = "Adaptive immunity", y = "Nasopharyngeal", fill = "NES") + scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-4,4)) + geom_text() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", fill = NA, linewidth=0.5), 
        axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank()) 

heatmap_np <- plot_grid(heatmap_np_innate, heatmap_np_adaptive, labels=NULL, nrow=1, rel_widths = c(0.59,0.36), align="h") 

# PAX Heatmap

# Load the required datasets
fgsea_pax_fever_cibersort <- read_xlsx("Statistical_Analyses/4_Symptoms/modules_pax_fever_cibersort.xlsx")
fgsea_pax_fever_cibersort$Variable <- "Fever"
fgsea_pax_fever_cibersort <- fgsea_pax_fever_cibersort[,vars_to_keep]
fgsea_pax_cough_cibersort <- read_xlsx("Statistical_Analyses/4_Symptoms/modules_pax_cough_cibersort.xlsx")
fgsea_pax_cough_cibersort$Variable <- "Cough"
fgsea_pax_cough_cibersort <- fgsea_pax_cough_cibersort[,vars_to_keep]
fgsea_pax_rhinorrhea_cibersort <- read_xlsx("Statistical_Analyses/4_Symptoms/modules_pax_rhinorrhea_cibersort.xlsx")
fgsea_pax_rhinorrhea_cibersort$Variable <- "Rhinorrhea"
fgsea_pax_rhinorrhea_cibersort <- fgsea_pax_rhinorrhea_cibersort[,vars_to_keep]
fgsea_pax_congestion_cibersort <- read_xlsx("Statistical_Analyses/4_Symptoms/modules_pax_congestion_cibersort.xlsx")
fgsea_pax_congestion_cibersort$Variable <- "Congestion"
fgsea_pax_congestion_cibersort <- fgsea_pax_congestion_cibersort[,vars_to_keep]
fgsea_pax_headache_cibersort <- read_xlsx("Statistical_Analyses/4_Symptoms/modules_pax_headache_cibersort.xlsx")
fgsea_pax_headache_cibersort$Variable <- "Headache"
fgsea_pax_headache_cibersort <- fgsea_pax_headache_cibersort[,vars_to_keep]
fgsea_pax_abd_pain_cibersort <- read_xlsx("Statistical_Analyses/4_Symptoms/modules_pax_abd_pain_cibersort.xlsx")
fgsea_pax_abd_pain_cibersort$Variable <- "Abdominal pain"
fgsea_pax_abd_pain_cibersort <- fgsea_pax_abd_pain_cibersort[,vars_to_keep]
fgsea_pax_anosmia_cibersort <- read_xlsx("Statistical_Analyses/4_Symptoms/modules_pax_anosmia_cibersort.xlsx")
fgsea_pax_anosmia_cibersort$Variable <- "Loss of smell"
fgsea_pax_anosmia_cibersort <- fgsea_pax_anosmia_cibersort[,vars_to_keep]
fgsea_pax_dysgeusia_cibersort <- read_xlsx("Statistical_Analyses/4_Symptoms/modules_pax_dysgeusia_cibersort.xlsx")
fgsea_pax_dysgeusia_cibersort$Variable <- "Loss of taste"
fgsea_pax_dysgeusia_cibersort <- fgsea_pax_dysgeusia_cibersort[,vars_to_keep]
fgsea_pax_myalgias_cibersort <- read_xlsx("Statistical_Analyses/4_Symptoms/modules_pax_myalgias_cibersort.xlsx")
fgsea_pax_myalgias_cibersort$Variable <- "Myalgias"
fgsea_pax_myalgias_cibersort <- fgsea_pax_myalgias_cibersort[,vars_to_keep]
fgsea_pax <- rbind(fgsea_pax_fever_cibersort, fgsea_pax_cough_cibersort, fgsea_pax_rhinorrhea_cibersort, fgsea_pax_congestion_cibersort, fgsea_pax_headache_cibersort, fgsea_pax_abd_pain_cibersort, 
                  fgsea_pax_anosmia_cibersort, fgsea_pax_dysgeusia_cibersort, fgsea_pax_myalgias_cibersort)
fgsea_pax$NES_if_sig <- fgsea_pax$NES
fgsea_pax$NES_if_sig[fgsea_pax$padj>=0.05] <- NA
fgsea_pax$Variable <- factor(fgsea_pax$Variable, levels=c("Loss of taste","Loss of smell","Rhinorrhea","Congestion","Myalgias","Abdominal pain","Headache","Cough","Fever"))

fgsea_pax_innate <- subset(fgsea_pax, pathway %in% innate_pathways)
fgsea_pax_innate$pathway <- str_to_sentence(fgsea_pax_innate$pathway)
fgsea_pax_innate$pathway[fgsea_pax_innate$pathway=="Type i interferon signaling"] <- "Type I interferon signaling"
fgsea_pax_innate$pathway[fgsea_pax_innate$pathway=="Type ii interferon signaling"] <- "Type II interferon signaling"
fgsea_pax_innate$pathway[fgsea_pax_innate$pathway=="Type iii interferon signaling"] <- "Type III interferon signaling"
fgsea_pax_innate$pathway[fgsea_pax_innate$pathway=="Tnf signaling"] <- "TNF signaling"
fgsea_pax_innate$pathway[fgsea_pax_innate$pathway=="Tlr signaling"] <- "Toll-like receptor signaling"
fgsea_pax_innate$pathway[fgsea_pax_innate$pathway=="Nlr signaling"] <- "Nod-like receptor signaling"
fgsea_pax_innate$pathway[fgsea_pax_innate$pathway=="Tlr signaling"] <- "Toll-like receptor signaling"
fgsea_pax_innate$pathway[fgsea_pax_innate$pathway=="Rna sensing"] <- "RNA sensing"
fgsea_pax_innate$pathway[fgsea_pax_innate$pathway=="Nk activity"] <- "NK cell activity"
fgsea_pax_innate$pathway[fgsea_pax_innate$pathway=="Leukotriene and prostaglandin inflammation"] <- "Leukotriene/prostaglandin inflammation"
fgsea_pax_innate$pathway <- factor(fgsea_pax_innate$pathway, levels=c("Innate immune cell activation","Interferon response","Type I interferon signaling","Type II interferon signaling",
                                                                    "Type III interferon signaling","TNF signaling","Toll-like receptor signaling","Nod-like receptor signaling",
                                                                    "Chemokine signaling","RNA sensing","Phagocytosis","Myeloid inflammation","Myeloid activation","Inflammasomes", "NK cell activity",
                                                                    "Complement system","Coagulation","Leukotriene/prostaglandin inflammation"))
heatmap_pax_innate <- fgsea_pax_innate %>% ggplot(aes(x = pathway, y = Variable, fill = NES, label=round(NES_if_sig,2))) + geom_tile() + 
  labs(x = "Innate immunity", y = "Peripheral blood", fill = "NES") + scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-4,4)) + geom_text() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", fill = NA, linewidth=0.5), 
        axis.title.y = element_text(size=12, face="bold", margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=12, face="bold", margin = ggplot2::margin(t = 0, r = 0, b = 5, l = 0)),
        axis.text.y = element_text(size=10), axis.text.x = element_text(size=10, color = "black", angle=45, hjust=1))

fgsea_pax_adaptive <- subset(fgsea_pax, pathway %in% adaptive_pathways)
fgsea_pax_adaptive$pathway <- str_to_sentence(fgsea_pax_adaptive$pathway)
fgsea_pax_adaptive$pathway[fgsea_pax_adaptive$pathway=="Mhc class i antigen presentation"] <- "MHC class I presentation"
fgsea_pax_adaptive$pathway[fgsea_pax_adaptive$pathway=="Mhc class ii antigen presentation"] <- "MHC class II presentation"
fgsea_pax_adaptive$pathway[fgsea_pax_adaptive$pathway=="Bcr signaling"] <- "B cell receptor signaling"
fgsea_pax_adaptive$pathway[fgsea_pax_adaptive$pathway=="Nf-kappab signaling"] <- "NF-kappaB signaling"
fgsea_pax_adaptive$pathway[fgsea_pax_adaptive$pathway=="Tcr signaling"] <- "T cell receptor signaling"
fgsea_pax_adaptive$pathway[fgsea_pax_adaptive$pathway=="Jak-stat signaling"] <- "JAK-STAT signaling"
fgsea_pax_adaptive$pathway <- factor(fgsea_pax_adaptive$pathway, levels=c("Adaptive immune response","MHC class I presentation","MHC class II presentation","T cell receptor signaling",
                                                                          "B cell receptor signaling", "NF-kappaB signaling","JAK-STAT signaling","Mononuclear cell migration",
                                                                          "Lymphocyte trafficking", "Immune memory","Immune exhaustion","Treg differentiation"))
heatmap_pax_adaptive <- fgsea_pax_adaptive %>% ggplot(aes(x = pathway, y = Variable, fill = NES, label=round(NES_if_sig,2))) + geom_tile() + 
  labs(x = "Adaptive immunity", y = "Peripheral blood", fill = "NES") + scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-4,4)) + geom_text() +
  theme_classic() + 
  theme(legend.position = "none",
        panel.background = element_rect(colour = "black", fill = NA, linewidth=0.5), 
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=12, face="bold", margin = ggplot2::margin(t = 0, r = 0, b = 5, l = 0)),
        axis.text.y = element_blank(), axis.text.x = element_text(size=10, color = "black", angle=45, hjust=1))

heatmap_pax <- plot_grid(heatmap_pax_innate, heatmap_pax_adaptive, labels=NULL, nrow=1, rel_widths = c(0.59,0.36), align="h") 

heatmaps <- plot_grid(heatmap_np, heatmap_pax, labels=NULL, ncol=1, rel_heights = c(0.37,0.628), align="v") 

legend_plot <- fgsea_pax %>% ggplot(aes(x = pathway, y = Variable, fill = NES, label=round(NES_if_sig,2))) + geom_tile() + 
  labs(x = NULL, y = NULL, fill = "NES") + scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-4,4)) + geom_text() +
  theme_classic() + 
  theme(legend.position=c(0.25,0.65), legend.text = element_text(size=12))
legend <- get_legend(legend_plot)

png(file="Statistical_Analyses/Figures/Figure_4.png", width = 18, height = 8, units = 'in', res = 1200)
plot_grid(heatmaps, legend, labels=NULL, ncol=2, rel_widths = c(1,0.1)) 
dev.off()