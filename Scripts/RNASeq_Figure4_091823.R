# BRAVE Kids RNA Sequencing Analysis
# Aditya Mohan (MD/PhD candidate)  / Matthew Kelly, MD, MPH 
# Figure 4 
# Last update: Sept. 18, 2023

remove(list=ls())
setwd("G:/My Drive/Research/BRAVE_Kids/RNA_Sequencing/Statistical_Analyses") 
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
fgsea_np_fever <- read_xlsx("4_Symptoms/modules_np_fever.xlsx")
fgsea_np_fever$Variable <- "Fever"
fgsea_np_fever <- fgsea_np_fever[,vars_to_keep]
fgsea_np_cough <- read_xlsx("4_Symptoms/modules_np_cough.xlsx")
fgsea_np_cough$Variable <- "Cough"
fgsea_np_cough <- fgsea_np_cough[,vars_to_keep]
fgsea_np_rhinorrhea <- read_xlsx("4_Symptoms/modules_np_rhinorrhea.xlsx")
fgsea_np_rhinorrhea$Variable <- "Rhinorrhea"
fgsea_np_rhinorrhea <- fgsea_np_rhinorrhea[,vars_to_keep]
fgsea_np_congestion <- read_xlsx("4_Symptoms/modules_np_congestion.xlsx")
fgsea_np_congestion$Variable <- "Congestion"
fgsea_np_congestion <- fgsea_np_congestion[,vars_to_keep]
fgsea_np_headache <- read_xlsx("4_Symptoms/modules_np_headache.xlsx")
fgsea_np_headache$Variable <- "Headache"
fgsea_np_headache <- fgsea_np_headache[,vars_to_keep]
fgsea_np_abd_pain <- read_xlsx("4_Symptoms/modules_np_abd_pain.xlsx")
fgsea_np_abd_pain$Variable <- "Abdominal pain"
fgsea_np_abd_pain <- fgsea_np_abd_pain[,vars_to_keep]
fgsea_np_anosmia <- read_xlsx("4_Symptoms/modules_np_anosmia.xlsx")
fgsea_np_anosmia$Variable <- "Loss of smell"
fgsea_np_anosmia <- fgsea_np_anosmia[,vars_to_keep]
fgsea_np_dysgeusia <- read_xlsx("4_Symptoms/modules_np_dysgeusia.xlsx")
fgsea_np_dysgeusia$Variable <- "Loss of taste"
fgsea_np_dysgeusia <- fgsea_np_dysgeusia[,vars_to_keep]
fgsea_np_myalgias <- read_xlsx("4_Symptoms/modules_np_myalgias.xlsx")
fgsea_np_myalgias$Variable <- "Myalgias"
fgsea_np_myalgias <- fgsea_np_myalgias[,vars_to_keep]
fgsea_np <- rbind(fgsea_np_fever, fgsea_np_cough, fgsea_np_rhinorrhea, fgsea_np_congestion, fgsea_np_headache, fgsea_np_abd_pain, 
                  fgsea_np_anosmia, fgsea_np_dysgeusia, fgsea_np_myalgias)
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
fgsea_pax_fever <- read_xlsx("4_Symptoms/modules_pax_fever.xlsx")
fgsea_pax_fever$Variable <- "Fever"
fgsea_pax_fever <- fgsea_pax_fever[,vars_to_keep]
fgsea_pax_cough <- read_xlsx("4_Symptoms/modules_pax_cough.xlsx")
fgsea_pax_cough$Variable <- "Cough"
fgsea_pax_cough <- fgsea_pax_cough[,vars_to_keep]
fgsea_pax_rhinorrhea <- read_xlsx("4_Symptoms/modules_pax_rhinorrhea.xlsx")
fgsea_pax_rhinorrhea$Variable <- "Rhinorrhea"
fgsea_pax_rhinorrhea <- fgsea_pax_rhinorrhea[,vars_to_keep]
fgsea_pax_congestion <- read_xlsx("4_Symptoms/modules_pax_congestion.xlsx")
fgsea_pax_congestion$Variable <- "Congestion"
fgsea_pax_congestion <- fgsea_pax_congestion[,vars_to_keep]
fgsea_pax_headache <- read_xlsx("4_Symptoms/modules_pax_headache.xlsx")
fgsea_pax_headache$Variable <- "Headache"
fgsea_pax_headache <- fgsea_pax_headache[,vars_to_keep]
fgsea_pax_abd_pain <- read_xlsx("4_Symptoms/modules_pax_abd_pain.xlsx")
fgsea_pax_abd_pain$Variable <- "Abdominal pain"
fgsea_pax_abd_pain <- fgsea_pax_abd_pain[,vars_to_keep]
fgsea_pax_anosmia <- read_xlsx("4_Symptoms/modules_pax_anosmia.xlsx")
fgsea_pax_anosmia$Variable <- "Loss of smell"
fgsea_pax_anosmia <- fgsea_pax_anosmia[,vars_to_keep]
fgsea_pax_dysgeusia <- read_xlsx("4_Symptoms/modules_pax_dysgeusia.xlsx")
fgsea_pax_dysgeusia$Variable <- "Loss of taste"
fgsea_pax_dysgeusia <- fgsea_pax_dysgeusia[,vars_to_keep]
fgsea_pax_myalgias <- read_xlsx("4_Symptoms/modules_pax_myalgias.xlsx")
fgsea_pax_myalgias$Variable <- "Myalgias"
fgsea_pax_myalgias <- fgsea_pax_myalgias[,vars_to_keep]
fgsea_pax <- rbind(fgsea_pax_fever, fgsea_pax_cough, fgsea_pax_rhinorrhea, fgsea_pax_congestion, fgsea_pax_headache, fgsea_pax_abd_pain, 
                  fgsea_pax_anosmia, fgsea_pax_dysgeusia, fgsea_pax_myalgias)
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
fgsea_pax_adaptive$pathway[fgsea_pax_adaptive$pathway=="Bcr signaling"] <- "BCR signaling"
fgsea_pax_adaptive$pathway[fgsea_pax_adaptive$pathway=="Nf-kappab signaling"] <- "NF-kappaB signaling"
fgsea_pax_adaptive$pathway[fgsea_pax_adaptive$pathway=="Tcr signaling"] <- "TCR signaling"
fgsea_pax_adaptive$pathway[fgsea_pax_adaptive$pathway=="Jak-stat signaling"] <- "JAK-STAT signaling"
fgsea_pax_adaptive$pathway <- factor(fgsea_pax_adaptive$pathway, levels=c("Adaptive immune response","MHC class I presentation","MHC class II presentation","TCR signaling","BCR signaling",
                                                                        "NF-kappaB signaling","JAK-STAT signaling","Mononuclear cell migration","Lymphocyte trafficking",
                                                                        "Immune memory","Immune exhaustion","Treg differentiation"))
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

# NO TITLE

#png(file="Figures/Figure_4.png", width = 18, height = 8, units = 'in', res = 1200)
#plot_grid(heatmaps, legend, labels=NULL, ncol=2, rel_widths = c(1,0.1)) 
#dev.off()

# ADDITION OF TITLE

fig4_plot <- plot_grid(heatmaps, legend, labels=NULL, ncol=2, rel_widths = c(1,0.1)) 
fig4_title <- ggdraw() + draw_label("Figure 4. Differential expression of immune modules in SARS-CoV-2 infection by symptom presence", 
                                    size=13, fontface='bold')

png(file="Figures/Figure_4.png", width = 18, height = 8.4, units = 'in', res = 1200)
plot_grid(fig4_title, NULL, fig4_plot, labels=NULL, ncol=1, rel_heights=c(0.03,0.02,0.97), align="v") 
dev.off()