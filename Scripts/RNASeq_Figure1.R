# BRAVE Kids RNA Sequencing Analysis
# Aditya Mohan (MD/PhD candidate)  / Matthew Kelly, MD, MPH 
# Figure 1 - comparison of cell populations and gene module expression in upper respiratory & peripheral blood samples in healthy individuals by age
# Analyses of gene module expression in peripheral blood samples adjusted for imputed cell proportions
# Last update: June 8, 2024

remove(list=ls())
setwd("______________________________") 
set.seed(1234)
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
library(patchwork)
library(rstatix)
library(data.table)
library(stringr)
library(RColorBrewer)
library(phyloseq)

# Specify expressions for figures
label_log2FC <- expression(bold(log[2] ~ fold ~ change))
label_adjp <- expression(bold(-log[10] ~ p[adj]))

# Upload files with CIBERSORTx data
phy.rnaseq.np <- readRDS("phy.rnaseq.np.rds")
phy.rnaseq.np.neg <- subset_samples(phy.rnaseq.np, corona=="Negative")
metadata_np_neg <- data.frame(sample_data(phy.rnaseq.np.neg))
metadata_np_neg_pops <- metadata_np_neg[,c("age_cat","B_cells","Plasma_cells", "T_cells_CD8","T_cells_CD4","T_cells_gamma_delta","NK_cells",
                                             "Mono_Macrophages","Dendritic_cells","Mast_cells","Eosinophils","Neutrophils")]
metadata_np_neg_pops$age_cat <- factor(metadata_np_neg_pops$age_cat, levels=c("0-5 years", "6-13 years", "14-20 years"))
metadata_np_neg_pops <- gather(metadata_np_neg_pops, cell_type, value, B_cells:Neutrophils, factor_key=TRUE)
metadata_np_neg_pops$cell_type2[metadata_np_neg_pops$cell_type=="B_cells"] <- "B cells"
metadata_np_neg_pops$cell_type2[metadata_np_neg_pops$cell_type=="Plasma_cells"] <- "Plasma cells"
metadata_np_neg_pops$cell_type2[metadata_np_neg_pops$cell_type=="T_cells_CD8"] <- "CD4+ T cells"
metadata_np_neg_pops$cell_type2[metadata_np_neg_pops$cell_type=="T_cells_CD4"] <- "CD8+ T cells"
metadata_np_neg_pops$cell_type2[metadata_np_neg_pops$cell_type=="T_cells_gamma_delta"] <- "Gamma delta T cells"
metadata_np_neg_pops$cell_type2[metadata_np_neg_pops$cell_type=="NK_cells"] <- "NK cells"
metadata_np_neg_pops$cell_type2[metadata_np_neg_pops$cell_type=="Mono_Macrophages"] <- "Mono/macrophages"
metadata_np_neg_pops$cell_type2[metadata_np_neg_pops$cell_type=="Dendritic_cells"] <- "Dendritic cells"
metadata_np_neg_pops$cell_type2[metadata_np_neg_pops$cell_type=="Mast_cells"] <- "Mast cells"
metadata_np_neg_pops$cell_type2[metadata_np_neg_pops$cell_type=="Eosinophils"] <- "Eosinophils"
metadata_np_neg_pops$cell_type2[metadata_np_neg_pops$cell_type=="Neutrophils"] <- "Neutrophils"
tapply(metadata_np_neg_pops$value, metadata_np_neg_pops$cell_type2, summary) 
# Exclude cell types identified in <25% of samples from figure
metadata_np_neg_pops <- subset(metadata_np_neg_pops, cell_type2!="Gamma delta T cells")

phy.rnaseq.pax <- readRDS("phy.rnaseq.pax.rds")
phy.rnaseq.pax.neg <- subset_samples(phy.rnaseq.pax, corona=="Negative")
metadata_pax_neg <- data.frame(sample_data(phy.rnaseq.pax.neg))
metadata_pax_neg_pops <- metadata_pax_neg[,c("age_cat","B_cells","Plasma_cells", "T_cells_CD8","T_cells_CD4","T_cells_gamma_delta","NK_cells",
                                           "Mono_Macrophages","Dendritic_cells","Mast_cells","Eosinophils","Neutrophils")]
metadata_pax_neg_pops$age_cat <- factor(metadata_pax_neg_pops$age_cat, levels=c("0-5 years", "6-13 years", "14-20 years", "Adult"))
metadata_pax_neg_pops <- gather(metadata_pax_neg_pops, cell_type, value, B_cells:Neutrophils, factor_key=TRUE)
metadata_pax_neg_pops$cell_type2[metadata_pax_neg_pops$cell_type=="B_cells"] <- "B cells"
metadata_pax_neg_pops$cell_type2[metadata_pax_neg_pops$cell_type=="Plasma_cells"] <- "Plasma cells"
metadata_pax_neg_pops$cell_type2[metadata_pax_neg_pops$cell_type=="T_cells_CD8"] <- "CD4+ T cells"
metadata_pax_neg_pops$cell_type2[metadata_pax_neg_pops$cell_type=="T_cells_CD4"] <- "CD8+ T cells"
metadata_pax_neg_pops$cell_type2[metadata_pax_neg_pops$cell_type=="T_cells_gamma_delta"] <- "Gamma delta T cells"
metadata_pax_neg_pops$cell_type2[metadata_pax_neg_pops$cell_type=="NK_cells"] <- "NK cells"
metadata_pax_neg_pops$cell_type2[metadata_pax_neg_pops$cell_type=="Mono_Macrophages"] <- "Mono/macrophages"
metadata_pax_neg_pops$cell_type2[metadata_pax_neg_pops$cell_type=="Dendritic_cells"] <- "Dendritic cells"
metadata_pax_neg_pops$cell_type2[metadata_pax_neg_pops$cell_type=="Mast_cells"] <- "Mast cells"
metadata_pax_neg_pops$cell_type2[metadata_pax_neg_pops$cell_type=="Eosinophils"] <- "Eosinophils"
metadata_pax_neg_pops$cell_type2[metadata_pax_neg_pops$cell_type=="Neutrophils"] <- "Neutrophils"
tapply(metadata_pax_neg_pops$value, metadata_pax_neg_pops$cell_type2, summary) 
# Exclude cell types identified in <25% of samples from figure
metadata_pax_neg_pops <- subset(metadata_pax_neg_pops, cell_type2!="Dendritic cells" & cell_type2!="Eosinophils" & cell_type2!="Gamma delta T cells")

# Upload FGSEA output files
fgsea_np_neg_0to5_vs_6to13_nocibersort <- data.frame(read_excel("Statistical_Analyses/1_COVID_Neg_by_Age/modules_np_neg_0to5_vs_6to13_nocibersort.xlsx"))
fgsea_np_neg_0to5_vs_14to20_nocibersort <- data.frame(read_excel("Statistical_Analyses/1_COVID_Neg_by_Age/modules_np_neg_0to5_vs_14to20_nocibersort.xlsx"))
fgsea_np_neg_6to13_vs_14to20_nocibersort <- data.frame(read_excel("Statistical_Analyses/1_COVID_Neg_by_Age/modules_np_neg_6to13_vs_14to20_nocibersort.xlsx"))
fgsea_pax_neg_0to5_vs_6to13_cibersort <- data.frame(read_excel("Statistical_Analyses/1_COVID_Neg_by_Age/modules_pax_neg_0to5_vs_6to13_cibersort.xlsx"))
fgsea_pax_neg_0to5_vs_14to20_cibersort <- data.frame(read_excel("Statistical_Analyses/1_COVID_Neg_by_Age/modules_pax_neg_0to5_vs_14to20_cibersort.xlsx"))
fgsea_pax_neg_6to13_vs_14to20_cibersort <- data.frame(read_excel("Statistical_Analyses/1_COVID_Neg_by_Age/modules_pax_neg_6to13_vs_14to20_cibersort.xlsx"))
fgsea_pax_neg_0to5_vs_adult_cibersort <- data.frame(read_excel("Statistical_Analyses/1_COVID_Neg_by_Age/modules_pax_neg_0to5_vs_adult_cibersort.xlsx"))
fgsea_pax_neg_6to13_vs_adult_cibersort <- data.frame(read_excel("Statistical_Analyses/1_COVID_Neg_by_Age/modules_pax_neg_6to13_vs_adult_cibersort.xlsx"))
fgsea_pax_neg_14to20_vs_adult_cibersort <- data.frame(read_excel("Statistical_Analyses/1_COVID_Neg_by_Age/modules_pax_neg_14to20_vs_adult_cibersort.xlsx"))

# CIBERSORTx cell populations barplot

# *, p<0.05
# **, p<0.01
# ***, p<0.001
# ****, p<0.0001

# UPPER RESPIRATORY SAMPLES (COVID-NEGATIVE BY AGE)
# Adjust comparisons for multiple testing using BH (n=10)
# B cells: 0.74
# Plasma cells: >0.99
# CD8+ T cells: >0.99
# CD4+ T cells: >0.99
# NK cells: >0.99
# Mono/macrophages: >0.99
# Dendritic cells: 0.27
# Mast cells: >0.99
# Eosinophils: >0.99
# Neutrophils: >0.99

cibersort_np_neg <- ggplot(metadata_np_neg_pops, aes(x=cell_type2, y=value, fill=age_cat)) + geom_boxplot() + ylab("Imputed immune cell proportions") + ylim(0,0.51) +
  theme(legend.title = element_blank(), legend.position = "right", legend.text = element_text(size=12), legend.box.spacing = unit(3, "pt"), legend.key=element_rect(fill="white"),
        panel.background = element_rect(colour = "black", fill = NA, linewidth=0.5), panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey70"),
        panel.grid.minor.y = element_line(linewidth = 0.2, colour = "grey70"),
        axis.title.y = element_text(size=13, face="bold", margin = ggplot2::margin(t = 0, r = 5, b = 0, l = 0)),
        axis.title.x = element_blank(), axis.text.y = element_text(size=11.5), axis.text.x = element_text(size=11.5, color = "black", angle=35, hjust=1), 
        plot.title = element_text(size=14, hjust = 0.5, face="bold")) + scale_fill_brewer(palette="Paired") +
  ggtitle("Upper respiratory") # No significant differences in cell populations

# PERIPHERAL BLOOD SAMPLES (COVID-NEGATIVE BY AGE)
# Adjust comparisons for multiple testing using BH (n=8)
# B cells: 9.03E-14
# Plasma cells: 0.005
# CD8+ T cells: 0.002
# CD4+ T cells: >0.99
# NK cells: 0.69
# Mono/macrophages: 1.95E-5
# Mast cells: 0.15
# Neutrophils: 0.00094

cibersort_pax_neg <- ggplot(metadata_pax_neg_pops, aes(x=cell_type2, y=value, fill=age_cat)) + geom_boxplot() + ylab("Imputed immune cell proportions") + ylim(0,0.41) +
  theme(legend.title = element_blank(), legend.position = "right", legend.text = element_text(size=12), legend.box.spacing = unit(3, "pt"), legend.key=element_rect(fill="white"),
        panel.background = element_rect(colour = "black", fill = NA, linewidth=0.5), panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey70"),
        panel.grid.minor.y = element_line(linewidth = 0.2, colour = "grey70"),
        axis.title.y = element_text(size=13, face="bold", margin = ggplot2::margin(t = 0, r = 5, b = 0, l = 0)),
        axis.title.x = element_blank(), axis.text.y = element_text(size=11.5), axis.text.x = element_text(size=11.5, color = "black", angle=35, hjust=1), 
        plot.title = element_text(size=14, hjust = 0.5, face="bold")) + scale_fill_brewer(palette="Paired") +
  ggtitle("Peripheral blood") + 
  annotate("text", x=1, y=0.41, size=6, label= "****") + # B cells (beta regression, p=9.03E-14) 
  annotate("text", x=8, y=0.41, size=6, label= "**") + # Plasma cells (beta regression, 0.005)
  annotate("text", x=3, y=0.41, size=6, label= "**") + # CD8 T cells (beta regression, p=0.002) 
  annotate("text", x=5, y=0.41, size=6, label= "****") + # Mono/macrophages (beta regression, p=1.95E-5)
  annotate("text", x=6, y=0.41, size=6, label= "***") # Neutrophils (beta regression, p=0.00094) 

fig1_cibersort <- plot_grid(NULL, cibersort_np_neg, NULL, cibersort_pax_neg, ncol=1, labels=c("","a","","b"), rel_heights=c(0.05,3,0.05,3), align="v") 
fig1_cibersort

# FGSEA modules dotplot - PEDIATRIC AGE GROUPS

modules_np_peds <- rbindlist(list(fgsea_np_neg_0to5_vs_6to13_nocibersort, fgsea_np_neg_0to5_vs_14to20_nocibersort, fgsea_np_neg_6to13_vs_14to20_nocibersort))
modules_np_peds$pathway <- str_to_sentence(modules_np_peds$pathway)
modules_np_peds$pathway[modules_np_peds$pathway=="Bcr signaling"] <- "B cell receptor signaling"
modules_np_peds$pathway[modules_np_peds$pathway=="Nk activity"] <- "NK cell activity"
modules_np_peds$pathway[modules_np_peds$pathway=="Tcr signaling"] <- "T cell receptor signaling"
modules_np_peds$pathway[modules_np_peds$pathway=="Type i interferon signaling"] <- "Type I interferon signaling"
modules_np_peds$pathway[modules_np_peds$pathway=="Type ii interferon signaling"] <- "Type II interferon signaling"
modules_np_peds$pathway[modules_np_peds$pathway=="Innate immune cell activation"] <- "Innate immune activation"
modules_np_peds$pathway[modules_np_peds$pathway=="Tnf signaling"] <- "TNF signaling"
modules_np_peds$pathway[modules_np_peds$pathway=="Tlr signaling"] <- "Toll-like receptor signaling"
modules_np_peds$pathway[modules_np_peds$pathway=="Treg differentiation"] <- "Treg cell differentiation"
modules_np_peds$Sample[modules_np_peds$Sample=="6-13 vs. 14-20 yr"] <- " 6-13 vs. 14-20 yr"
modules_np_peds$Sample <- factor(modules_np_peds$Sample, levels=c("0-5 vs. 6-13 yr","0-5 vs. 14-20 yr", " 6-13 vs. 14-20 yr"))

if (is.element('Up-regulated', modules_np_peds$Expression) & is.element('Unchanged', modules_np_peds$Expression) & is.element('Down-regulated', modules_np_peds$Expression)) {
  cen = c("lightsteelblue2", "gray80", "#F9AAAE")
} else if ( is.element('Up-regulated', modules_np_peds$Expression) & is.element('Unchanged', modules_np_peds$Expression)) {
  cen = c("gray80", "#F9AAAE")
} else if ( is.element('Unchanged', modules_np_peds$Expression) & is.element('Down-regulated', modules_np_peds$Expression)) {
  cen = c("lightsteelblue2", "gray80")
} else if ( is.element('Up-regulated', modules_np_peds$Expression) & is.element('Down-regulated', modules_np_peds$Expression)) {
  cen = c("lightsteelblue2", "#F9AAAE")
} else {cen = c("gray80")}

innate <- c("Coagulation", "Complement system", "NK cell activity", "Inflammasomes", "Innate immune activation", 
            "Myeloid inflammation", "Myeloid activation", "Phagocytosis", "Toll-like receptor signaling", "TNF signaling",
            "Type II interferon signaling", "Type I interferon signaling", "Interferon response")
modules_np_peds_innate <- subset(modules_np_peds, pathway %in% innate)
modules_np_peds_innate$pathway <- factor(modules_np_peds_innate$pathway, levels=c("Coagulation", "Complement system", "NK cell activity", "Inflammasomes", "Innate immune activation", 
                                                                        "Myeloid inflammation", "Myeloid activation", "Phagocytosis", "Toll-like receptor signaling", "TNF signaling",
                                                                        "Type II interferon signaling", "Type I interferon signaling", "Interferon response"))

innate_plot_np <- ggplot(modules_np_peds_innate, aes(x = pathway, y = Sample)) +
  geom_point(aes(size = (abs(NES))^3, fill = Expression), shape = 21) +
  geom_text(aes(label = formatC(NES, format='f', digits = 2)), parse = FALSE) +
  scale_fill_manual(values = cen) +
  labs(x= "Innate immunity", y = "") + 
  theme(legend.key=element_blank(), legend.position = "none", panel.background = element_rect(colour = "black", fill = NA, linewidth=0.8),
        axis.title.y = element_text(size=14, face="bold", margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12, color = "black"), axis.text.x = element_blank(), axis.ticks = element_blank(),
        plot.title=element_text(size=14, color="black", hjust=0.5, face="bold")) + 
  scale_size_continuous(range = c(13,15)) + ggtitle("Upper respiratory") +
  coord_flip() 
innate_plot_np

adaptive <- c("Treg cell differentiation", "Lymphocyte trafficking", "Mononuclear cell migration", 
              "T cell receptor signaling", "B cell receptor signaling", "Adaptive immune response")
modules_np_peds_adaptive <- subset(modules_np_peds, pathway %in% adaptive)
modules_np_peds_adaptive$pathway <- factor(modules_np_peds_adaptive$pathway, levels=c("Treg cell differentiation", "Lymphocyte trafficking", "Mononuclear cell migration", 
                                                                            "T cell receptor signaling", "B cell receptor signaling", "Adaptive immune response"))

adaptive_plot_np <- ggplot(modules_np_peds_adaptive, aes(x = pathway, y = Sample)) +
  geom_point(aes(size = (abs(NES))^3, fill = Expression), shape = 21) +
  geom_text(aes(label = formatC(NES, format='f', digits = 2)), parse = FALSE) +
  scale_fill_manual(values = cen) +
  labs(x= "Adaptive immunity") + 
  theme(legend.key=element_blank(), legend.position = "none", panel.background = element_rect(colour = "black", fill = NA, linewidth=0.8),
        axis.title.y = element_text(size=14, face="bold", margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12, color = "black"), 
        axis.text.x = element_text(size=12, color = "black", angle=45, hjust=1), 
        axis.ticks.y = element_blank()) + 
  scale_size_continuous(range = c(12.5,14)) +
  coord_flip() 
adaptive_plot_np

modules_pax_peds <- rbindlist(list(fgsea_pax_neg_0to5_vs_6to13_cibersort, fgsea_pax_neg_0to5_vs_14to20_cibersort, fgsea_pax_neg_6to13_vs_14to20_cibersort))
modules_pax_peds$pathway <- str_to_sentence(modules_pax_peds$pathway)
modules_pax_peds$pathway[modules_pax_peds$pathway=="Bcr signaling"] <- "B cell receptor signaling"
modules_pax_peds$pathway[modules_pax_peds$pathway=="Nk activity"] <- "NK cell activity"
modules_pax_peds$pathway[modules_pax_peds$pathway=="Tcr signaling"] <- "T cell receptor signaling"
modules_pax_peds$pathway[modules_pax_peds$pathway=="Type i interferon signaling"] <- "Type I interferon signaling"
modules_pax_peds$pathway[modules_pax_peds$pathway=="Type ii interferon signaling"] <- "Type II interferon signaling"
modules_pax_peds$pathway[modules_pax_peds$pathway=="Innate immune cell activation"] <- "Innate immune activation"
modules_pax_peds$pathway[modules_pax_peds$pathway=="Tnf signaling"] <- "TNF signaling"
modules_pax_peds$pathway[modules_pax_peds$pathway=="Tlr signaling"] <- "Toll-like receptor signaling"
modules_pax_peds$pathway[modules_pax_peds$pathway=="Treg differentiation"] <- "Treg cell differentiation"
modules_pax_peds$Sample[modules_pax_peds$Sample=="6-13 vs. 14-20 yr"] <- " 6-13 vs. 14-20 yr"
modules_pax_peds$Sample <- factor(modules_pax_peds$Sample, levels=c("0-5 vs. 6-13 yr","0-5 vs. 14-20 yr", " 6-13 vs. 14-20 yr"))

if (is.element('Up-regulated', modules_pax_peds$Expression) & is.element('Unchanged', modules_pax_peds$Expression) & is.element('Down-regulated', modules_pax_peds$Expression)) {
  cen = c("lightsteelblue2", "gray80", "#F9AAAE")
} else if ( is.element('Up-regulated', modules_pax_peds$Expression) & is.element('Unchanged', modules_pax_peds$Expression)) {
  cen = c("gray80", "#F9AAAE")
} else if ( is.element('Unchanged', modules_pax_peds$Expression) & is.element('Down-regulated', modules_pax_peds$Expression)) {
  cen = c("lightsteelblue2", "gray80")
} else if ( is.element('Up-regulated', modules_pax_peds$Expression) & is.element('Down-regulated', modules_pax_peds$Expression)) {
  cen = c("lightsteelblue2", "#F9AAAE")
} else {cen = c("gray80")}

modules_pax_peds_innate <- subset(modules_pax_peds, pathway %in% innate)
modules_pax_peds_innate$pathway <- factor(modules_pax_peds_innate$pathway, levels=c("Coagulation", "Complement system", "NK cell activity", "Inflammasomes", "Innate immune activation", 
                                                                          "Myeloid inflammation", "Myeloid activation", "Phagocytosis", "Toll-like receptor signaling", "TNF signaling",
                                                                          "Type II interferon signaling", "Type I interferon signaling", "Interferon response"))

innate_plot_pax <- ggplot(modules_pax_peds_innate, aes(x = pathway, y = Sample)) +
  geom_point(aes(size = (abs(NES))^3, fill = Expression), shape = 21) +
  geom_text(aes(label = formatC(NES, format='f', digits = 2)), parse = FALSE) +
  scale_fill_manual(values = cen) +
  labs(x= "Innate immunity", y = "") + 
  theme(legend.key=element_blank(), legend.position = "none", panel.background = element_rect(colour = "black", fill = NA, linewidth=0.8),
        axis.title.y = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_blank(), axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        plot.title=element_text(size=14, color="black", hjust=0.5, face="bold")) + 
  scale_size_continuous(range = c(13,15)) + ggtitle("Peripheral blood") +
  coord_flip() 
innate_plot_pax

modules_pax_peds_adaptive <- subset(modules_pax_peds, pathway %in% adaptive)
modules_pax_peds_adaptive$pathway <- factor(modules_pax_peds_adaptive$pathway, levels=c("Treg cell differentiation", "Lymphocyte trafficking", "Mononuclear cell migration", 
                                                                              "T cell receptor signaling", "B cell receptor signaling", "Adaptive immune response"))

adaptive_plot_pax <- ggplot(modules_pax_peds_adaptive, aes(x = pathway, y = Sample)) +
  geom_point(aes(size = (abs(NES))^3, fill = Expression), shape = 21) +
  geom_text(aes(label = formatC(NES, format='f', digits = 2)), parse = FALSE) +
  scale_fill_manual(values = cen) +
  labs(x= "Adaptive immunity", y = "") + 
  theme(legend.key=element_blank(), legend.position = "none", panel.background = element_rect(colour = "black", fill = NA, linewidth=0.8),
        axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=12, color = "black", angle=45, hjust=1)) + 
  scale_size_continuous(range = c(12.5,14)) +
  coord_flip() 
adaptive_plot_pax

module_plot_np <- plot_grid(innate_plot_np, NULL, adaptive_plot_np, labels=NULL, ncol=1, rel_heights=c(4.05,-0.06,2.5), align="v") 
module_plot_pax <- plot_grid(innate_plot_pax, NULL, adaptive_plot_pax, labels=NULL, ncol=1, rel_heights=c(4.05,-0.06,2.5), align="v") 
fig1_modules_peds <- plot_grid(module_plot_np, NULL, module_plot_pax, labels=NULL, ncol=3, rel_widths=c(4,-0.05, 1.8), align="vh") 
title_peds <- ggdraw() + draw_label("                                            Pediatric age groups", size=14, fontface='bold')
fig1_modules_peds_title <- plot_grid(title_peds, fig1_modules_peds, labels=NULL, ncol=1, rel_heights=c(0.03,0.97), align="v")

# FGSEA modules dotplot - COMPARISONS TO ADULTS (PERIPHERAL BLOOD ONLY)

modules_pax_peds_vs_adults <- rbindlist(list(fgsea_pax_neg_0to5_vs_adult_cibersort, fgsea_pax_neg_6to13_vs_adult_cibersort, fgsea_pax_neg_14to20_vs_adult_cibersort))
modules_pax_peds_vs_adults$pathway <- str_to_sentence(modules_pax_peds_vs_adults$pathway)
modules_pax_peds_vs_adults$pathway[modules_pax_peds_vs_adults$pathway=="Bcr signaling"] <- "B cell receptor signaling"
modules_pax_peds_vs_adults$pathway[modules_pax_peds_vs_adults$pathway=="Nk activity"] <- "NK cell activity"
modules_pax_peds_vs_adults$pathway[modules_pax_peds_vs_adults$pathway=="Tcr signaling"] <- "T cell receptor signaling"
modules_pax_peds_vs_adults$pathway[modules_pax_peds_vs_adults$pathway=="Type i interferon signaling"] <- "Type I interferon signaling"
modules_pax_peds_vs_adults$pathway[modules_pax_peds_vs_adults$pathway=="Type ii interferon signaling"] <- "Type II interferon signaling"
modules_pax_peds_vs_adults$pathway[modules_pax_peds_vs_adults$pathway=="Innate immune cell activation"] <- "Innate immune activation"
modules_pax_peds_vs_adults$pathway[modules_pax_peds_vs_adults$pathway=="Tnf signaling"] <- "TNF signaling"
modules_pax_peds_vs_adults$pathway[modules_pax_peds_vs_adults$pathway=="Tlr signaling"] <- "Toll-like receptor signaling"
modules_pax_peds_vs_adults$pathway[modules_pax_peds_vs_adults$pathway=="Treg differentiation"] <- "Treg cell differentiation"
modules_pax_peds_vs_adults$Sample <- factor(modules_pax_peds_vs_adults$Sample, levels=c("0-5 yr vs. Adult","6-13 yr vs. Adult", "14-20 yr vs. Adult"))

if (is.element('Up-regulated', modules_pax_peds_vs_adults$Expression) & is.element('Unchanged', modules_pax_peds_vs_adults$Expression) & is.element('Down-regulated', modules_pax_peds_vs_adults$Expression)) {
  cen = c("lightsteelblue2", "gray80", "#F9AAAE")
} else if ( is.element('Up-regulated', modules_pax_peds_vs_adults$Expression) & is.element('Unchanged', modules_pax_peds_vs_adults$Expression)) {
  cen = c("gray80", "#F9AAAE")
} else if ( is.element('Unchanged', modules_pax_peds_vs_adults$Expression) & is.element('Down-regulated', modules_pax_peds_vs_adults$Expression)) {
  cen = c("lightsteelblue2", "gray80")
} else if ( is.element('Up-regulated', modules_pax_peds_vs_adults$Expression) & is.element('Down-regulated', modules_pax_peds_vs_adults$Expression)) {
  cen = c("lightsteelblue2", "#F9AAAE")
} else {cen = c("gray80")}

modules_pax_peds_vs_adults_innate <- subset(modules_pax_peds_vs_adults, pathway %in% innate)
modules_pax_peds_vs_adults_innate$pathway <- factor(modules_pax_peds_vs_adults_innate$pathway, levels=c("Coagulation", "Complement system", "NK cell activity", "Inflammasomes", "Innate immune activation", 
                                                                          "Myeloid inflammation", "Myeloid activation", "Phagocytosis", "Toll-like receptor signaling", "TNF signaling",
                                                                          "Type II interferon signaling", "Type I interferon signaling", "Interferon response"))

innate_plot_pax <- ggplot(modules_pax_peds_vs_adults_innate, aes(x = pathway, y = Sample)) +
  geom_point(aes(size = (abs(NES))^3, fill = Expression), shape = 21) +
  geom_text(aes(label = formatC(NES, format='f', digits = 2)), parse = FALSE) +
  scale_fill_manual(values = cen) +
  labs(x= "Innate immunity", y = "") + 
  theme(legend.key=element_blank(), legend.position = "none", panel.background = element_rect(colour = "black", fill = NA, linewidth=0.8),
        axis.title.y = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_blank(), axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        plot.title=element_text(size=14, color="black", hjust=0.5, face="bold")) + 
  scale_size_continuous(range = c(13,15)) + ggtitle("Peripheral blood") +
  coord_flip() 
innate_plot_pax

modules_pax_peds_vs_adults_adaptive <- subset(modules_pax_peds_vs_adults, pathway %in% adaptive)
modules_pax_peds_vs_adults_adaptive$pathway <- factor(modules_pax_peds_vs_adults_adaptive$pathway, levels=c("Treg cell differentiation", "Lymphocyte trafficking", "Mononuclear cell migration", 
                                                                              "T cell receptor signaling", "B cell receptor signaling", "Adaptive immune response"))

adaptive_plot_pax <- ggplot(modules_pax_peds_vs_adults_adaptive, aes(x = pathway, y = Sample)) +
  geom_point(aes(size = (abs(NES))^3, fill = Expression), shape = 21) +
  geom_text(aes(label = formatC(NES, format='f', digits = 2)), parse = FALSE) +
  scale_fill_manual(values = cen) +
  labs(x= "Adaptive immunity", y = "") + 
  theme(legend.key=element_blank(), legend.position = "none", panel.background = element_rect(colour = "black", fill = NA, linewidth=0.8),
        axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=12, color = "black", angle=45, hjust=1)) + 
  scale_size_continuous(range = c(12.5,14)) +
  coord_flip() 
adaptive_plot_pax

fig1_modules_adult <- plot_grid(innate_plot_pax, NULL, adaptive_plot_pax, labels=NULL, ncol=1, rel_heights=c(4.05,-0.06,2.5), align="v") 
title_adult <- ggdraw() + draw_label("Pediatric vs. adult", size=14, fontface='bold')
fig1_modules_adult_title <- plot_grid(title_adult, fig1_modules_adult, labels=NULL, ncol=1, rel_heights=c(0.03,0.97), align="v")

png(file="Statistical_Analyses/Figures/Figure_1.png", width = 15, height = 11, units = 'in', res = 1200)
plot_grid(fig1_cibersort, NULL, fig1_modules_peds_title, NULL, fig1_modules_adult_title, labels=c("","","c","d",""), ncol=5, rel_widths=c(8,0.1,7.8,0.15,2.5), align="h") 
dev.off()