# BRAVE Kids RNA Sequencing Analysis
# Aditya Mohan (MD/PhD candidate)  / Matthew Kelly, MD, MPH 
# Figure 1
# Last update: Oct. 6, 2023

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
phy.rnaseq.np.Neg <- subset_samples(phy.rnaseq.np, corona=="Negative")
metadata_np_neg <- data.frame(sample_data(phy.rnaseq.np.Neg))
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
phy.rnaseq.pax.Neg <- subset_samples(phy.rnaseq.pax, corona=="Negative")
metadata_pax_neg <- data.frame(sample_data(phy.rnaseq.pax.Neg))
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
fgsea_np_Neg_0to5_6to13 <- data.frame(read_excel("Statistical_Analyses/1_COVID_Neg_by_Age/modules_np_Neg_0to5_6to13.xlsx"))
fgsea_np_Neg_0to5_14to20 <- data.frame(read_excel("Statistical_Analyses/1_COVID_Neg_by_Age/modules_np_Neg_0to5_14to20.xlsx"))
fgsea_np_Neg_6to13_14to20 <- data.frame(read_excel("Statistical_Analyses/1_COVID_Neg_by_Age/modules_np_Neg_6to13_14to20.xlsx"))
fgsea_pax_Neg_0to5_6to13 <- data.frame(read_excel("Statistical_Analyses/1_COVID_Neg_by_Age/modules_pax_Neg_0to5_6to13.xlsx"))
fgsea_pax_Neg_0to5_14to20 <- data.frame(read_excel("Statistical_Analyses/1_COVID_Neg_by_Age/modules_pax_Neg_0to5_14to20.xlsx"))
fgsea_pax_Neg_6to13_14to20 <- data.frame(read_excel("Statistical_Analyses/1_COVID_Neg_by_Age/modules_pax_Neg_6to13_14to20.xlsx"))
fgsea_pax_Neg_0to5_Adult <- data.frame(read_excel("Statistical_Analyses/1_COVID_Neg_by_Age/modules_pax_Neg_0to5_Adult.xlsx"))
fgsea_pax_Neg_6to13_Adult <- data.frame(read_excel("Statistical_Analyses/1_COVID_Neg_by_Age/modules_pax_Neg_6to13_Adult.xlsx"))
fgsea_pax_Neg_14to20_Adult <- data.frame(read_excel("Statistical_Analyses/1_COVID_Neg_by_Age/modules_pax_Neg_14to20_Adult.xlsx"))

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

# ANALYSES ADJUSTING FOR CIBERSORT CELL POPULATIONS

# FGSEA modules dotplot - PEDIATRIC AGE GROUPS

modules_NP <- rbindlist(list(fgsea_np_Neg_0to5_6to13, fgsea_np_Neg_0to5_14to20, fgsea_np_Neg_6to13_14to20))
modules_NP$pathway <- str_to_sentence(modules_NP$pathway)
modules_NP$pathway[modules_NP$pathway=="Bcr signaling"] <- "B cell receptor signaling"
modules_NP$pathway[modules_NP$pathway=="Nk activity"] <- "NK cell activity"
modules_NP$pathway[modules_NP$pathway=="Tcr signaling"] <- "T cell receptor signaling"
modules_NP$pathway[modules_NP$pathway=="Type i interferon signaling"] <- "Type I interferon signaling"
modules_NP$pathway[modules_NP$pathway=="Type ii interferon signaling"] <- "Type II interferon signaling"
modules_NP$Sample[modules_NP$Sample=="6-13 vs. 14-20 yr"] <- " 6-13 vs. 14-20 yr"
modules_NP$Sample <- factor(modules_NP$Sample, levels=c("0-5 vs. 6-13 yr","0-5 vs. 14-20 yr", " 6-13 vs. 14-20 yr"))

if (is.element('Up-regulated', modules_NP$Expression) & is.element('Unchanged', modules_NP$Expression) & is.element('Down-regulated', modules_NP$Expression)) {
  cen = c("lightsteelblue2", "gray80", "#F9AAAE")
} else if ( is.element('Up-regulated', modules_NP$Expression) & is.element('Unchanged', modules_NP$Expression)) {
  cen = c("gray80", "#F9AAAE")
} else if ( is.element('Unchanged', modules_NP$Expression) & is.element('Down-regulated', modules_NP$Expression)) {
  cen = c("lightsteelblue2", "gray80")
} else if ( is.element('Up-regulated', modules_NP$Expression) & is.element('Down-regulated', modules_NP$Expression)) {
  cen = c("lightsteelblue2", "#F9AAAE")
} else {cen = c("gray80")}

innate <- c("Coagulation", "Complement system", "NK cell activity", "Inflammasomes", "Innate immune activation", "Myeloid activation", "Myeloid inflammation", "Phagocytosis", 
            "Type II interferon signaling", "Type I interferon signaling", "Interferon response")
modules_NP_innate <- subset(modules_NP, pathway %in% innate)
modules_NP_innate$pathway <- factor(modules_NP_innate$pathway, levels=c("Coagulation", "Complement system", "NK cell activity", "Inflammasomes", "Innate immune activation", "Myeloid activation", 
                                                                        "Myeloid inflammation", "Phagocytosis", "Type II interferon signaling",
                                                                        "Type I interferon signaling", "Interferon response"))

innate_plot_np <- ggplot(modules_NP_innate, aes(x = pathway, y = Sample)) +
  geom_point(aes(size = (abs(NES))^3, fill = Expression), shape = 21) +
  geom_text(aes(label = formatC(NES, format='f', digits = 2)), parse = FALSE) +
  scale_fill_manual(values = cen) +
  labs(x= "Innate immunity", y = "") + 
  theme(legend.key=element_blank(), legend.position = "none", panel.background = element_rect(colour = "black", fill = NA, linewidth=0.8),
        axis.title.y = element_text(size=14, face="bold", margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12, color = "black"), axis.text.x = element_blank(), axis.ticks = element_blank(),
        plot.title=element_text(size=14, color="black", hjust=0.5, face="bold")) + 
  scale_size_continuous(range = c(13.5, 18.5)) + ggtitle("Upper respiratory") +
  coord_flip() 
innate_plot_np

adaptive <- c("Lymphocyte trafficking", "Mononuclear cell migration", "T cell receptor signaling", "B cell receptor signaling", 
              "Adaptive immune response")
modules_NP_adaptive <- subset(modules_NP, pathway %in% adaptive)
modules_NP_adaptive$pathway <- factor(modules_NP_adaptive$pathway, levels=c("Lymphocyte trafficking", "Mononuclear cell migration",  
                                                                            "T cell receptor signaling", "B cell receptor signaling", "Adaptive immune response"))

adaptive_plot_np <- ggplot(modules_NP_adaptive, aes(x = pathway, y = Sample)) +
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
  scale_size_continuous(range = c(13.5, 18.5)) +
  coord_flip() 
adaptive_plot_np

modules_Pax <- rbindlist(list(fgsea_pax_Neg_0to5_6to13, fgsea_pax_Neg_0to5_14to20, fgsea_pax_Neg_6to13_14to20))
modules_Pax$pathway <- str_to_sentence(modules_Pax$pathway)
modules_Pax$pathway[modules_Pax$pathway=="Bcr signaling"] <- "B cell receptor signaling"
modules_Pax$pathway[modules_Pax$pathway=="Nk activity"] <- "NK cell activity"
modules_Pax$pathway[modules_Pax$pathway=="Tcr signaling"] <- "T cell receptor signaling"
modules_Pax$pathway[modules_Pax$pathway=="Type i interferon signaling"] <- "Type I interferon signaling"
modules_Pax$pathway[modules_Pax$pathway=="Type ii interferon signaling"] <- "Type II interferon signaling"
modules_Pax$Sample[modules_Pax$Sample=="6-13 vs. 14-20 yr"] <- " 6-13 vs. 14-20 yr"
modules_Pax$Sample <- factor(modules_Pax$Sample, levels=c("0-5 vs. 6-13 yr","0-5 vs. 14-20 yr", " 6-13 vs. 14-20 yr"))

if (is.element('Up-regulated', modules_Pax$Expression) & is.element('Unchanged', modules_Pax$Expression) & is.element('Down-regulated', modules_Pax$Expression)) {
  cen = c("lightsteelblue2", "gray80", "#F9AAAE")
} else if ( is.element('Up-regulated', modules_Pax$Expression) & is.element('Unchanged', modules_Pax$Expression)) {
  cen = c("gray80", "#F9AAAE")
} else if ( is.element('Unchanged', modules_Pax$Expression) & is.element('Down-regulated', modules_Pax$Expression)) {
  cen = c("lightsteelblue2", "gray80")
} else if ( is.element('Up-regulated', modules_Pax$Expression) & is.element('Down-regulated', modules_Pax$Expression)) {
  cen = c("lightsteelblue2", "#F9AAAE")
} else {cen = c("gray80")}

modules_Pax_innate <- subset(modules_Pax, pathway %in% innate)
modules_Pax_innate$pathway <- factor(modules_Pax_innate$pathway, levels=c("Coagulation", "Complement system", "NK cell activity", "Inflammasomes", "Innate immune activation", "Myeloid activation", 
                                                                          "Myeloid inflammation", "Phagocytosis", "Type II interferon signaling",
                                                                          "Type I interferon signaling", "Interferon response"))

innate_plot_pax <- ggplot(modules_Pax_innate, aes(x = pathway, y = Sample)) +
  geom_point(aes(size = (abs(NES))^3, fill = Expression), shape = 21) +
  geom_text(aes(label = formatC(NES, format='f', digits = 2)), parse = FALSE) +
  scale_fill_manual(values = cen) +
  labs(x= "Innate immunity", y = "") + 
  theme(legend.key=element_blank(), legend.position = "none", panel.background = element_rect(colour = "black", fill = NA, linewidth=0.8),
        axis.title.y = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_blank(), axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        plot.title=element_text(size=14, color="black", hjust=0.5, face="bold")) + 
  scale_size_continuous(range = c(13.5, 18.5)) + ggtitle("Peripheral blood") +
  coord_flip() 
innate_plot_pax

modules_Pax_adaptive <- subset(modules_Pax, pathway %in% adaptive)
modules_Pax_adaptive$pathway <- factor(modules_Pax_adaptive$pathway, levels=c("Lymphocyte trafficking", "Mononuclear cell migration", 
                                                                              "T cell receptor signaling", "B cell receptor signaling", "Adaptive immune response"))

adaptive_plot_pax <- ggplot(modules_Pax_adaptive, aes(x = pathway, y = Sample)) +
  geom_point(aes(size = (abs(NES))^3, fill = Expression), shape = 21) +
  geom_text(aes(label = formatC(NES, format='f', digits = 2)), parse = FALSE) +
  scale_fill_manual(values = cen) +
  labs(x= "Adaptive immunity", y = "") + 
  theme(legend.key=element_blank(), legend.position = "none", panel.background = element_rect(colour = "black", fill = NA, linewidth=0.8),
        axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=12, color = "black", angle=45, hjust=1)) + 
  scale_size_continuous(range = c(13.5, 18.5)) +
  coord_flip() 
adaptive_plot_pax

module_plot_np <- plot_grid(innate_plot_np, NULL, adaptive_plot_np, labels=NULL, ncol=1, rel_heights=c(3.8,-0.06,2.55), align="v") 
module_plot_pax <- plot_grid(innate_plot_pax, NULL, adaptive_plot_pax, labels=NULL, ncol=1, rel_heights=c(3.8,-0.06,2.55), align="v") 
fig1_modules_peds <- plot_grid(module_plot_np, NULL, module_plot_pax, labels=NULL, ncol=3, rel_widths=c(3.6,-0.05, 1.8), align="vh") 
title_peds <- ggdraw() + draw_label("                                            Pediatric age groups", size=14, fontface='bold')
fig1_modules_peds_title <- plot_grid(title_peds, fig1_modules_peds, labels=NULL, ncol=1, rel_heights=c(0.03,0.97), align="v")

# FGSEA modules dotplot - COMPARISONS TO ADULTS (PERIPHERAL BLOOD ONLY)

modules_Pax <- rbindlist(list(fgsea_pax_Neg_0to5_Adult, fgsea_pax_Neg_6to13_Adult, fgsea_pax_Neg_14to20_Adult))
modules_Pax$pathway <- str_to_sentence(modules_Pax$pathway)
modules_Pax$pathway[modules_Pax$pathway=="Bcr signaling"] <- "B cell receptor signaling"
modules_Pax$pathway[modules_Pax$pathway=="Nk activity"] <- "NK cell activity"
modules_Pax$pathway[modules_Pax$pathway=="Tcr signaling"] <- "T cell receptor signaling"
modules_Pax$pathway[modules_Pax$pathway=="Type i interferon signaling"] <- "Type I interferon signaling"
modules_Pax$pathway[modules_Pax$pathway=="Type ii interferon signaling"] <- "Type II interferon signaling"
modules_Pax$Sample <- factor(modules_Pax$Sample, levels=c("0-5 yr vs. Adult","6-13 yr vs. Adult", "14-20 yr vs. Adult"))

if (is.element('Up-regulated', modules_Pax$Expression) & is.element('Unchanged', modules_Pax$Expression) & is.element('Down-regulated', modules_Pax$Expression)) {
  cen = c("lightsteelblue2", "gray80", "#F9AAAE")
} else if ( is.element('Up-regulated', modules_Pax$Expression) & is.element('Unchanged', modules_Pax$Expression)) {
  cen = c("gray80", "#F9AAAE")
} else if ( is.element('Unchanged', modules_Pax$Expression) & is.element('Down-regulated', modules_Pax$Expression)) {
  cen = c("lightsteelblue2", "gray80")
} else if ( is.element('Up-regulated', modules_Pax$Expression) & is.element('Down-regulated', modules_Pax$Expression)) {
  cen = c("lightsteelblue2", "#F9AAAE")
} else {cen = c("gray80")}

modules_Pax_innate <- subset(modules_Pax, pathway %in% innate)
modules_Pax_innate$pathway <- factor(modules_Pax_innate$pathway, levels=c("Coagulation", "Complement system", "NK cell activity", "Inflammasomes", "Innate immune activation", "Myeloid activation", 
                                                                          "Myeloid inflammation", "Phagocytosis", "Type II interferon signaling",
                                                                          "Type I interferon signaling", "Interferon response"))

innate_plot_pax <- ggplot(modules_Pax_innate, aes(x = pathway, y = Sample)) +
  geom_point(aes(size = (abs(NES))^3, fill = Expression), shape = 21) +
  geom_text(aes(label = formatC(NES, format='f', digits = 2)), parse = FALSE) +
  scale_fill_manual(values = cen) +
  labs(x= "Innate immunity", y = "") + 
  theme(legend.key=element_blank(), legend.position = "none", panel.background = element_rect(colour = "black", fill = NA, linewidth=0.8),
        axis.title.y = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_blank(), axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        plot.title=element_text(size=14, color="black", hjust=0.5, face="bold")) + 
  scale_size_continuous(range = c(13.5, 18.5)) + ggtitle("Peripheral blood") +
  coord_flip() 
innate_plot_pax

modules_Pax_adaptive <- subset(modules_Pax, pathway %in% adaptive)
modules_Pax_adaptive$pathway <- factor(modules_Pax_adaptive$pathway, levels=c("Lymphocyte trafficking", "Mononuclear cell migration", 
                                                                              "T cell receptor signaling", "B cell receptor signaling", "Adaptive immune response"))

adaptive_plot_pax <- ggplot(modules_Pax_adaptive, aes(x = pathway, y = Sample)) +
  geom_point(aes(size = (abs(NES))^3, fill = Expression), shape = 21) +
  geom_text(aes(label = formatC(NES, format='f', digits = 2)), parse = FALSE) +
  scale_fill_manual(values = cen) +
  labs(x= "Adaptive immunity", y = "") + 
  theme(legend.key=element_blank(), legend.position = "none", panel.background = element_rect(colour = "black", fill = NA, linewidth=0.8),
        axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=12, color = "black", angle=45, hjust=1)) + 
  scale_size_continuous(range = c(13.5, 18.5)) +
  coord_flip() 
adaptive_plot_pax

fig1_modules_adult <- plot_grid(innate_plot_pax, NULL, adaptive_plot_pax, labels=NULL, ncol=1, rel_heights=c(3.8,-0.06,2.55), align="v") 
title_adult <- ggdraw() + draw_label("Pediatric vs. adult", size=14, fontface='bold')
fig1_modules_adult_title <- plot_grid(title_adult, fig1_modules_adult, labels=NULL, ncol=1, rel_heights=c(0.03,0.97), align="v")

# NO TITLE

#png(file="Statistical_Analyses/Figures/Figure_1.png", width = 16.5, height = 11, units = 'in', res = 1200)
#plot_grid(fig1_cibersort, NULL, fig1_modules_peds_title, NULL, fig1_modules_adult_title, labels=c("a","","b","c",""), ncol=5, rel_widths=c(8,0.1,8,0.18,2.7), align="h") 
#dev.off()

# ADDITION OF TITLE

fig1_plot <- plot_grid(fig1_cibersort, NULL, fig1_modules_peds_title, NULL, fig1_modules_adult_title, labels=c("","","c","d",""), ncol=5, rel_widths=c(8,0.1,8,0.18,2.7), align="h") 
fig1_title <- ggdraw() + draw_label("Figure 1. Transcriptional profiles within the upper respiratory tract and peripheral blood of healthy children, adolescents, and adults", 
                                    size=15, fontface='bold')

png(file="Statistical_Analyses/Figures/Figure_1.png", width = 16.5, height = 11.75, units = 'in', res = 1200)
plot_grid(fig1_title, NULL, fig1_plot, labels=NULL, ncol=1, rel_heights=c(0.03,0.02,0.97), align="v") 
dev.off()