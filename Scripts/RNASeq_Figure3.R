# BRAVE Kids RNA Sequencing Analysis
# Aditya Mohan (MD/PhD candidate) / Matthew Kelly, MD, MPH 
# Figure 3 - cell populations and gene module expression among SARS-CoV-2-infected children, adolescents, and adults
# Analyses of gene module expression in peripheral blood samples adjusted for imputed cell proportions
# Last update: February 22, 2025

remove(list=ls())
setwd("___________________") 
set.seed(1234)
getRversion()

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
library(tibble)
library(msigdbr)
library(reshape2)
library(data.table)
library(RColorBrewer)
library(stringr)
library(phyloseq)

# Upload files with CIBERSORTx data
phy.rnaseq.np <- readRDS("phy.rnaseq.np.rds")
phy.rnaseq.np <- subset_samples(phy.rnaseq.np, corona=="Positive")
metadata_np_pos <- data.frame(sample_data(phy.rnaseq.np))
metadata_np_pos_pops <- metadata_np_pos[,c("age_cat","B_cells","Plasma_cells", "T_cells_CD8","T_cells_CD4","T_cells_gamma_delta","NK_cells",
                                           "Mono_Macrophages","Dendritic_cells","Mast_cells","Eosinophils","Neutrophils")]
metadata_np_pos_pops$age_cat <- factor(metadata_np_pos_pops$age_cat, levels=c("0-5 years", "6-13 years", "14-20 years"))
metadata_np_pos_pops <- gather(metadata_np_pos_pops, cell_type, value, B_cells:Neutrophils, factor_key=TRUE)
metadata_np_pos_pops$cell_type2[metadata_np_pos_pops$cell_type=="B_cells"] <- "B cells"
metadata_np_pos_pops$cell_type2[metadata_np_pos_pops$cell_type=="Plasma_cells"] <- "Plasma cells"
metadata_np_pos_pops$cell_type2[metadata_np_pos_pops$cell_type=="T_cells_CD4"] <- "CD4+ T cells"
metadata_np_pos_pops$cell_type2[metadata_np_pos_pops$cell_type=="T_cells_CD8"] <- "CD8+ T cells"
metadata_np_pos_pops$cell_type2[metadata_np_pos_pops$cell_type=="T_cells_gamma_delta"] <- "Gamma delta T cells"
metadata_np_pos_pops$cell_type2[metadata_np_pos_pops$cell_type=="NK_cells"] <- "NK cells"
metadata_np_pos_pops$cell_type2[metadata_np_pos_pops$cell_type=="Mono_Macrophages"] <- "Mono/macrophages"
metadata_np_pos_pops$cell_type2[metadata_np_pos_pops$cell_type=="Dendritic_cells"] <- "Dendritic cells"
metadata_np_pos_pops$cell_type2[metadata_np_pos_pops$cell_type=="Mast_cells"] <- "Mast cells"
metadata_np_pos_pops$cell_type2[metadata_np_pos_pops$cell_type=="Eosinophils"] <- "Eosinophils"
metadata_np_pos_pops$cell_type2[metadata_np_pos_pops$cell_type=="Neutrophils"] <- "Neutrophils"
tapply(metadata_np_pos_pops$value, metadata_np_pos_pops$cell_type2, summary) 
# Exclude cell types identified in <25% of samples from figure
metadata_np_pos_pops <- subset(metadata_np_pos_pops, cell_type2!="Eosinophils")
metadata_np_pos_pops <- subset(metadata_np_pos_pops, cell_type2!="Gamma delta T cells")

phy.rnaseq.pax <- readRDS("phy.rnaseq.pax.rds")
phy.rnaseq.pax.Pos <- subset_samples(phy.rnaseq.pax, corona=="Positive")
metadata_pax_pos <- data.frame(sample_data(phy.rnaseq.pax.Pos))
metadata_pax_pos_pops <- metadata_pax_pos[,c("age_cat","B_cells","Plasma_cells", "T_cells_CD8","T_cells_CD4","T_cells_gamma_delta","NK_cells",
                                             "Mono_Macrophages","Dendritic_cells","Mast_cells","Eosinophils","Neutrophils")]
metadata_pax_pos_pops$age_cat <- factor(metadata_pax_pos_pops$age_cat, levels=c("0-5 years", "6-13 years", "14-20 years", "Adult"))
metadata_pax_pos_pops <- gather(metadata_pax_pos_pops, cell_type, value, B_cells:Neutrophils, factor_key=TRUE)
metadata_pax_pos_pops$cell_type2[metadata_pax_pos_pops$cell_type=="B_cells"] <- "B cells"
metadata_pax_pos_pops$cell_type2[metadata_pax_pos_pops$cell_type=="Plasma_cells"] <- "Plasma cells"
metadata_pax_pos_pops$cell_type2[metadata_pax_pos_pops$cell_type=="T_cells_CD4"] <- "CD4+ T cells"
metadata_pax_pos_pops$cell_type2[metadata_pax_pos_pops$cell_type=="T_cells_CD8"] <- "CD8+ T cells"
metadata_pax_pos_pops$cell_type2[metadata_pax_pos_pops$cell_type=="T_cells_gamma_delta"] <- "Gamma delta T cells"
metadata_pax_pos_pops$cell_type2[metadata_pax_pos_pops$cell_type=="NK_cells"] <- "NK cells"
metadata_pax_pos_pops$cell_type2[metadata_pax_pos_pops$cell_type=="Mono_Macrophages"] <- "Mono/macrophages"
metadata_pax_pos_pops$cell_type2[metadata_pax_pos_pops$cell_type=="Dendritic_cells"] <- "Dendritic cells"
metadata_pax_pos_pops$cell_type2[metadata_pax_pos_pops$cell_type=="Mast_cells"] <- "Mast cells"
metadata_pax_pos_pops$cell_type2[metadata_pax_pos_pops$cell_type=="Eosinophils"] <- "Eosinophils"
metadata_pax_pos_pops$cell_type2[metadata_pax_pos_pops$cell_type=="Neutrophils"] <- "Neutrophils"
tapply(metadata_pax_pos_pops$value, metadata_pax_pos_pops$cell_type2, summary) 
# Exclude cell types identified in <25% of samples from figure
metadata_pax_pos_pops <- subset(metadata_pax_pos_pops, cell_type2!="Dendritic cells" & cell_type2!="Eosinophils" & cell_type2!="Gamma delta T cells")

# ANALYSES ADJUSTING FOR CIBERSORT CELL POPULATIONS

# Upload FGSEA output files
fsgea_np_0to5_pos_neg_nocibersort <- data.frame(read_excel("Statistical_Analyses/3_COVID_Pos_by_Age/modules_np_0to5_pos_neg_nocibersort.xlsx"))
fsgea_np_6to13_pos_neg_nocibersort <- data.frame(read_excel("Statistical_Analyses/3_COVID_Pos_by_Age/modules_np_6to13_pos_neg_nocibersort.xlsx"))
fsgea_np_14to20_pos_neg_nocibersort <- data.frame(read_excel("Statistical_Analyses/3_COVID_Pos_by_Age/modules_np_14to20_pos_neg_nocibersort.xlsx"))
fsgea_pax_0to5_pos_neg_cibersort <- data.frame(read_excel("Statistical_Analyses/3_COVID_Pos_by_Age/modules_pax_0to5_pos_neg_cibersort.xlsx"))
fsgea_pax_6to13_pos_neg_cibersort <- data.frame(read_excel("Statistical_Analyses/3_COVID_Pos_by_Age/modules_pax_6to13_pos_neg_cibersort.xlsx"))
fsgea_pax_14to20_pos_neg_cibersort <- data.frame(read_excel("Statistical_Analyses/3_COVID_Pos_by_Age/modules_pax_14to20_pos_neg_cibersort.xlsx"))
fsgea_pax_adult_pos_neg_cibersort <- data.frame(read_excel("Statistical_Analyses/3_COVID_Pos_by_Age/modules_pax_adult_pos_neg_cibersort.xlsx"))

# CIBERSORTx cell populations barplot

# *, p<0.05
# **, p<0.01
# ***, p<0.001
# ****, p<0.0001

# UPPER RESPIRATORY SAMPLES (COVID-POSITIVE BY AGE)
# Adjust comparisons for multiple testing using BH (n=9)
# B cells: 0.003
# Plasma cells: >0.99
# CD8+ T cells: 0.11
# CD4+ T cells: >0.99
# NK cells: >0.99
# Mono/macrophages: 0.34
# Dendritic cells: >0.99
# Mast cells: >0.99
# Neutrophils: >0.99

cibersort_np_pos <- ggplot(metadata_np_pos_pops, aes(x=cell_type2, y=value, fill=age_cat)) + geom_boxplot() + ylab("Imputed immune cell proportions") + 
  theme(legend.title = element_blank(), legend.position = "right", legend.text = element_text(size=12), legend.box.spacing = unit(3, "pt"), legend.key=element_rect(colour="white"),
        panel.background = element_rect(colour = "black", fill = NA, linewidth=0.5), panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey70"),
        panel.grid.minor.y = element_line(linewidth = 0.2, colour = "grey70"),
        axis.title.y = element_text(size=13, face="bold", margin = ggplot2::margin(t = 0, r = 5, b = 0, l = 0)),
        axis.title.x = element_blank(), axis.text.y = element_text(size=11.5), axis.text.x = element_text(size=11.5, color = "black", angle=35, hjust=1), 
        plot.title = element_text(size=14, hjust = 0.5, face="bold")) + scale_fill_brewer(palette="Paired") +
  scale_y_continuous(limits=c(0,0.61), breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6)) +
  ggtitle("Upper respiratory") + 
  annotate("text", x=1, y=0.61, size=6, label= "**") # B cells (beta regression, p=0.003) 

# PERIPHERAL BLOOD SAMPLES (COVID-POSITIVE BY AGE)
# Adjust comparisons for multiple testing using BH (n=9)
# B cells: 6.66E-19
# Plasma cells: 0.0011
# CD8+ T cells: 0.0003
# CD4+ T cells: 0.09
# NK cells: 0.02
# Mono/macrophages: 1.52E-11
# Dendritic cells: >0.99 
# Mast cells: >0.99
# Neutrophils: 2.08E-05

cibersort_pax_pos <- ggplot(metadata_pax_pos_pops, aes(x=cell_type2, y=value, fill=age_cat)) + geom_boxplot() + ylab("Imputed immune cell proportions") + 
  theme(legend.title = element_blank(), legend.position = "right", legend.text = element_text(size=12), legend.box.spacing = unit(3, "pt"), legend.key=element_rect(colour="white"),
        panel.background = element_rect(colour = "black", fill = NA, linewidth=0.5), panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey70"),
        panel.grid.minor.y = element_line(linewidth = 0.2, colour = "grey70"),
        axis.title.y = element_text(size=13, face="bold", margin = ggplot2::margin(t = 0, r = 5, b = 0, l = 0)),
        axis.title.x = element_blank(), axis.text.y = element_text(size=11.5), axis.text.x = element_text(size=11.5, color = "black", angle=35, hjust=1), 
        plot.title = element_text(size=14, hjust = 0.5, face="bold")) + scale_fill_brewer(palette="Paired") +
  scale_y_continuous(limits=c(0,0.51), breaks=c(0,0.1,0.2,0.3,0.4,0.5)) +
  ggtitle("Peripheral blood") + 
  annotate("text", x=1, y=0.51, size=6, label= "****") + # B cells (beta regression, p=6.66E-19) 
  annotate("text", x=8, y=0.51, size=6, label= "**") + # Plasma cells (beta regression, p=0.0011)
  annotate("text", x=3, y=0.51, size=6, label= "***") + # CD8+ T cells (beta regression, p=0.0003)
  annotate("text", x=7, y=0.51, size=6, label= "*") + # NK cells (beta regression, p=0.02) 
  annotate("text", x=5, y=0.51, size=6, label= "****") + # Mono/macrophages (beta regression, p=1.52E-11) 
  annotate("text", x=6, y=0.51, size=6, label= "****") # Neutrophils (beta regression, p=2.08E-05) 

fig3_cibersort <- plot_grid(NULL, cibersort_np_pos, NULL, cibersort_pax_pos, ncol=1, labels=c("","a","","b"), rel_heights=c(0.05,3,0.05,3), align="v") 
fig3_cibersort

# FGSEA modules dotplot

modules_np <- rbindlist(list(fsgea_np_0to5_pos_neg_nocibersort, fsgea_np_6to13_pos_neg_nocibersort, fsgea_np_14to20_pos_neg_nocibersort))
modules_np$pathway <- str_to_sentence(modules_np$pathway)
modules_np$pathway[modules_np$pathway=="Bcr signaling"] <- "B cell receptor signaling"
modules_np$pathway[modules_np$pathway=="Nk activity"] <- "NK cell activity"
modules_np$pathway[modules_np$pathway=="Tcr signaling"] <- "T cell receptor signaling"
modules_np$pathway[modules_np$pathway=="Type i interferon signaling"] <- "Type I interferon signaling"
modules_np$pathway[modules_np$pathway=="Type ii interferon signaling"] <- "Type II interferon signaling"
modules_np$pathway[modules_np$pathway=="Innate immune cell activation"] <- "Innate immune activation"
modules_np$pathway[modules_np$pathway=="Tnf signaling"] <- "TNF signaling"
modules_np$pathway[modules_np$pathway=="Tlr signaling"] <- "Toll-like receptor signaling"
modules_np$pathway[modules_np$pathway=="Treg differentiation"] <- "Treg cell differentiation"
modules_np$Sample <- factor(modules_np$Sample, levels=c("0-5 yr", "6-13 yr", "14-20 yr"))

if (is.element('Up-regulated', modules_np$Expression) & is.element('Unchanged', modules_np$Expression) & is.element('Down-regulated', modules_np$Expression)) {
  cen = c("lightsteelblue2", "gray80", "#F9AAAE")
} else if ( is.element('Up-regulated', modules_np$Expression) & is.element('Unchanged', modules_np$Expression)) {
  cen = c("gray80", "#F9AAAE")
} else if ( is.element('Unchanged', modules_np$Expression) & is.element('Down-regulated', modules_np$Expression)) {
  cen = c("lightsteelblue2", "gray80")
} else if ( is.element('Up-regulated', modules_np$Expression) & is.element('Down-regulated', modules_np$Expression)) {
  cen = c("lightsteelblue2", "#F9AAAE")
} else {cen = c("gray80")}

innate <- c("Coagulation", "Complement system", "NK cell activity", "Inflammasomes", "Innate immune activation", 
            "Myeloid inflammation", "Myeloid activation", "Phagocytosis", "Toll-like receptor signaling", "TNF signaling",
            "Type II interferon signaling", "Type I interferon signaling", "Interferon response")
modules_np_innate <- subset(modules_np, pathway %in% innate)
modules_np_innate$pathway <- factor(modules_np_innate$pathway, levels=c("Coagulation", "Complement system", "NK cell activity", "Inflammasomes", "Innate immune activation", 
                                                                        "Myeloid inflammation", "Myeloid activation", "Phagocytosis", "Toll-like receptor signaling", "TNF signaling",
                                                                        "Type II interferon signaling", "Type I interferon signaling", "Interferon response"))

innate_plot_np <- ggplot(modules_np_innate, aes(x = pathway, y = Sample)) +
  geom_point(aes(size = (abs(NES))^3, fill = Expression), shape = 21) +
  geom_text(aes(label = formatC(NES, format='f', digits = 2)), parse = FALSE) +
  scale_fill_manual(values = cen) +
  labs(x= "Innate immunity", y = "") + 
  theme(legend.key=element_blank(), legend.position = "none", panel.background = element_rect(colour = "black", fill = NA, linewidth=0.8),
        axis.title.y = element_text(size=14, face="bold", margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12, color = "black"), axis.text.x = element_blank(), axis.ticks = element_blank(),
        plot.title=element_text(size=14, color="black", hjust=0.5, face="bold")) + 
  scale_size_continuous(range = c(13,16)) + ggtitle("Upper respiratory") +
  coord_flip() 
innate_plot_np

adaptive <- c("Treg cell differentiation", "Lymphocyte trafficking", "Mononuclear cell migration", 
              "T cell receptor signaling", "B cell receptor signaling", "Adaptive immune response")
modules_np_adaptive <- subset(modules_np, pathway %in% adaptive)
modules_np_adaptive$pathway <- factor(modules_np_adaptive$pathway, levels=c("Treg cell differentiation", "Lymphocyte trafficking", "Mononuclear cell migration", 
                                                                            "T cell receptor signaling", "B cell receptor signaling", "Adaptive immune response"))

adaptive_plot_np <- ggplot(modules_np_adaptive, aes(x = pathway, y = Sample)) +
  geom_point(aes(size = (abs(NES))^3, fill = Expression), shape = 21) +
  geom_text(aes(label = formatC(NES, format='f', digits = 2)), parse = FALSE) +
  scale_fill_manual(values = cen) + # All modules are up-regulated
  labs(x= "Adaptive immunity") + 
  theme(legend.key=element_blank(), legend.position = "none", panel.background = element_rect(colour = "black", fill = NA, linewidth=0.8),
        axis.title.y = element_text(size=14, face="bold", margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12, color = "black"), 
        axis.text.x = element_text(size=12, color = "black", angle=45, hjust=1), 
        axis.ticks.y = element_blank()) + 
  scale_size_continuous(range = c(13,16)) +
  coord_flip() 
adaptive_plot_np

modules_pax <- rbindlist(list(fsgea_pax_0to5_pos_neg_cibersort, fsgea_pax_6to13_pos_neg_cibersort, fsgea_pax_14to20_pos_neg_cibersort, fsgea_pax_adult_pos_neg_cibersort))
modules_pax$pathway <- str_to_sentence(modules_pax$pathway)
modules_pax$pathway[modules_pax$pathway=="Bcr signaling"] <- "B cell receptor signaling"
modules_pax$pathway[modules_pax$pathway=="Nk activity"] <- "NK cell activity"
modules_pax$pathway[modules_pax$pathway=="Tcr signaling"] <- "T cell receptor signaling"
modules_pax$pathway[modules_pax$pathway=="Type i interferon signaling"] <- "Type I interferon signaling"
modules_pax$pathway[modules_pax$pathway=="Type ii interferon signaling"] <- "Type II interferon signaling"
modules_pax$pathway[modules_pax$pathway=="Innate immune cell activation"] <- "Innate immune activation"
modules_pax$pathway[modules_pax$pathway=="Tnf signaling"] <- "TNF signaling"
modules_pax$pathway[modules_pax$pathway=="Tlr signaling"] <- "Toll-like receptor signaling"
modules_pax$pathway[modules_pax$pathway=="Treg differentiation"] <- "Treg cell differentiation"
modules_pax$Sample <- factor(modules_pax$Sample, levels=c("0-5 yr","6-13 yr", "14-20 yr", "Adult"))

if (is.element('Up-regulated', modules_pax$Expression) & is.element('Unchanged', modules_pax$Expression) & is.element('Down-regulated', modules_pax$Expression)) {
  cen = c("lightsteelblue2", "gray80", "#F9AAAE")
} else if ( is.element('Up-regulated', modules_pax$Expression) & is.element('Unchanged', modules_pax$Expression)) {
  cen = c("gray80", "#F9AAAE")
} else if ( is.element('Unchanged', modules_pax$Expression) & is.element('Down-regulated', modules_pax$Expression)) {
  cen = c("lightsteelblue2", "gray80")
} else if ( is.element('Up-regulated', modules_pax$Expression) & is.element('Down-regulated', modules_pax$Expression)) {
  cen = c("lightsteelblue2", "#F9AAAE")
} else {cen = c("gray80")}

modules_pax_innate <- subset(modules_pax, pathway %in% innate)
modules_pax_innate$pathway <- factor(modules_pax_innate$pathway, levels=c("Coagulation", "Complement system", "NK cell activity", "Inflammasomes", "Innate immune activation", 
                                                                          "Myeloid inflammation", "Myeloid activation", "Phagocytosis", "Toll-like receptor signaling", "TNF signaling",
                                                                          "Type II interferon signaling", "Type I interferon signaling", "Interferon response"))

innate_plot_pax <- ggplot(modules_pax_innate, aes(x = pathway, y = Sample)) +
  geom_point(aes(size = (abs(NES))^3, fill = Expression), shape = 21) +
  geom_text(aes(label = formatC(NES, format='f', digits = 2)), parse = FALSE) +
  scale_fill_manual(values = cen) +
  labs(x= "Innate immunity", y = "") + 
  theme(legend.key=element_blank(), legend.position = "none", panel.background = element_rect(colour = "black", fill = NA, linewidth=0.8),
        axis.title.y = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_blank(), axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        plot.title=element_text(size=14, color="black", hjust=0.5, face="bold")) + 
  scale_size_continuous(range = c(13,16)) + ggtitle("Peripheral blood") +
  coord_flip() 
innate_plot_pax

modules_pax_adaptive <- subset(modules_pax, pathway %in% adaptive)
modules_pax_adaptive$pathway <- factor(modules_pax_adaptive$pathway, levels=c("Treg cell differentiation", "Lymphocyte trafficking", "Mononuclear cell migration", 
                                                                              "T cell receptor signaling", "B cell receptor signaling", "Adaptive immune response"))

adaptive_plot_pax <- ggplot(modules_pax_adaptive, aes(x = pathway, y = Sample)) +
  geom_point(aes(size = (abs(NES))^3, fill = Expression), shape = 21) +
  geom_text(aes(label = formatC(NES, format='f', digits = 2)), parse = FALSE) +
  scale_fill_manual(values = cen) +
  labs(x= "Adaptive immunity", y = "") + 
  theme(legend.key=element_blank(), legend.position = "none", panel.background = element_rect(colour = "black", fill = NA, linewidth=0.8),
        axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=12, color = "black", angle=45, hjust=1)) + 
  scale_size_continuous(range = c(13,16)) +
  coord_flip() 
adaptive_plot_pax

module_plot_np <- plot_grid(innate_plot_np, NULL, adaptive_plot_np, labels=NULL, ncol=1, rel_heights=c(4.05,-0.06,2.22), align="v") 
module_plot_pax <- plot_grid(innate_plot_pax, NULL, adaptive_plot_pax, labels=NULL, ncol=1, rel_heights=c(4.05,-0.06,2.22), align="v") 
fig3_modules <- plot_grid(module_plot_np, NULL, module_plot_pax, labels=NULL, ncol=3, rel_widths=c(3.2,-0.05, 1.9), align="vh") 

png(file="Statistical_Analyses/Figures/Figure_3.png", width = 13.5, height = 11, units = 'in', res = 1200)
plot_grid(fig3_cibersort, fig3_modules, labels=c("","c"), ncol=2, rel_widths=c(0.95,1), align="h") 
dev.off()

# Save files as a Source Data file
source_data <- list('Fig3a'=metadata_np_pos_pops, 'Fig3b'=metadata_pax_pos_pops, 'Fig3c_URT_0to5'=fsgea_np_0to5_pos_neg_nocibersort,
                    'Fig3c_URT_6to13'=fsgea_np_6to13_pos_neg_nocibersort, 'Fig3c_URT_14to20'=fsgea_np_14to20_pos_neg_nocibersort, 
                    'Fig3c_BLD_0to5'=fsgea_pax_0to5_pos_neg_cibersort, 'Fig3c_BLD_6to13'=fsgea_pax_6to13_pos_neg_cibersort,
                    'Fig3c_BLD_14to20'=fsgea_pax_14to20_pos_neg_cibersort, 'Fig3c_BLD_adult'=fsgea_pax_adult_pos_neg_cibersort)
openxlsx::write.xlsx(source_data, file="Statistical_Analyses/Source_Data/Figure_3.xlsx")