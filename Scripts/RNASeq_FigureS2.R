# BRAVE Kids RNA Sequencing Analysis
# Aditya Mohan (MD/PhD candidate)  / Matthew Kelly, MD, MPH 
# Figure S2 - CIBERSORT plots by COVID status
# Last update: July 27, 2023

remove(list=ls())
setwd("______________________") 
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
library(tibble)
library(msigdbr)
library(reshape2)
library(data.table)

# Upload files with CIBERSORTx data
phy.rnaseq.np <- readRDS("phy.rnaseq.np.rds")
metadata_np <- data.frame(sample_data(phy.rnaseq.np))
metadata_np_pops <- metadata_np[,c("corona","B_cells","Plasma_cells", "T_cells_CD8","T_cells_CD4","T_cells_gamma_delta","NK_cells",
                                           "Mono_Macrophages","Dendritic_cells","Mast_cells","Eosinophils","Neutrophils")]
metadata_np_pops <- gather(metadata_np_pops, cell_type, value, B_cells:Neutrophils, factor_key=TRUE)
metadata_np_pops$cell_type2[metadata_np_pops$cell_type=="B_cells"] <- "B cells"
metadata_np_pops$cell_type2[metadata_np_pops$cell_type=="Plasma_cells"] <- "Plasma cells"
metadata_np_pops$cell_type2[metadata_np_pops$cell_type=="T_cells_CD8"] <- "CD4+ T cells"
metadata_np_pops$cell_type2[metadata_np_pops$cell_type=="T_cells_CD4"] <- "CD8+ T cells"
metadata_np_pops$cell_type2[metadata_np_pops$cell_type=="T_cells_gamma_delta"] <- "Gamma delta T cells"
metadata_np_pops$cell_type2[metadata_np_pops$cell_type=="NK_cells"] <- "NK cells"
metadata_np_pops$cell_type2[metadata_np_pops$cell_type=="Mono_Macrophages"] <- "Mono/macrophages"
metadata_np_pops$cell_type2[metadata_np_pops$cell_type=="Dendritic_cells"] <- "Dendritic cells"
metadata_np_pops$cell_type2[metadata_np_pops$cell_type=="Mast_cells"] <- "Mast cells"
metadata_np_pops$cell_type2[metadata_np_pops$cell_type=="Eosinophils"] <- "Eosinophils"
metadata_np_pops$cell_type2[metadata_np_pops$cell_type=="Neutrophils"] <- "Neutrophils"
tapply(metadata_np_pops$value, metadata_np_pops$cell_type2, summary) 
# Exclude cell types identified in <25% of samples from figure
metadata_np_pops <- subset(metadata_np_pops, cell_type2!="Gamma delta T cells")

phy.rnaseq.pax <- readRDS("phy.rnaseq.pax.rds")
metadata_pax <- data.frame(sample_data(phy.rnaseq.pax))
metadata_pax_pops <- metadata_pax[,c("corona","B_cells","Plasma_cells", "T_cells_CD8","T_cells_CD4","T_cells_gamma_delta","NK_cells",
                                             "Mono_Macrophages","Dendritic_cells","Mast_cells","Eosinophils","Neutrophils")]
metadata_pax_pops <- gather(metadata_pax_pops, cell_type, value, B_cells:Neutrophils, factor_key=TRUE)
metadata_pax_pops$cell_type2[metadata_pax_pops$cell_type=="B_cells"] <- "B cells"
metadata_pax_pops$cell_type2[metadata_pax_pops$cell_type=="Plasma_cells"] <- "Plasma cells"
metadata_pax_pops$cell_type2[metadata_pax_pops$cell_type=="T_cells_CD8"] <- "CD4+ T cells"
metadata_pax_pops$cell_type2[metadata_pax_pops$cell_type=="T_cells_CD4"] <- "CD8+ T cells"
metadata_pax_pops$cell_type2[metadata_pax_pops$cell_type=="T_cells_gamma_delta"] <- "Gamma delta T cells"
metadata_pax_pops$cell_type2[metadata_pax_pops$cell_type=="NK_cells"] <- "NK cells"
metadata_pax_pops$cell_type2[metadata_pax_pops$cell_type=="Mono_Macrophages"] <- "Mono/macrophages"
metadata_pax_pops$cell_type2[metadata_pax_pops$cell_type=="Dendritic_cells"] <- "Dendritic cells"
metadata_pax_pops$cell_type2[metadata_pax_pops$cell_type=="Mast_cells"] <- "Mast cells"
metadata_pax_pops$cell_type2[metadata_pax_pops$cell_type=="Eosinophils"] <- "Eosinophils"
metadata_pax_pops$cell_type2[metadata_pax_pops$cell_type=="Neutrophils"] <- "Neutrophils"
tapply(metadata_pax_pops$value, metadata_pax_pops$cell_type2, summary) 
# Exclude cell types identified in <25% of samples from figure
metadata_pax_pops <- subset(metadata_pax_pops, cell_type2!="Dendritic cells" & cell_type2!="Eosinophils" & cell_type2!="Gamma delta T cells")

# UPPER RESPIRATORY SAMPLES (COVID-NEGATIVE VS. COVID-POSITIVE)
# Adjust comparisons for multiple testing using BH (n=11)
# B cells: 0.52
# Plasma cells: 0.29
# CD8+ T cells: >0.99
# CD4+ T cells: >0.99
# Gamma Delta T cells: >0.99
# NK cells: 0.33
# Mono/macrophages: 0.03
# Dendritic cells: 0.02
# Mast cells: >0.99
# Eosinophils: >0.99
# Neutrophils: >0.99 

metadata_np_pops$corona2[metadata_np_pops$corona=="Negative"] <- "Uninfected"
metadata_np_pops$corona2[metadata_np_pops$corona=="Positive"] <- "Infected"
cibersort_np <- ggplot(metadata_np_pops, aes(x=cell_type2, y=value, fill=corona2)) + geom_boxplot() + ylab("Imputed immune cell proportions") + ylim(0,0.61) +
  theme(legend.title = element_blank(), legend.position = "right", legend.text = element_text(size=9), legend.box.spacing = unit(3, "pt"), legend.key=element_rect(fill="white"),
        panel.background = element_rect(colour = "black", fill = NA, linewidth=0.5), panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey70"),
        panel.grid.minor.y = element_line(linewidth = 0.2, colour = "grey70"),
        axis.title.y = element_text(size=10, face="bold", margin = ggplot2::margin(t = 0, r = 5, b = 0, l = 0)),
        axis.title.x = element_blank(), axis.text.y = element_text(size=9), axis.text.x = element_text(size=9, color = "black", angle=35, hjust=1), 
        plot.title = element_text(size=11, hjust = 0.5, face="bold")) + scale_fill_manual(values = c("cornflowerblue","indianred1")) + 
  ggtitle("Upper respiratory") + 
  annotate("text", x=7, y=0.61, size=5, label= "*") + # Mono/macrophages (age-adjusted beta regression, p=0.03)
  annotate("text", x=4, y=0.61, size=5, label= "*") # Dendritic cells (age-adjusted beta regression, p=0.02)

# PERIPHERAL BLOOD SAMPLES
# Adjust comparisons for multiple testing using BH (n=9)
# B cells: 0.94 
# Plasma cells: 0.03
# CD8+ T cells: >0.99
# CD4+ T cells: 0.15
# NK cells: 0.96
# Mono/macrophages: >0.99
# Dendritic cells: >0.99
# Mast cells: >0.99
# Neutrophils: 0.28

metadata_pax_pops$corona2[metadata_pax_pops$corona=="Negative"] <- "Uninfected"
metadata_pax_pops$corona2[metadata_pax_pops$corona=="Positive"] <- "Infected"
cibersort_pax <- ggplot(metadata_pax_pops, aes(x=cell_type2, y=value, fill=corona2)) + geom_boxplot() + ylab("Imputed immune cell proportions") + ylim(0,0.51) +
  theme(legend.title = element_blank(), legend.position = "right", legend.text = element_text(size=9), legend.box.spacing = unit(3, "pt"), legend.key=element_rect(fill="white"),
        panel.background = element_rect(colour = "black", fill = NA, linewidth=0.5), panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey70"),
        panel.grid.minor.y = element_line(linewidth = 0.2, colour = "grey70"),
        axis.title.y = element_text(size=10, face="bold", margin = ggplot2::margin(t = 0, r = 5, b = 0, l = 0)),
        axis.title.x = element_blank(), axis.text.y = element_text(size=9), axis.text.x = element_text(size=9, color = "black", angle=35, hjust=1), 
        plot.title = element_text(size=11, hjust = 0.5, face="bold")) + scale_fill_manual(values = c("cornflowerblue","indianred1")) + 
  ggtitle("Peripheral blood") + 
  annotate("text", x=8, y=0.51, size=5, label= "*") # Plasma cells (age-adjusted beta regression, p=0.03)

png(file="Figures/Figure_S2.png", width = 6, height = 8, units = 'in', res = 1200)
plot_grid(cibersort_np, cibersort_pax, ncol=1) 
dev.off()