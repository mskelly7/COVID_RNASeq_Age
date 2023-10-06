# BRAVE Kids RNA Sequencing Analysis
# Aditya Mohan (MD/PhD candidate)  / Matthew Kelly, MD, MPH 
# Figure 6: Correlation between URT immune responses in acute infection and neutralizing antibody responses
# Last update: Oct. 6, 2023

remove(list=ls())
setwd("____________________________________") 
options(rstudio.help.showDataPreview = FALSE)
set.seed(1234)
version

library(readr)
library(circlize)
library(data.table)
library(Hmisc)
library(readxl)
library(writexl)
library(readr)
library(dplyr)
library(ggplot2)
library(gridExtra) 
library(cowplot) 
library(forcats) 
library(tibble)
library(reshape2)
library(phyloseq)
library(plyr)
library(ggtext)
library(GSVA)
library(tabletools)

# Load the required datasets
phy.rnaseq.np <- readRDS("phy.rnaseq.np.rds")
phy.rnaseq.np.neut <- subset_samples(phy.rnaseq.np, corona=="Positive" & !is.na(neut_ID50_2mo))
nsamples(phy.rnaseq.np.neut)
genes_np_neut <- data.frame(tax_table(phy.rnaseq.np.neut))
counts_np_neut <- data.frame(otu_table(phy.rnaseq.np.neut))
metadata_np_neut <- data.frame(sample_data(phy.rnaseq.np.neut))
neuts_only <- metadata_np_neut[,c("alias_sequencing_id","neut_ID50_2mo")]

# Evaluate for correlations between imputed immune cell proportions and neutralizing antibody responses at 2 months

metadata_np_neut <- metadata_np_neut[,c("neut_ID50_2mo","B_cells","T_cells_CD4","T_cells_CD8","Dendritic_cells","Eosinophils","Mast_cells","Mono_Macrophages","Neutrophils","NK_cells","Plasma_cells")]
corr_metadata_np_neut <- rcorr(as.matrix(metadata_np_neut),type="pearson")
diag(corr_metadata_np_neut$P) <- 0
r_mat_neut <- corr_metadata_np_neut$r
r_mat_neut <- r_mat_neut[1,2:11]
p_mat_neut <- corr_metadata_np_neut$P
p_mat_neut <- p_mat_neut[1,2:11]
corr_cells <- data.frame(t(rbind(r_mat_neut, p_mat_neut)))
# No significant correlations between imputed immune cell proportions and 2-month antibody response

# Evaluate for correlations between immune module expression and neutralizing antibody responses at 2 months

genes_np_neut <- data.frame(tax_table(phy.rnaseq.np.neut))
counts_np_neut <- data.frame(otu_table(phy.rnaseq.np.neut))
data_np_neut <- merge(genes_np_neut, counts_np_neut, by="row.names")
data_np_neut <- data_np_neut[,-c(1,3,4)] # Drop Row.names, Chromosome, Length
nrow(data_np_neut)
# Sum counts for genes that have multiple EnsemblIDs
data_np_neut <- ddply(data_np_neut, "Gene", numcolwise(sum))
nrow(data_np_neut)
row.names(data_np_neut) <- data_np_neut$Gene
data_np_neut$Gene <- NULL
data_np_neut <- as.matrix(data_np_neut)
data_np_neut_log2 <- log2(data_np_neut)

modules <- read_excel("Statistical_Analyses/modules_61.xlsx", sheet="modules_61", na="NA")
modules_names <- colnames(modules)
modules <- as.list(as.data.frame(modules))
modules <- lapply(modules, function(x) x[x!=""])
ssgsea_np <- gsva(data_np_neut_log2, modules, method="ssgsea")
ssgsea_np <- as.data.frame(t(ssgsea_np))
colnames(ssgsea_np) <- paste(colnames(ssgsea_np),"NSB",sep=".")
ssgsea_np_neut <- merge(ssgsea_np, neuts_only, by="row.names")
row.names(ssgsea_np_neut) <- ssgsea_np_neut$alias_sequencing_id
ssgsea_np_neut <- ssgsea_np_neut[,-c(1,63)]
# Correlations with 2-month antibody responses
corr_modules <- rcorr(as.matrix(ssgsea_np_neut),type="pearson")
corr_modules_padj <- rcorr_padjust(corr_modules, method = "BH")
diag(corr_modules_padj$P) <- 0
r_mat_neut <- corr_modules_padj$r
r_mat_neut <- r_mat_neut[62,1:61]
p_mat_neut <- corr_modules_padj$P
p_mat_neut <- p_mat_neut[62,1:61]
corr_modules_padj <- data.frame(t(rbind(r_mat_neut, p_mat_neut)))
write.csv(corr_modules_padj, "Statistical_Analyses/6_Neutralizing_Antibodies/NP_vs_Neuts_correlations_all.csv")
corr_modules_padj_sig <- subset(corr_modules_padj, p_mat_neut<0.05)
corr_modules_padj_sig

# Create correlation plots for modules of interest

names(ssgsea_np_neut) <- sub(" ", ".", names(ssgsea_np_neut))
names(ssgsea_np_neut) <- sub(" ", ".", names(ssgsea_np_neut))
names(ssgsea_np_neut) <- sub(" ", ".", names(ssgsea_np_neut))
names(ssgsea_np_neut) <- sub(" ", ".", names(ssgsea_np_neut))
names(ssgsea_np_neut) <- sub("-", ".", names(ssgsea_np_neut))
names(ssgsea_np_neut) <- sub("-", ".", names(ssgsea_np_neut))
names(ssgsea_np_neut) <- sub("-", ".", names(ssgsea_np_neut))

#Interferon Response Genes.NSB 0.5671815 0.002871133
summary(ssgsea_np_neut$Interferon.Response.Genes.NSB)
Interferon.Response.Genes.NSB.plot <- ggplot(data = ssgsea_np_neut, aes(x = neut_ID50_2mo, y = Interferon.Response.Genes.NSB)) +
  geom_point() +
  geom_smooth(method = "lm", color="red", linewidth=0.5, se = FALSE) +
  xlim(0,1500) + scale_y_continuous(breaks=c(0.40,0.45,0.50,0.55,0.60,0.65,0.70), limits=c(0.40,0.70), labels = scales::number_format(accuracy = 0.01, decimal.mark = '.')) + 
  labs(title="Interferon response genes") + ylab("NES") + xlab("") +
  annotate(geom="text", label="ρ=0.56,      =0.02", x=1200, y=0.426, size=4.5) +
  annotate(geom="text", label=deparse(bquote(p[adj])), x=1232, y=0.4225, size=4.5, parse=TRUE) +
  theme_bw() +
  theme(plot.title = ggtext::element_markdown(size=14, hjust=0.5, face="bold"), 
        axis.title.x = element_text(size = 14, vjust=-0.3), axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) 
Interferon.Response.Genes.NSB.plot

#DNA.Sensing.NSB 0.5599415 0.003347087
summary(ssgsea_np_neut$DNA.Sensing.NSB)
DNA.Sensing.NSB.plot <- ggplot(data = ssgsea_np_neut, aes(x = neut_ID50_2mo, y = DNA.Sensing.NSB)) +
  geom_point() +
  geom_smooth(method = "lm", color="red", linewidth=0.5, se = FALSE) +
  xlim(0,1500) + scale_y_continuous(breaks=c(0.30,0.35,0.40,0.45,0.50,0.55), limits=c(0.30,0.55), labels = scales::number_format(accuracy = 0.01, decimal.mark = '.')) + 
  labs(title="DNA sensing") + ylab("") + xlab("") +
  annotate(geom="text", label="ρ=0.56,      =0.003", x=1200, y=0.325, size=4.5) +
  annotate(geom="text", label=deparse(bquote(p[adj])), x=1212, y=0.3221, size=4.5, parse=TRUE) +
  theme_bw() +
  theme(plot.title = ggtext::element_markdown(size=14, hjust=0.5, face="bold"), 
        axis.title.x = element_text(size = 14, margin = margin(t = 0, r = 0, b = 10, l = 0)), axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) 
DNA.Sensing.NSB.plot

#RNA Sensing.NSB 0.5279824 0.006292867
summary(ssgsea_np_neut$RNA.Sensing.NSB)
RNA.Sensing.NSB.plot <- ggplot(data = ssgsea_np_neut, aes(x = neut_ID50_2mo, y = RNA.Sensing.NSB)) +
  geom_point() +
  geom_smooth(method = "lm", color="red", linewidth=0.5, se = FALSE) +
  xlim(0,1500) + scale_y_continuous(breaks=c(0.30,0.35,0.40,0.45,0.50,0.55), limits=c(0.30,0.55), labels = scales::number_format(accuracy = 0.01, decimal.mark = '.')) + 
  labs(title="RNA sensing") + ylab("") + xlab("") +
  annotate(geom="text", label="ρ=0.53,      =0.006", x=1200, y=0.325, size=4.5) +
  annotate(geom="text", label=deparse(bquote(p[adj])), x=1212, y=0.3221, size=4.5, parse=TRUE) +
  theme_bw() +
  theme(plot.title = ggtext::element_markdown(size=14, hjust=0.5, face="bold"), 
        axis.title.x = element_text(size = 14, vjust=-0.3), axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) 
RNA.Sensing.NSB.plot

#MHC.Class.I.Antigen.Presentation.NSB 0.4499562 0.023487391
summary(ssgsea_np_neut$MHC.Class.I.Antigen.Presentation.NSB)
MHC.Class.I.Antigen.Presentation.NSB.plot <- ggplot(data = ssgsea_np_neut, aes(x = neut_ID50_2mo, y = MHC.Class.I.Antigen.Presentation.NSB)) +
  geom_point() +
  geom_smooth(method = "lm", color="red", linewidth=0.5, se = FALSE) +
  xlim(0,1500) + scale_y_continuous(breaks=c(0.35,0.40,0.45,0.50,0.55,0.60), limits=c(0.35,0.60), labels = scales::number_format(accuracy = 0.01, decimal.mark = '.')) + 
  labs(title="MHC class I antigen presentation") + ylab("NES") + xlab("") +
  annotate(geom="text", label="ρ=0.45,      =0.02", x=1200, y=0.375, size=4.5) +
  annotate(geom="text", label=deparse(bquote(p[adj])), x=1232, y=0.3721, size=4.5, parse=TRUE) +
  theme_bw() +
  theme(plot.title = ggtext::element_markdown(size=14, hjust=0.5, face="bold"), 
        axis.title.x = element_text(size = 14, vjust=-0.3), axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) 
MHC.Class.I.Antigen.Presentation.NSB.plot

#T.cell.Costimulation.NSB 0.4669279 0.018041577
summary(ssgsea_np_neut$T.cell.Costimulation.NSB)
T.cell.Costimulation.NSB.plot <- ggplot(data = ssgsea_np_neut, aes(x = neut_ID50_2mo, y = T.cell.Costimulation.NSB)) +
  geom_point() +
  geom_smooth(method = "lm", color="red", linewidth=0.5, se = FALSE) +
  xlim(0,1500) + scale_y_continuous(breaks=c(-0.40,-0.30,-0.20,-0.10,0.00,0.10), limits=c(-0.40,0.10), labels = scales::number_format(accuracy = 0.01, decimal.mark = '.')) + 
  labs(title="T cell costimulation") + ylab("") + xlab("Serum neutralization activity (ID50)") +
  annotate(geom="text", label="ρ=0.47,      =0.02", x=1200, y=-0.35, size=4.5) +
  annotate(geom="text", label=deparse(bquote(p[adj])), x=1232, y=-0.357, size=4.5, parse=TRUE) +
  theme_bw() +
  theme(plot.title = ggtext::element_markdown(size=14, hjust=0.5, face="bold"), 
        axis.title.x = element_text(size = 14, vjust=-0.3), axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) 
T.cell.Costimulation.NSB.plot

#Immune.Exhaustion.NSB  0.4669279 0.018041577
summary(ssgsea_np_neut$Immune.Exhaustion.NSB)
Immune.Exhaustion.NSB.plot <- ggplot(data = ssgsea_np_neut, aes(x = neut_ID50_2mo, y = Immune.Exhaustion.NSB)) +
  geom_point() +
  geom_smooth(method = "lm", color="red", linewidth=0.5, se = FALSE) +
  xlim(0,1500) + scale_y_continuous(breaks=c(-0.35,-0.25,-0.15,-0.05,0.05,0.15), limits=c(-0.35,0.15), labels = scales::number_format(accuracy = 0.01, decimal.mark = '.')) + 
  labs(title="Immune exhaustion") + ylab("") + xlab("") + 
  annotate(geom="text", label="ρ=0.47,      =0.02", x=1200, y=-0.30, size=4.5) +
  annotate(geom="text", label=deparse(bquote(p[adj])), x=1227, y=-0.307, size=4.5, parse=TRUE) +
  theme_bw() +
  theme(plot.title = ggtext::element_markdown(size=14, hjust=0.5, face="bold"), 
        axis.title.x = element_text(size = 14, vjust=-0.3), axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) 
Immune.Exhaustion.NSB.plot

# NO TITLE

#png(file="Statistical_Analyses/Figures/Figure_6.png", width = 16, height = 10.5, units = 'in', res = 1200)
#plot_grid(Interferon.Response.Genes.NSB.plot, RNA.Sensing.NSB.plot, DNA.Sensing.NSB.plot, MHC.Class.I.Antigen.Presentation.NSB.plot, T.cell.Costimulation.NSB.plot, Immune.Exhaustion.NSB.plot,
#          labels=NULL, ncol=3, align="vh") 
#dev.off()

# ADDITION OF TITLE

fig6_plot <- plot_grid(Interferon.Response.Genes.NSB.plot, RNA.Sensing.NSB.plot, DNA.Sensing.NSB.plot, MHC.Class.I.Antigen.Presentation.NSB.plot, T.cell.Costimulation.NSB.plot, Immune.Exhaustion.NSB.plot,
                                 labels=NULL, ncol=3, align="vh") 
fig6_title <- ggdraw() + draw_label("Figure 6. Correlations between upper respiratory immune modules and convalescent neutralizing antibody responses", 
                                    size=13.5, fontface='bold')

png(file="Statistical_Analyses/Figures/Figure_6.png", width = 16, height = 10.5, units = 'in', res = 1200)
plot_grid(fig6_title, NULL, fig6_plot, labels=NULL, ncol=1, rel_heights=c(0.03,0.02,0.97), align="v") 
dev.off()

#Interferon.Response.Genes.NSB         0.5671815 0.002871133
#RNA.Sensing.NSB                       0.5279824 0.006292867
#DNA.Sensing.NSB                       0.5599415 0.003347087
#MHC.Class.I.Antigen.Presentation.NSB  0.4499562 0.023487391
#T.cell.Costimulation.NSB              0.4605413 0.019904315
#Immune.Exhaustion.NSB                 0.4669279 0.018041577