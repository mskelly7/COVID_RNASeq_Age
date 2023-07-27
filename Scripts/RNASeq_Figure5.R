# BRAVE Kids RNA Sequencing Analysis
# Aditya Mohan (MD/PhD candidate)  / Matthew Kelly, MD, MPH 
# Figure 5: heatmap with NP vs. PAX module correlations 
# Last update: July 27, 2023

remove(list=ls())
setwd("________________________") 
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
library(GSVA)

# Load the required datasets
phy.rnaseq.np <- readRDS("phy.rnaseq.np.rds")
phy.rnaseq.np.intersect <- subset_samples(phy.rnaseq.np, corona=="Positive")
nsamples(phy.rnaseq.np.intersect)
metadata_np_intersect <- data.frame(sample_data(phy.rnaseq.np.intersect))
np_ids <- metadata_np_intersect$alias_study_id

phy.rnaseq.pax <- readRDS("phy.rnaseq.pax.rds")
phy.rnaseq.pax.intersect <- subset_samples(phy.rnaseq.pax, corona=="Positive")
nsamples(phy.rnaseq.pax.intersect)
metadata_pax_intersect <- data.frame(sample_data(phy.rnaseq.pax.intersect))
pax_ids <- metadata_pax_intersect$alias_study_id
final_ids <- intersect(np_ids, pax_ids)

phy.rnaseq.np <- readRDS("phy.rnaseq.np.rds")
phy.rnaseq.np.Pos <- subset_samples(phy.rnaseq.np, corona=="Positive" & alias_study_id %in% final_ids)
nsamples(phy.rnaseq.np.Pos)
genes_np_Pos <- data.frame(tax_table(phy.rnaseq.np.Pos))
counts_np_Pos <- data.frame(otu_table(phy.rnaseq.np.Pos))
data_np_Pos <- merge(genes_np_Pos, counts_np_Pos, by="row.names")
data_np_Pos <- data_np_Pos[,-c(1,3,4)] # Drop Row.names, Chromosome, Length
nrow(data_np_Pos)
# Sum counts for genes that have multiple EnsemblIDs
data_np_Pos <- ddply(data_np_Pos, "Gene", numcolwise(sum))
nrow(data_np_Pos)
row.names(data_np_Pos) <- data_np_Pos$Gene
data_np_Pos$Gene <- NULL
data_np_Pos <- as.matrix(data_np_Pos)
data_np_Pos_log2 <- log2(data_np_Pos)

phy.rnaseq.pax <- readRDS("G:/My Drive/Research/BRAVE_Kids/RNA_Sequencing/phy.rnaseq.pax.rds")
phy.rnaseq.pax.Pos <- subset_samples(phy.rnaseq.pax, corona=="Positive" & alias_study_id %in% final_ids)
nsamples(phy.rnaseq.pax.Pos)
genes_pax_Pos <- data.frame(tax_table(phy.rnaseq.pax.Pos))
counts_pax_Pos <- data.frame(otu_table(phy.rnaseq.pax.Pos))
names(counts_pax_Pos) <- gsub(x = names(counts_pax_Pos), pattern = "X025CBE.M0.PAX", replacement = "025CBE.M0.PAX")
names(counts_pax_Pos) <- gsub(x = names(counts_pax_Pos), pattern = "X0719FA.M0.PAX", replacement = "0719FA.M0.PAX")
names(counts_pax_Pos) <- gsub(x = names(counts_pax_Pos), pattern = "X0BF51C.M0.PAX", replacement = "0BF51C.M0.PAX")
names(counts_pax_Pos) <- gsub(x = names(counts_pax_Pos), pattern = "X0BF51C.M6.PAX", replacement = "0BF51C.M6.PAX")
names(counts_pax_Pos) <- gsub(x = names(counts_pax_Pos), pattern = "X0E1ECA.M0.PAX", replacement = "0E1ECA.M0.PAX")
names(counts_pax_Pos) <- gsub(x = names(counts_pax_Pos), pattern = "X18EF5E.M0.PAX", replacement = "18EF5E.M0.PAX")
names(counts_pax_Pos) <- gsub(x = names(counts_pax_Pos), pattern = "X1AC6B8.M0.PAX", replacement = "1AC6B8.M0.PAX")
names(counts_pax_Pos) <- gsub(x = names(counts_pax_Pos), pattern = "X1E06CA.M0.PAX", replacement = "1E06CA.M0.PAX")
names(counts_pax_Pos) <- gsub(x = names(counts_pax_Pos), pattern = "X234BBF.M0.PAX", replacement = "234BBF.M0.PAX")
names(counts_pax_Pos) <- gsub(x = names(counts_pax_Pos), pattern = "X322088.M0.PAX", replacement = "322088.M0.PAX")
names(counts_pax_Pos) <- gsub(x = names(counts_pax_Pos), pattern = "X3501A9.M0.PAX", replacement = "3501A9.M0.PAX")
names(counts_pax_Pos) <- gsub(x = names(counts_pax_Pos), pattern = "X3B441F.M0.PAX", replacement = "3B441F.M0.PAX")
names(counts_pax_Pos) <- gsub(x = names(counts_pax_Pos), pattern = "X496933.M0.PAX", replacement = "496933.M0.PAX")
names(counts_pax_Pos) <- gsub(x = names(counts_pax_Pos), pattern = "X58A568.M0.PAX", replacement = "58A568.M0.PAX")
names(counts_pax_Pos) <- gsub(x = names(counts_pax_Pos), pattern = "X6AA2F7.M0.PAX", replacement = "6AA2F7.M0.PAX")
names(counts_pax_Pos) <- gsub(x = names(counts_pax_Pos), pattern = "X71549F.M0.PAX", replacement = "71549F.M0.PAX")
names(counts_pax_Pos) <- gsub(x = names(counts_pax_Pos), pattern = "X72C971.M0.PAX", replacement = "72C971.M0.PAX")
names(counts_pax_Pos) <- gsub(x = names(counts_pax_Pos), pattern = "X7D880A.M0.PAX", replacement = "7D880A.M0.PAX")
names(counts_pax_Pos) <- gsub(x = names(counts_pax_Pos), pattern = "X826E5F.M0.PAX", replacement = "826E5F.M0.PAX")
names(counts_pax_Pos) <- gsub(x = names(counts_pax_Pos), pattern = "X830937.M0.PAX", replacement = "830937.M0.PAX")
names(counts_pax_Pos) <- gsub(x = names(counts_pax_Pos), pattern = "X8867DF.M0.PAX", replacement = "8867DF.M0.PAX")
names(counts_pax_Pos) <- gsub(x = names(counts_pax_Pos), pattern = "X88732B.M0.PAX", replacement = "88732B.M0.PAX")
names(counts_pax_Pos) <- gsub(x = names(counts_pax_Pos), pattern = "X88B192.M0.PAX", replacement = "88B192.M0.PAX")
data_pax_Pos <- merge(genes_pax_Pos, counts_pax_Pos, by="row.names")
data_pax_Pos <- data_pax_Pos[,-c(1,3,4)] # Drop Row.names, Chromosome, Length
nrow(data_pax_Pos)
# Sum counts for genes that have multiple EnsemblIDs
data_pax_Pos <- ddply(data_pax_Pos, "Gene", numcolwise(sum))
nrow(data_pax_Pos)
row.names(data_pax_Pos) <- data_pax_Pos$Gene
data_pax_Pos$Gene <- NULL
data_pax_Pos <- as.matrix(data_pax_Pos)
data_pax_Pos_log2 <- log2(data_pax_Pos)

# Upload module information 
modules <- read_excel("modules_61.xlsx", sheet="modules_61", na="NA")
pathways_to_keep <- c("Innate Immune Cell Activation","Interferon Response","Type I Interferon Signaling","Type II Interferon Signaling","Type III Interferon Signaling", 
                      "TNF Signaling","Chemokine Signaling","TLR Signaling","NLR Signaling","RNA Sensing","Phagocytosis","Myeloid Activation",
                      "Myeloid Inflammation","Inflammasomes","Complement System","Coagulation","Leukotriene and Prostaglandin Inflammation",
                      "Adaptive Immune Response","MHC Class I Antigen Presentation","MHC Class II Antigen Presentation","Mononuclear Cell Migration","Lymphocyte Trafficking",
                      "BCR Signaling","NF-kappaB Signaling","TCR Signaling","NK Activity","JAK-STAT Signaling","Immune Memory","Immune Exhaustion","Treg Differentiation")
modules <- modules[,names(modules) %in% pathways_to_keep]
colnames(modules)[which(names(modules) == "Innate Immune Cell Activation")] <- "Innate immune cell activation"
colnames(modules)[which(names(modules) == "Interferon Response")] <- "Interferon response"
colnames(modules)[which(names(modules) == "Type I Interferon Signaling")] <- "Type I interferon signaling"
colnames(modules)[which(names(modules) == "Type II Interferon Signaling")] <- "Type II interferon signaling"
colnames(modules)[which(names(modules) == "Type III Interferon Signaling")] <- "Type III interferon signaling"
colnames(modules)[which(names(modules) == "TNF Signaling")] <- "TNF signaling"
colnames(modules)[which(names(modules) == "Chemokine Signaling")] <- "Chemokine signaling"
colnames(modules)[which(names(modules) == "TLR Signaling")] <- "Toll-like receptor signaling"
colnames(modules)[which(names(modules) == "NLR Signaling")] <- "Nod-like receptor signaling"
colnames(modules)[which(names(modules) == "RNA Sensing")] <- "RNA sensing"
colnames(modules)[which(names(modules) == "Myeloid Activation")] <- "Myeloid activation"
colnames(modules)[which(names(modules) == "Myeloid Inflammation")] <- "Myeloid inflammation"
colnames(modules)[which(names(modules) == "Complement System")] <- "Complement system"
colnames(modules)[which(names(modules) == "Leukotriene and Prostaglandin Inflammation")] <- "Leukotriene/prostaglandin inflammation"
colnames(modules)[which(names(modules) == "Adaptive Immune Response")] <- "Adaptive immune response"
colnames(modules)[which(names(modules) == "MHC Class I Antigen Presentation")] <- "MHC class I presentation"
colnames(modules)[which(names(modules) == "MHC Class II Antigen Presentation")] <- "MHC class II presentation"
colnames(modules)[which(names(modules) == "Mononuclear Cell Migration")] <- "Mononuclear cell migration"
colnames(modules)[which(names(modules) == "Lymphocyte Trafficking")] <- "Lymphocyte trafficking"
colnames(modules)[which(names(modules) == "BCR Signaling")] <- "BCR signaling"
colnames(modules)[which(names(modules) == "NF-kappaB Signaling")] <- "NF-kappaB signaling"
colnames(modules)[which(names(modules) == "TCR Signaling")] <- "TCR signaling"
colnames(modules)[which(names(modules) == "NK Activity")] <- "NK cell activity"
colnames(modules)[which(names(modules) == "JAK-STAT Signaling")] <- "JAK-STAT signaling"
colnames(modules)[which(names(modules) == "Immune Memory")] <- "Immune memory"
colnames(modules)[which(names(modules) == "Immune Exhaustion")] <- "Immune exhaustion"
colnames(modules)[which(names(modules) == "Treg Differentiation")] <- "Treg differentiation"
modules_names <- colnames(modules)
modules <- as.list(as.data.frame(modules))
modules <- lapply(modules, function(x) x[x!=""])

# flattenCorrMatrix function (http://www.sthda.com/english/wiki/correlation-matrix-formatting-and-visualization)
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

# ALL AGES

ssgsea_np <- gsva(data_np_Pos_log2, modules, method="ssgsea")
ssgsea_np <- as.data.frame(t(ssgsea_np))
pathway_order <- c("Innate immune cell activation","Interferon response","Type I interferon signaling","Type II interferon signaling","Type III interferon signaling","TNF signaling",
                   "Toll-like receptor signaling","Nod-like receptor signaling","Chemokine signaling","RNA sensing","Phagocytosis","Myeloid inflammation","Myeloid activation","Inflammasomes",
                   "Complement system","Coagulation","Leukotriene/prostaglandin inflammation","Adaptive immune response","MHC class I presentation","MHC class II presentation",
                   "TCR signaling","BCR signaling","NF-kappaB signaling","JAK-STAT signaling","NK cell activity","Mononuclear cell migration","Lymphocyte trafficking",
                   "Immune memory","Immune exhaustion","Treg differentiation")
ssgsea_np <- ssgsea_np[,pathway_order]
colnames(ssgsea_np) <- paste(colnames(ssgsea_np),"NSB",sep=".")
row.names(ssgsea_np) <- gsub(".NSB", "", row.names(ssgsea_np))

ssgsea_pax <- gsva(data_pax_Pos_log2, modules, method="ssgsea")
ssgsea_pax  <- (ssgsea_pax - rowMeans(ssgsea_pax))/(rowSds(as.matrix(ssgsea_pax)))[row(ssgsea_pax)]
ssgsea_pax <- as.data.frame(t(ssgsea_pax))
ssgsea_pax <- ssgsea_pax[,pathway_order]
colnames(ssgsea_pax) <- paste(colnames(ssgsea_pax),"PAX",sep=".")
row.names(ssgsea_pax) <- gsub(".PAX", "", row.names(ssgsea_pax))

ssgsea_all <- merge(ssgsea_np, ssgsea_pax, by="row.names")
row.names(ssgsea_all) <- ssgsea_all[,1]
ssgsea_all$Row.names <- NULL
corr_mat <- rcorr(as.matrix(ssgsea_all),type="pearson")
diag(corr_mat$P) <- 0
r_mat <- corr_mat$r
r_mat <- r_mat[1:30,31:60]
row.names(r_mat) <- gsub(".NSB","",row.names(r_mat))
colnames(r_mat) <- gsub(".PAX","",row.names(r_mat))
p_mat <- corr_mat$P
p_mat <- p_mat[1:30,31:60]
row.names(p_mat) <- gsub(".NSB","",row.names(p_mat))
colnames(p_mat) <- gsub(".PAX","",row.names(p_mat))
corr_tab <- flattenCorrMatrix(corr_mat$r, corr_mat$P)
corr_tab <- dplyr::filter(corr_tab, grepl('.NSB', row))
corr_tab <- dplyr::filter(corr_tab, grepl('.PAX', column))
write.csv(corr_tab, "5_NP_vs_PAX_Correlations/NP_vs_PAX_correlations_all.csv")

png(file="Figures/Figure_5.png", width = 7.5, height = 5, units = 'in', res = 1200)
corrplot(r_mat, type="full", order="original", p.mat = p_mat, sig.level = 0.05, insig = "blank", tl.cex = 0.5, tl.srt = 45, tl.col="black", 
         col = COL2('PuOr', 10))
mtext("Peripheral blood",side=3,line=3,cex=0.75,adj=0.58)
mtext("Upper respiratory",side=2,line=1,cex=0.75,adj=0.40)
dev.off()