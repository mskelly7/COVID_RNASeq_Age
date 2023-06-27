# BRAVE Kids RNA Sequencing Analysis
# Aditya Mohan (MD/PhD candidate) / Matthew Kelly, MD, MPH 
# DESeq2 models for differentially expressed genes
# Models are adjusted for sequencing batch, sex, cell populations (PC1-PC5), and (as applicable) age category

remove(list=ls())
setwd("G:/My Drive/Research/BRAVE_Kids/RNA_Sequencing/Statistical_Analyses") 
set.seed(1234)
version

library(readr)
library(dplyr)
library(DESeq2)
library(readxl)
library(writexl)
library(tidyr) 
library(ggpubr)
library(plyr) 
library(data.table)

# Load the required datasets
phy.rnaseq.np <- readRDS("G:/My Drive/Research/BRAVE_Kids/RNA_Sequencing/phy.rnaseq.np.rds")
metadata_np <- data.frame(sample_data(phy.rnaseq.np))
metadata_np$batch_num <- as.factor(metadata_np$batch_num)
metadata_np$age_cat <- as.factor(metadata_np$age_cat)
metadata_np$sex <- as.factor(metadata_np$sex)
metadata_np$corona <- as.factor(metadata_np$corona)
genes_np <- data.frame(tax_table(phy.rnaseq.np))
counts_np <- data.frame(otu_table(phy.rnaseq.np))
# Create datasets for NP samples among COVID-
metadata_np_Neg <- subset(metadata_np, corona=="Negative")
ids_np_Neg <- metadata_np_Neg$alias_sequencing_id
counts_np_Neg <- counts_np %>% dplyr::select(all_of(ids_np_Neg))
# Create datasets for NP samples among COVID+
metadata_np_Pos <- subset(metadata_np, corona=="Positive")
ids_np_Pos <- metadata_np_Pos$alias_sequencing_id
counts_np_Pos <- counts_np %>% dplyr::select(all_of(ids_np_Pos))
# Create datasets for NP samples among 0-5 years
metadata_np_0to5 <- subset(metadata_np, age_cat=="0-5 years")
ids_np_0to5 <- metadata_np_0to5$alias_sequencing_id
counts_np_0to5 <- counts_np %>% dplyr::select(all_of(ids_np_0to5))
# Create datasets for NP samples among 6-13 years
metadata_np_6to13 <- subset(metadata_np, age_cat=="6-13 years")
ids_np_6to13 <- metadata_np_6to13$alias_sequencing_id
counts_np_6to13 <- counts_np %>% dplyr::select(all_of(ids_np_6to13))
# Create datasets for NP samples among 14-20 years
metadata_np_14to20 <- subset(metadata_np, age_cat=="14-20 years")
ids_np_14to20 <- metadata_np_14to20$alias_sequencing_id
counts_np_14to20 <- counts_np %>% dplyr::select(all_of(ids_np_14to20))

phy.rnaseq.pax <- readRDS("G:/My Drive/Research/BRAVE_Kids/RNA_Sequencing/phy.rnaseq.pax.rds")
metadata_pax <- data.frame(sample_data(phy.rnaseq.pax))
metadata_pax$batch_num <- as.factor(metadata_pax$batch_num)
metadata_pax$age_cat <- as.factor(metadata_pax$age_cat)
metadata_pax$sex <- as.factor(metadata_pax$sex)
metadata_pax$corona <- as.factor(metadata_pax$corona)
genes_pax <- data.frame(tax_table(phy.rnaseq.pax))
counts_pax <- data.frame(otu_table(phy.rnaseq.pax))
names(counts_pax) <- gsub(x = names(counts_pax), pattern = "X025CBE.M0.PAX", replacement = "025CBE.M0.PAX")
names(counts_pax) <- gsub(x = names(counts_pax), pattern = "X0719FA.M0.PAX", replacement = "0719FA.M0.PAX")
names(counts_pax) <- gsub(x = names(counts_pax), pattern = "X0BF51C.M0.PAX", replacement = "0BF51C.M0.PAX")
names(counts_pax) <- gsub(x = names(counts_pax), pattern = "X0BF51C.M6.PAX", replacement = "0BF51C.M6.PAX")
names(counts_pax) <- gsub(x = names(counts_pax), pattern = "X0E1ECA.M0.PAX", replacement = "0E1ECA.M0.PAX")
names(counts_pax) <- gsub(x = names(counts_pax), pattern = "X18EF5E.M0.PAX", replacement = "18EF5E.M0.PAX")
names(counts_pax) <- gsub(x = names(counts_pax), pattern = "X1AC6B8.M0.PAX", replacement = "1AC6B8.M0.PAX")
names(counts_pax) <- gsub(x = names(counts_pax), pattern = "X1E06CA.M0.PAX", replacement = "1E06CA.M0.PAX")
names(counts_pax) <- gsub(x = names(counts_pax), pattern = "X234BBF.M0.PAX", replacement = "234BBF.M0.PAX")
names(counts_pax) <- gsub(x = names(counts_pax), pattern = "X322088.M0.PAX", replacement = "322088.M0.PAX")
names(counts_pax) <- gsub(x = names(counts_pax), pattern = "X3501A9.M0.PAX", replacement = "3501A9.M0.PAX")
names(counts_pax) <- gsub(x = names(counts_pax), pattern = "X3B441F.M0.PAX", replacement = "3B441F.M0.PAX")
names(counts_pax) <- gsub(x = names(counts_pax), pattern = "X496933.M0.PAX", replacement = "496933.M0.PAX")
names(counts_pax) <- gsub(x = names(counts_pax), pattern = "X58A568.M0.PAX", replacement = "58A568.M0.PAX")
names(counts_pax) <- gsub(x = names(counts_pax), pattern = "X6AA2F7.M0.PAX", replacement = "6AA2F7.M0.PAX")
names(counts_pax) <- gsub(x = names(counts_pax), pattern = "X71549F.M0.PAX", replacement = "71549F.M0.PAX")
names(counts_pax) <- gsub(x = names(counts_pax), pattern = "X72C971.M0.PAX", replacement = "72C971.M0.PAX")
names(counts_pax) <- gsub(x = names(counts_pax), pattern = "X7D880A.M0.PAX", replacement = "7D880A.M0.PAX")
names(counts_pax) <- gsub(x = names(counts_pax), pattern = "X826E5F.M0.PAX", replacement = "826E5F.M0.PAX")
names(counts_pax) <- gsub(x = names(counts_pax), pattern = "X830937.M0.PAX", replacement = "830937.M0.PAX")
names(counts_pax) <- gsub(x = names(counts_pax), pattern = "X8867DF.M0.PAX", replacement = "8867DF.M0.PAX")
names(counts_pax) <- gsub(x = names(counts_pax), pattern = "X88732B.M0.PAX", replacement = "88732B.M0.PAX")
names(counts_pax) <- gsub(x = names(counts_pax), pattern = "X88B192.M0.PAX", replacement = "88B192.M0.PAX")
# Create datasets for PAX samples among COVID-
metadata_pax_Neg <- subset(metadata_pax, corona=="Negative")
ids_pax_Neg <- metadata_pax_Neg$alias_sequencing_id
counts_pax_Neg <- counts_pax %>% dplyr::select(all_of(ids_pax_Neg))
# Create datasets for PAX samples among COVID+
metadata_pax_Pos <- subset(metadata_pax, corona=="Positive")
ids_pax_Pos <- metadata_pax_Pos$alias_sequencing_id
counts_pax_Pos <- counts_pax %>% dplyr::select(all_of(ids_pax_Pos))
# Create datasets for PAX samples among 0-5 years
metadata_pax_0to5 <- subset(metadata_pax, age_cat=="0-5 years")
ids_pax_0to5 <- metadata_pax_0to5$alias_sequencing_id
counts_pax_0to5 <- counts_pax %>% dplyr::select(all_of(ids_pax_0to5))
# Create datasets for PAX samples among 6-13 years
metadata_pax_6to13 <- subset(metadata_pax, age_cat=="6-13 years")
ids_pax_6to13 <- metadata_pax_6to13$alias_sequencing_id
counts_pax_6to13 <- counts_pax %>% dplyr::select(all_of(ids_pax_6to13))
# Create datasets for PAX samples among 14-20 years
metadata_pax_14to20 <- subset(metadata_pax, age_cat=="14-20 years")
ids_pax_14to20 <- metadata_pax_14to20$alias_sequencing_id
counts_pax_14to20 <- counts_pax %>% dplyr::select(all_of(ids_pax_14to20))
# Create datasets for PAX samples among Adult
metadata_pax_adult <- subset(metadata_pax, age_cat=="Adult")
ids_pax_adult <- metadata_pax_adult$alias_sequencing_id
counts_pax_adult <- counts_pax %>% dplyr::select(all_of(ids_pax_adult))

############################################################################
### MODEL 1 - DIFFERENCES IN NP GENE EXPRESSION BY AGE IN COVID NEGATIVE ###
############################################################################

nrow(metadata_np_Neg)
all(rownames(metadata_np_Neg) == colnames(counts_np_Neg))
dds_np_Neg <- DESeqDataSetFromMatrix(countData = counts_np_Neg, colData = metadata_np_Neg, design = ~ batch_num + age_cat + sex + PC1 + PC2 + PC3 + PC4 + PC5)

# Require a normalized count of at least 10 in two or more samples
dds_np_Neg_norm <- estimateSizeFactors(dds_np_Neg)
dds_np_Neg_norm <- counts(dds_np_Neg_norm, normalized=TRUE)
filter <- rowSums(dds_np_Neg_norm >= 10) >= (nrow(metadata_np_Neg))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_np_Neg <- dds_np_Neg[filter,]
dds_np_Neg <- DESeq(dds_np_Neg)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_np_Neg <- dds_np_Neg[which(mcols(dds_np_Neg)$betaConv),]
resultsNames(dds_np_Neg) # lists the coefficients

# Output results for 0-5 years vs. 6-13 years = genes with logFC>0 are overexpressed in children 0-5 years
dds_np_Neg_0to5_6to13 <- results(dds_np_Neg, contrast=c("age_cat","0-5 years","6-13 years"))
dds_np_Neg_0to5_6to13 <- data.frame(EnsemblID=row.names(dds_np_Neg_0to5_6to13), dds_np_Neg_0to5_6to13)
dds_np_Neg_0to5_6to13 <- merge(dds_np_Neg_0to5_6to13, genes_np, by="row.names")
dds_np_Neg_0to5_6to13 <- dds_np_Neg_0to5_6to13[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_np_Neg_0to5_6to13 <- subset(dds_np_Neg_0to5_6to13, !is.na(padj))
dds_np_Neg_0to5_6to13 <- dds_np_Neg_0to5_6to13[order(dds_np_Neg_0to5_6to13$padj),]
dds_np_Neg_0to5_6to13 <- dds_np_Neg_0to5_6to13 %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_np_Neg_0to5_6to13 <- dds_np_Neg_0to5_6to13[order(dds_np_Neg_0to5_6to13$Expression, dds_np_Neg_0to5_6to13$log2FoldChange),]
write.csv(dds_np_Neg_0to5_6to13, "1_COVID_Neg_by_Age/genes_np_Neg_0to5_6to13.csv", row.names=FALSE)

# Output results for 0-5 years vs. 14-20 years = genes with logFC>0 are overexpressed in children 0-5 years
# results(dds, contrast=c("condition","C","B")) = genes with logFC > 0 are overexpressed in C
dds_np_Neg_0to5_14to20 <- results(dds_np_Neg, contrast=c("age_cat","0-5 years","14-20 years"))
dds_np_Neg_0to5_14to20 <- data.frame(EnsemblID=row.names(dds_np_Neg_0to5_14to20), dds_np_Neg_0to5_14to20)
dds_np_Neg_0to5_14to20 <- merge(dds_np_Neg_0to5_14to20, genes_np, by="row.names")
dds_np_Neg_0to5_14to20 <- dds_np_Neg_0to5_14to20[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_np_Neg_0to5_14to20 <- subset(dds_np_Neg_0to5_14to20, !is.na(padj))
dds_np_Neg_0to5_14to20 <- dds_np_Neg_0to5_14to20[order(dds_np_Neg_0to5_14to20$padj),]
dds_np_Neg_0to5_14to20 <- dds_np_Neg_0to5_14to20 %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_np_Neg_0to5_14to20 <- dds_np_Neg_0to5_14to20[order(dds_np_Neg_0to5_14to20$Expression, dds_np_Neg_0to5_14to20$log2FoldChange),]
write.csv(dds_np_Neg_0to5_14to20, "1_COVID_Neg_by_Age/genes_np_Neg_0to5_14to20.csv", row.names=FALSE)

# Output results for 6-13 years vs. 14-20 years = genes with logFC>0 are overexpressed in children 6-13 years
# results(dds, contrast=c("condition","C","B")) = genes with logFC > 0 are overexpressed in C
dds_np_Neg_6to13_14to20 <- results(dds_np_Neg, contrast=c("age_cat","6-13 years","14-20 years"))
dds_np_Neg_6to13_14to20 <- data.frame(EnsemblID=row.names(dds_np_Neg_6to13_14to20), dds_np_Neg_6to13_14to20)
dds_np_Neg_6to13_14to20 <- merge(dds_np_Neg_6to13_14to20, genes_np, by="row.names")
dds_np_Neg_6to13_14to20 <- dds_np_Neg_6to13_14to20[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_np_Neg_6to13_14to20 <- subset(dds_np_Neg_6to13_14to20, !is.na(padj))
dds_np_Neg_6to13_14to20 <- dds_np_Neg_6to13_14to20[order(dds_np_Neg_6to13_14to20$padj),]
dds_np_Neg_6to13_14to20 <- dds_np_Neg_6to13_14to20 %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_np_Neg_6to13_14to20 <- dds_np_Neg_6to13_14to20[order(dds_np_Neg_6to13_14to20$Expression, dds_np_Neg_6to13_14to20$log2FoldChange),]
write.csv(dds_np_Neg_6to13_14to20, "1_COVID_Neg_by_Age/genes_np_Neg_6to13_14to20.csv", row.names=FALSE)

#############################################################################
### MODEL 2 - DIFFERENCES IN PAX GENE EXPRESSION BY AGE IN COVID NEGATIVE ###
#############################################################################

nrow(metadata_pax_Neg)
all(rownames(metadata_pax_Neg) == colnames(counts_pax_Neg))
dds_pax_Neg <- DESeqDataSetFromMatrix(countData = counts_pax_Neg, colData = metadata_pax_Neg, design = ~ batch_num + age_cat + sex + PC1 + PC2 + PC3 + PC4 + PC5)

# Require a normalized count of at least 10 in two or more samples
dds_pax_Neg_norm <- estimateSizeFactors(dds_pax_Neg)
dds_pax_Neg_norm <- counts(dds_pax_Neg_norm, normalized=TRUE)
filter <- rowSums(dds_pax_Neg_norm >= 10) >= (nrow(metadata_pax_Neg))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_pax_Neg <- dds_pax_Neg[filter,]
dds_pax_Neg <- DESeq(dds_pax_Neg)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_pax_Neg <- dds_pax_Neg[which(mcols(dds_pax_Neg)$betaConv),]
resultsNames(dds_pax_Neg) # lists the coefficients

# Output results for 0-5 years vs. 6-13 years = genes with logFC>0 are overexpressed in children 0-5 years
dds_pax_Neg_0to5_6to13 <- results(dds_pax_Neg, contrast=c("age_cat","0-5 years","6-13 years"))
dds_pax_Neg_0to5_6to13 <- data.frame(EnsemblID=row.names(dds_pax_Neg_0to5_6to13), dds_pax_Neg_0to5_6to13)
dds_pax_Neg_0to5_6to13 <- merge(dds_pax_Neg_0to5_6to13, genes_pax, by="row.names")
dds_pax_Neg_0to5_6to13 <- dds_pax_Neg_0to5_6to13[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_pax_Neg_0to5_6to13 <- subset(dds_pax_Neg_0to5_6to13, !is.na(padj))
dds_pax_Neg_0to5_6to13 <- dds_pax_Neg_0to5_6to13[order(dds_pax_Neg_0to5_6to13$padj),]
dds_pax_Neg_0to5_6to13 <- dds_pax_Neg_0to5_6to13 %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_pax_Neg_0to5_6to13 <- dds_pax_Neg_0to5_6to13[order(dds_pax_Neg_0to5_6to13$Expression, dds_pax_Neg_0to5_6to13$log2FoldChange),]
write.csv(dds_pax_Neg_0to5_6to13, "1_COVID_Neg_by_Age/genes_pax_Neg_0to5_6to13.csv", row.names=FALSE)

# Output results for 0-5 years vs. 14-20 years = genes with logFC>0 are overexpressed in children 0-5 years
# results(dds, contrast=c("condition","C","B")) = genes with logFC > 0 are overexpressed in C
dds_pax_Neg_0to5_14to20 <- results(dds_pax_Neg, contrast=c("age_cat","0-5 years","14-20 years"))
dds_pax_Neg_0to5_14to20 <- data.frame(EnsemblID=row.names(dds_pax_Neg_0to5_14to20), dds_pax_Neg_0to5_14to20)
dds_pax_Neg_0to5_14to20 <- merge(dds_pax_Neg_0to5_14to20, genes_pax, by="row.names")
dds_pax_Neg_0to5_14to20 <- dds_pax_Neg_0to5_14to20[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_pax_Neg_0to5_14to20 <- subset(dds_pax_Neg_0to5_14to20, !is.na(padj))
dds_pax_Neg_0to5_14to20 <- dds_pax_Neg_0to5_14to20[order(dds_pax_Neg_0to5_14to20$padj),]
dds_pax_Neg_0to5_14to20 <- dds_pax_Neg_0to5_14to20 %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_pax_Neg_0to5_14to20 <- dds_pax_Neg_0to5_14to20[order(dds_pax_Neg_0to5_14to20$Expression, dds_pax_Neg_0to5_14to20$log2FoldChange),]
write.csv(dds_pax_Neg_0to5_14to20, "1_COVID_Neg_by_Age/genes_pax_Neg_0to5_14to20.csv", row.names=FALSE)

# Output results for 6-13 years vs. 14-20 years = genes with logFC>0 are overexpressed in children 6-13 years
# results(dds, contrast=c("condition","C","B")) = genes with logFC > 0 are overexpressed in C
dds_pax_Neg_6to13_14to20 <- results(dds_pax_Neg, contrast=c("age_cat","6-13 years","14-20 years"))
dds_pax_Neg_6to13_14to20 <- data.frame(EnsemblID=row.names(dds_pax_Neg_6to13_14to20), dds_pax_Neg_6to13_14to20)
dds_pax_Neg_6to13_14to20 <- merge(dds_pax_Neg_6to13_14to20, genes_pax, by="row.names")
dds_pax_Neg_6to13_14to20 <- dds_pax_Neg_6to13_14to20[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_pax_Neg_6to13_14to20 <- subset(dds_pax_Neg_6to13_14to20, !is.na(padj))
dds_pax_Neg_6to13_14to20 <- dds_pax_Neg_6to13_14to20[order(dds_pax_Neg_6to13_14to20$padj),]
dds_pax_Neg_6to13_14to20 <- dds_pax_Neg_6to13_14to20 %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_pax_Neg_6to13_14to20 <- dds_pax_Neg_6to13_14to20[order(dds_pax_Neg_6to13_14to20$Expression, dds_pax_Neg_6to13_14to20$log2FoldChange),]
write.csv(dds_pax_Neg_6to13_14to20, "1_COVID_Neg_by_Age/genes_pax_Neg_6to13_14to20.csv", row.names=FALSE)

##############################################################
### MODEL 3 - DIFFERENCES IN NP EXPRESSION BY COVID STATUS ###
##############################################################

nrow(metadata_np)
all(rownames(metadata_np) == colnames(counts_np))
# Would need to scale age variable if included as continuous in models
#metadata_np$age_scale <- scale(metadata_np$age, center = TRUE, scale=TRUE)
dds_np <- DESeqDataSetFromMatrix(countData = counts_np, colData = metadata_np, design = ~ batch_num + age_cat + sex + corona + PC1 + PC2 + PC3 + PC4 + PC5)

# Require a normalized count of at least 10 in two or more samples
dds_np_norm <- estimateSizeFactors(dds_np)
dds_np_norm <- counts(dds_np_norm, normalized=TRUE)
filter <- rowSums(dds_np_norm >= 10) >= (nrow(metadata_np))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_np <- dds_np[filter,]
dds_np <- DESeq(dds_np)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_np <- dds_np[which(mcols(dds_np)$betaConv),]
resultsNames(dds_np) # lists the coefficients

# Output results for COVID status
dds_np_Pos_Neg <- results(dds_np, contrast=c("corona","Positive","Negative"))
dds_np_Pos_Neg <- data.frame(EnsemblID=row.names(dds_np_Pos_Neg), dds_np_Pos_Neg)
dds_np_Pos_Neg <- merge(dds_np_Pos_Neg, genes_np, by="row.names")
dds_np_Pos_Neg <- dds_np_Pos_Neg[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_np_Pos_Neg <- subset(dds_np_Pos_Neg, !is.na(padj))
dds_np_Pos_Neg <- dds_np_Pos_Neg[order(dds_np_Pos_Neg$padj),]
dds_np_Pos_Neg <- dds_np_Pos_Neg %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_np_Pos_Neg <- dds_np_Pos_Neg[order(dds_np_Pos_Neg$Expression, dds_np_Pos_Neg$log2FoldChange),]
write.csv(dds_np_Pos_Neg, "2_COVID_Pos_vs_Neg/genes_np_Pos_Neg.csv", row.names=FALSE)

####################################################################
### MODEL 4 - DIFFERENCES IN PAX GENE EXPRESSION BY COVID STATUS ###
####################################################################

nrow(metadata_pax)
all(rownames(metadata_pax) == colnames(counts_pax))
# Would need to scale age variable if included as continuous in models
#metadata_pax$age_scale <- scale(metadata_pax$age, center = TRUE, scale=TRUE)
dds_pax <- DESeqDataSetFromMatrix(countData = counts_pax, colData = metadata_pax, design = ~ batch_num + age_cat + sex + corona + PC1 + PC2 + PC3 + PC4 + PC5)

# Require a normalized count of at least 10 in two or more samples
dds_pax_norm <- estimateSizeFactors(dds_pax)
dds_pax_norm <- counts(dds_pax_norm, normalized=TRUE)
filter <- rowSums(dds_pax_norm >= 10) >= (nrow(metadata_pax))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_pax <- dds_pax[filter,]
dds_pax <- DESeq(dds_pax)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_pax <- dds_pax[which(mcols(dds_pax)$betaConv),]
resultsNames(dds_pax) # lists the coefficients

# Output results for COVID status
dds_pax_Pos_Neg <- results(dds_pax, contrast=c("corona","Positive","Negative"))
dds_pax_Pos_Neg <- data.frame(EnsemblID=row.names(dds_pax_Pos_Neg), dds_pax_Pos_Neg)
dds_pax_Pos_Neg <- merge(dds_pax_Pos_Neg, genes_pax, by="row.names")
dds_pax_Pos_Neg <- dds_pax_Pos_Neg[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_pax_Pos_Neg <- subset(dds_pax_Pos_Neg, !is.na(padj))
dds_pax_Pos_Neg <- dds_pax_Pos_Neg[order(dds_pax_Pos_Neg$padj),]
dds_pax_Pos_Neg <- dds_pax_Pos_Neg %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_pax_Pos_Neg <- dds_pax_Pos_Neg[order(dds_pax_Pos_Neg$Expression, dds_pax_Pos_Neg$log2FoldChange),]
write.csv(dds_pax_Pos_Neg, "2_COVID_Pos_vs_Neg/genes_pax_Pos_Neg.csv", row.names=FALSE)

###########################################################################
### MODELS 5-7 - COMPARING NP SAMPLES BY COVID STATUS IN AGE CATEGORIES ###
###########################################################################

# 0-5 years

nrow(metadata_np_0to5)
all(rownames(metadata_np_0to5) == colnames(counts_np_0to5))
dds_np_0to5 <- DESeqDataSetFromMatrix(countData = counts_np_0to5, colData = metadata_np_0to5, design = ~ batch_num + sex + corona + PC1 + PC2 + PC3 + PC4 + PC5)

# Require a normalized count of at least 10 in two or more samples
dds_np_0to5_norm <- estimateSizeFactors(dds_np_0to5)
dds_np_0to5_norm <- counts(dds_np_0to5_norm, normalized=TRUE)
filter <- rowSums(dds_np_0to5_norm >= 10) >= (nrow(metadata_np_0to5))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_np_0to5 <- dds_np_0to5[filter,]
dds_np_0to5 <- DESeq(dds_np_0to5)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_np_0to5 <- dds_np_0to5[which(mcols(dds_np_0to5)$betaConv),]
resultsNames(dds_np_0to5) # lists the coefficients

# Output results for COVID status
dds_np_0to5_Pos_Neg <- results(dds_np_0to5, contrast=c("corona","Positive","Negative"))
dds_np_0to5_Pos_Neg <- data.frame(EnsemblID=row.names(dds_np_0to5_Pos_Neg), dds_np_0to5_Pos_Neg)
dds_np_0to5_Pos_Neg <- merge(dds_np_0to5_Pos_Neg, genes_np, by="row.names")
dds_np_0to5_Pos_Neg <- dds_np_0to5_Pos_Neg[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_np_0to5_Pos_Neg <- subset(dds_np_0to5_Pos_Neg, !is.na(padj))
dds_np_0to5_Pos_Neg <- dds_np_0to5_Pos_Neg[order(dds_np_0to5_Pos_Neg$padj),]
dds_np_0to5_Pos_Neg <- dds_np_0to5_Pos_Neg %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_np_0to5_Pos_Neg <- dds_np_0to5_Pos_Neg[order(dds_np_0to5_Pos_Neg$Expression, dds_np_0to5_Pos_Neg$log2FoldChange),]
write.csv(dds_np_0to5_Pos_Neg, "3_COVID_Pos_by_Age/genes_np_0to5_Pos_Neg.csv", row.names=FALSE)

# 6-13 years

nrow(metadata_np_6to13)
all(rownames(metadata_np_6to13) == colnames(counts_np_6to13))
dds_np_6to13 <- DESeqDataSetFromMatrix(countData = counts_np_6to13, colData = metadata_np_6to13, design = ~ batch_num + sex + corona + PC1 + PC2 + PC3 + PC4 + PC5)

# Require a normalized count of at least 10 in two or more samples
dds_np_6to13_norm <- estimateSizeFactors(dds_np_6to13)
dds_np_6to13_norm <- counts(dds_np_6to13_norm, normalized=TRUE)
filter <- rowSums(dds_np_6to13_norm >= 10) >= (nrow(metadata_np_6to13))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_np_6to13 <- dds_np_6to13[filter,]
dds_np_6to13 <- DESeq(dds_np_6to13)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_np_6to13 <- dds_np_6to13[which(mcols(dds_np_6to13)$betaConv),]
resultsNames(dds_np_6to13) # lists the coefficients

# Output results for COVID status
dds_np_6to13_Pos_Neg <- results(dds_np_6to13, contrast=c("corona","Positive","Negative"))
dds_np_6to13_Pos_Neg <- data.frame(EnsemblID=row.names(dds_np_6to13_Pos_Neg), dds_np_6to13_Pos_Neg)
dds_np_6to13_Pos_Neg <- merge(dds_np_6to13_Pos_Neg, genes_np, by="row.names")
dds_np_6to13_Pos_Neg <- dds_np_6to13_Pos_Neg[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_np_6to13_Pos_Neg <- subset(dds_np_6to13_Pos_Neg, !is.na(padj))
dds_np_6to13_Pos_Neg <- dds_np_6to13_Pos_Neg[order(dds_np_6to13_Pos_Neg$padj),]
dds_np_6to13_Pos_Neg <- dds_np_6to13_Pos_Neg %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_np_6to13_Pos_Neg <- dds_np_6to13_Pos_Neg[order(dds_np_6to13_Pos_Neg$Expression, dds_np_6to13_Pos_Neg$log2FoldChange),]
write.csv(dds_np_6to13_Pos_Neg, "3_COVID_Pos_by_Age/genes_np_6to13_Pos_Neg.csv", row.names=FALSE)

# 14-20 years

nrow(metadata_np_14to20)
all(rownames(metadata_np_14to20) == colnames(counts_np_14to20))
dds_np_14to20 <- DESeqDataSetFromMatrix(countData = counts_np_14to20, colData = metadata_np_14to20, design = ~ batch_num + sex + corona + PC1 + PC2 + PC3 + PC4 + PC5)

# Require a normalized count of at least 10 in two or more samples
dds_np_14to20_norm <- estimateSizeFactors(dds_np_14to20)
dds_np_14to20_norm <- counts(dds_np_14to20_norm, normalized=TRUE)
filter <- rowSums(dds_np_14to20_norm >= 10) >= (nrow(metadata_np_14to20))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_np_14to20 <- dds_np_14to20[filter,]
dds_np_14to20 <- DESeq(dds_np_14to20)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_np_14to20 <- dds_np_14to20[which(mcols(dds_np_14to20)$betaConv),]
resultsNames(dds_np_14to20) # lists the coefficients

# Output results for COVID status
dds_np_14to20_Pos_Neg <- results(dds_np_14to20, contrast=c("corona","Positive","Negative"))
dds_np_14to20_Pos_Neg <- data.frame(EnsemblID=row.names(dds_np_14to20_Pos_Neg), dds_np_14to20_Pos_Neg)
dds_np_14to20_Pos_Neg <- merge(dds_np_14to20_Pos_Neg, genes_np, by="row.names")
dds_np_14to20_Pos_Neg <- dds_np_14to20_Pos_Neg[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_np_14to20_Pos_Neg <- subset(dds_np_14to20_Pos_Neg, !is.na(padj))
dds_np_14to20_Pos_Neg <- dds_np_14to20_Pos_Neg[order(dds_np_14to20_Pos_Neg$padj),]
dds_np_14to20_Pos_Neg <- dds_np_14to20_Pos_Neg %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_np_14to20_Pos_Neg <- dds_np_14to20_Pos_Neg[order(dds_np_14to20_Pos_Neg$Expression, dds_np_14to20_Pos_Neg$log2FoldChange),]
write.csv(dds_np_14to20_Pos_Neg, "3_COVID_Pos_by_Age/genes_np_14to20_Pos_Neg.csv", row.names=FALSE)

#############################################################################
### MODELS 8-11 - COMPARING PAX SAMPLES BY COVID STATUS IN AGE CATEGORIES ###
#############################################################################

# 0-5 years

nrow(metadata_pax_0to5)
all(rownames(metadata_pax_0to5) == colnames(counts_pax_0to5))
dds_pax_0to5 <- DESeqDataSetFromMatrix(countData = counts_pax_0to5, colData = metadata_pax_0to5, design = ~ batch_num + sex + corona + PC1 + PC2 + PC3 + PC4 + PC5)

# Require a normalized count of at least 10 in two or more samples
dds_pax_0to5_norm <- estimateSizeFactors(dds_pax_0to5)
dds_pax_0to5_norm <- counts(dds_pax_0to5_norm, normalized=TRUE)
filter <- rowSums(dds_pax_0to5_norm >= 10) >= (nrow(metadata_pax_0to5))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_pax_0to5 <- dds_pax_0to5[filter,]
dds_pax_0to5 <- DESeq(dds_pax_0to5)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_pax_0to5 <- dds_pax_0to5[which(mcols(dds_pax_0to5)$betaConv),]
resultsNames(dds_pax_0to5) # lists the coefficients

# Output results for COVID status
dds_pax_0to5_Pos_Neg <- results(dds_pax_0to5, contrast=c("corona","Positive","Negative"))
dds_pax_0to5_Pos_Neg <- data.frame(EnsemblID=row.names(dds_pax_0to5_Pos_Neg), dds_pax_0to5_Pos_Neg)
dds_pax_0to5_Pos_Neg <- merge(dds_pax_0to5_Pos_Neg, genes_pax, by="row.names")
dds_pax_0to5_Pos_Neg <- dds_pax_0to5_Pos_Neg[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_pax_0to5_Pos_Neg <- subset(dds_pax_0to5_Pos_Neg, !is.na(padj))
dds_pax_0to5_Pos_Neg <- dds_pax_0to5_Pos_Neg[order(dds_pax_0to5_Pos_Neg$padj),]
dds_pax_0to5_Pos_Neg <- dds_pax_0to5_Pos_Neg %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_pax_0to5_Pos_Neg <- dds_pax_0to5_Pos_Neg[order(dds_pax_0to5_Pos_Neg$Expression, dds_pax_0to5_Pos_Neg$log2FoldChange),]
write.csv(dds_pax_0to5_Pos_Neg, "3_COVID_Pos_by_Age/genes_pax_0to5_Pos_Neg.csv", row.names=FALSE)

# 6-13 years

nrow(metadata_pax_6to13)
all(rownames(metadata_pax_6to13) == colnames(counts_pax_6to13))
dds_pax_6to13 <- DESeqDataSetFromMatrix(countData = counts_pax_6to13, colData = metadata_pax_6to13, design = ~ batch_num + sex + corona + PC1 + PC2 + PC3 + PC4 + PC5)

# Require a normalized count of at least 10 in two or more samples
dds_pax_6to13_norm <- estimateSizeFactors(dds_pax_6to13)
dds_pax_6to13_norm <- counts(dds_pax_6to13_norm, normalized=TRUE)
filter <- rowSums(dds_pax_6to13_norm >= 10) >= (nrow(metadata_pax_6to13))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_pax_6to13 <- dds_pax_6to13[filter,]
dds_pax_6to13 <- DESeq(dds_pax_6to13)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_pax_6to13 <- dds_pax_6to13[which(mcols(dds_pax_6to13)$betaConv),]
resultsNames(dds_pax_6to13) # lists the coefficients

# Output results for COVID status
dds_pax_6to13_Pos_Neg <- results(dds_pax_6to13, contrast=c("corona","Positive","Negative"))
dds_pax_6to13_Pos_Neg <- data.frame(EnsemblID=row.names(dds_pax_6to13_Pos_Neg), dds_pax_6to13_Pos_Neg)
dds_pax_6to13_Pos_Neg <- merge(dds_pax_6to13_Pos_Neg, genes_pax, by="row.names")
dds_pax_6to13_Pos_Neg <- dds_pax_6to13_Pos_Neg[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_pax_6to13_Pos_Neg <- subset(dds_pax_6to13_Pos_Neg, !is.na(padj))
dds_pax_6to13_Pos_Neg <- dds_pax_6to13_Pos_Neg[order(dds_pax_6to13_Pos_Neg$padj),]
dds_pax_6to13_Pos_Neg <- dds_pax_6to13_Pos_Neg %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_pax_6to13_Pos_Neg <- dds_pax_6to13_Pos_Neg[order(dds_pax_6to13_Pos_Neg$Expression, dds_pax_6to13_Pos_Neg$log2FoldChange),]
write.csv(dds_pax_6to13_Pos_Neg, "3_COVID_Pos_by_Age/genes_pax_6to13_Pos_Neg.csv", row.names=FALSE)

# 14-20 years

nrow(metadata_pax_14to20)
all(rownames(metadata_pax_14to20) == colnames(counts_pax_14to20))
dds_pax_14to20 <- DESeqDataSetFromMatrix(countData = counts_pax_14to20, colData = metadata_pax_14to20, design = ~ batch_num + sex + corona + PC1 + PC2 + PC3 + PC4 + PC5)

# Require a normalized count of at least 10 in two or more samples
dds_pax_14to20_norm <- estimateSizeFactors(dds_pax_14to20)
dds_pax_14to20_norm <- counts(dds_pax_14to20_norm, normalized=TRUE)
filter <- rowSums(dds_pax_14to20_norm >= 10) >= (nrow(metadata_pax_14to20))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_pax_14to20 <- dds_pax_14to20[filter,]
dds_pax_14to20 <- DESeq(dds_pax_14to20)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_pax_14to20 <- dds_pax_14to20[which(mcols(dds_pax_14to20)$betaConv),]
resultsNames(dds_pax_14to20) # lists the coefficients

# Output results for COVID status
dds_pax_14to20_Pos_Neg <- results(dds_pax_14to20, contrast=c("corona","Positive","Negative"))
dds_pax_14to20_Pos_Neg <- data.frame(EnsemblID=row.names(dds_pax_14to20_Pos_Neg), dds_pax_14to20_Pos_Neg)
dds_pax_14to20_Pos_Neg <- merge(dds_pax_14to20_Pos_Neg, genes_pax, by="row.names")
dds_pax_14to20_Pos_Neg <- dds_pax_14to20_Pos_Neg[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_pax_14to20_Pos_Neg <- subset(dds_pax_14to20_Pos_Neg, !is.na(padj))
dds_pax_14to20_Pos_Neg <- dds_pax_14to20_Pos_Neg[order(dds_pax_14to20_Pos_Neg$padj),]
dds_pax_14to20_Pos_Neg <- dds_pax_14to20_Pos_Neg %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_pax_14to20_Pos_Neg <- dds_pax_14to20_Pos_Neg[order(dds_pax_14to20_Pos_Neg$Expression, dds_pax_14to20_Pos_Neg$log2FoldChange),]
write.csv(dds_pax_14to20_Pos_Neg, "3_COVID_Pos_by_Age/genes_pax_14to20_Pos_Neg.csv", row.names=FALSE)

# Adult

nrow(metadata_pax_adult)
all(rownames(metadata_pax_adult) == colnames(counts_pax_adult))
dds_pax_adult <- DESeqDataSetFromMatrix(countData = counts_pax_adult, colData = metadata_pax_adult, design = ~ batch_num + sex + corona + PC1 + PC2 + PC3 + PC4 + PC5)

# Require a normalized count of at least 10 in two or more samples
dds_pax_adult_norm <- estimateSizeFactors(dds_pax_adult)
dds_pax_adult_norm <- counts(dds_pax_adult_norm, normalized=TRUE)
filter <- rowSums(dds_pax_adult_norm >= 10) >= (nrow(metadata_pax_adult))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_pax_adult <- dds_pax_adult[filter,]
dds_pax_adult <- DESeq(dds_pax_adult)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_pax_adult <- dds_pax_adult[which(mcols(dds_pax_adult)$betaConv),]
resultsNames(dds_pax_adult) # lists the coefficients

# Output results for COVID status
dds_pax_adult_Pos_Neg <- results(dds_pax_adult, contrast=c("corona","Positive","Negative"))
dds_pax_adult_Pos_Neg <- data.frame(EnsemblID=row.names(dds_pax_adult_Pos_Neg), dds_pax_adult_Pos_Neg)
dds_pax_adult_Pos_Neg <- merge(dds_pax_adult_Pos_Neg, genes_pax, by="row.names")
dds_pax_adult_Pos_Neg <- dds_pax_adult_Pos_Neg[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_pax_adult_Pos_Neg <- subset(dds_pax_adult_Pos_Neg, !is.na(padj))
dds_pax_adult_Pos_Neg <- dds_pax_adult_Pos_Neg[order(dds_pax_adult_Pos_Neg$padj),]
dds_pax_adult_Pos_Neg <- dds_pax_adult_Pos_Neg %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_pax_adult_Pos_Neg <- dds_pax_adult_Pos_Neg[order(dds_pax_adult_Pos_Neg$Expression, dds_pax_adult_Pos_Neg$log2FoldChange),]
write.csv(dds_pax_adult_Pos_Neg, "3_COVID_Pos_by_Age/genes_pax_adult_Pos_Neg.csv", row.names=FALSE)

###############################################################################
### MODELS FOR PATIENT CHARACTERISTICS & SYMPTOMS AMONG COVID+ - NP SAMPLES ###
###############################################################################

# Obesity
metadata_np_obesity <- subset(metadata_np_Pos, !is.na(obesity))
ids_np_obesity <- metadata_np_obesity$alias_sequencing_id
counts_np_obesity <- counts_np %>% dplyr::select(all_of(ids_np_obesity))
nrow(metadata_np_obesity)
all(rownames(metadata_np_obesity) == colnames(counts_np_obesity))
dds_np_obesity <- DESeqDataSetFromMatrix(countData = counts_np_obesity, colData = metadata_np_obesity, design = ~ batch_num + sex + age_cat + obesity + PC1 + PC2 + PC3 + PC4 + PC5)
# Require a normalized count of at least 10 in two or more samples
dds_np_obesity_norm <- estimateSizeFactors(dds_np_obesity)
dds_np_obesity_norm <- counts(dds_np_obesity_norm, normalized=TRUE)
filter <- rowSums(dds_np_obesity_norm >= 10) >= (nrow(metadata_np_obesity))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_np_obesity <- dds_np_obesity[filter,]
dds_np_obesity <- DESeq(dds_np_obesity)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_np_obesity <- dds_np_obesity[which(mcols(dds_np_obesity)$betaConv),]
resultsNames(dds_np_obesity) # lists the coefficients
dds_np_obesity <- results(dds_np_obesity, contrast=c("obesity","Y","N"))
dds_np_obesity <- data.frame(EnsemblID=row.names(dds_np_obesity), dds_np_obesity)
dds_np_obesity <- merge(dds_np_obesity, genes_np, by="row.names")
dds_np_obesity <- dds_np_obesity[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_np_obesity <- subset(dds_np_obesity, !is.na(padj))
dds_np_obesity <- dds_np_obesity[order(dds_np_obesity$padj),]
dds_np_obesity <- dds_np_obesity %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_np_obesity <- dds_np_obesity[order(dds_np_obesity$Expression, dds_np_obesity$log2FoldChange),]
write.csv(dds_np_obesity, "4_Symptoms/genes_np_obesity.csv", row.names=FALSE)

# Asthma
metadata_np_asthma <- subset(metadata_np_Pos, !is.na(asthma))
ids_np_asthma <- metadata_np_asthma$alias_sequencing_id
counts_np_asthma <- counts_np %>% dplyr::select(all_of(ids_np_asthma))
nrow(metadata_np_asthma)
all(rownames(metadata_np_asthma) == colnames(counts_np_asthma))
dds_np_asthma <- DESeqDataSetFromMatrix(countData = counts_np_asthma, colData = metadata_np_asthma, design = ~ batch_num + sex + age_cat + asthma + PC1 + PC2 + PC3 + PC4 + PC5)
# Require a normalized count of at least 10 in two or more samples
dds_np_asthma_norm <- estimateSizeFactors(dds_np_asthma)
dds_np_asthma_norm <- counts(dds_np_asthma_norm, normalized=TRUE)
filter <- rowSums(dds_np_asthma_norm >= 10) >= (nrow(metadata_np_asthma))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_np_asthma <- dds_np_asthma[filter,]
dds_np_asthma <- DESeq(dds_np_asthma)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_np_asthma <- dds_np_asthma[which(mcols(dds_np_asthma)$betaConv),]
resultsNames(dds_np_asthma) # lists the coefficients
dds_np_asthma <- results(dds_np_asthma, contrast=c("asthma","Y","N"))
dds_np_asthma <- data.frame(EnsemblID=row.names(dds_np_asthma), dds_np_asthma)
dds_np_asthma <- merge(dds_np_asthma, genes_np, by="row.names")
dds_np_asthma <- dds_np_asthma[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_np_asthma <- subset(dds_np_asthma, !is.na(padj))
dds_np_asthma <- dds_np_asthma[order(dds_np_asthma$padj),]
dds_np_asthma <- dds_np_asthma %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_np_asthma <- dds_np_asthma[order(dds_np_asthma$Expression, dds_np_asthma$log2FoldChange),]
write.csv(dds_np_asthma, "4_Symptoms/genes_np_asthma.csv", row.names=FALSE)

# Fever
metadata_np_fever <- subset(metadata_np_Pos, !is.na(fever))
ids_np_fever <- metadata_np_fever$alias_sequencing_id
counts_np_fever <- counts_np %>% dplyr::select(all_of(ids_np_fever))
nrow(metadata_np_fever)
all(rownames(metadata_np_fever) == colnames(counts_np_fever))
dds_np_fever <- DESeqDataSetFromMatrix(countData = counts_np_fever, colData = metadata_np_fever, design = ~ batch_num + sex + age_cat + fever + PC1 + PC2 + PC3 + PC4 + PC5)
# Require a normalized count of at least 10 in two or more samples
dds_np_fever_norm <- estimateSizeFactors(dds_np_fever)
dds_np_fever_norm <- counts(dds_np_fever_norm, normalized=TRUE)
filter <- rowSums(dds_np_fever_norm >= 10) >= (nrow(metadata_np_fever))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_np_fever <- dds_np_fever[filter,]
dds_np_fever <- DESeq(dds_np_fever)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_np_fever <- dds_np_fever[which(mcols(dds_np_fever)$betaConv),]
resultsNames(dds_np_fever) # lists the coefficients
dds_np_fever <- results(dds_np_fever, contrast=c("fever","Y","N"))
dds_np_fever <- data.frame(EnsemblID=row.names(dds_np_fever), dds_np_fever)
dds_np_fever <- merge(dds_np_fever, genes_np, by="row.names")
dds_np_fever <- dds_np_fever[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_np_fever <- subset(dds_np_fever, !is.na(padj))
dds_np_fever <- dds_np_fever[order(dds_np_fever$padj),]
dds_np_fever <- dds_np_fever %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_np_fever <- dds_np_fever[order(dds_np_fever$Expression, dds_np_fever$log2FoldChange),]
write.csv(dds_np_fever, "4_Symptoms/genes_np_fever.csv", row.names=FALSE)

# Cough
metadata_np_cough <- subset(metadata_np_Pos, !is.na(cough))
ids_np_cough <- metadata_np_cough$alias_sequencing_id
counts_np_cough <- counts_np %>% dplyr::select(all_of(ids_np_cough))
nrow(metadata_np_cough)
all(rownames(metadata_np_cough) == colnames(counts_np_cough))
dds_np_cough <- DESeqDataSetFromMatrix(countData = counts_np_cough, colData = metadata_np_cough, design = ~ batch_num + sex + age_cat + cough + PC1 + PC2 + PC3 + PC4 + PC5)
# Require a normalized count of at least 10 in two or more samples
dds_np_cough_norm <- estimateSizeFactors(dds_np_cough)
dds_np_cough_norm <- counts(dds_np_cough_norm, normalized=TRUE)
filter <- rowSums(dds_np_cough_norm >= 10) >= (nrow(metadata_np_cough))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_np_cough <- dds_np_cough[filter,]
dds_np_cough <- DESeq(dds_np_cough)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_np_cough <- dds_np_cough[which(mcols(dds_np_cough)$betaConv),]
resultsNames(dds_np_cough) # lists the coefficients
dds_np_cough <- results(dds_np_cough, contrast=c("cough","Y","N"))
dds_np_cough <- data.frame(EnsemblID=row.names(dds_np_cough), dds_np_cough)
dds_np_cough <- merge(dds_np_cough, genes_np, by="row.names")
dds_np_cough <- dds_np_cough[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_np_cough <- subset(dds_np_cough, !is.na(padj))
dds_np_cough <- dds_np_cough[order(dds_np_cough$padj),]
dds_np_cough <- dds_np_cough %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_np_cough <- dds_np_cough[order(dds_np_cough$Expression, dds_np_cough$log2FoldChange),]
write.csv(dds_np_cough, "4_Symptoms/genes_np_cough.csv", row.names=FALSE)

# Rhinorrhea
metadata_np_rhinorrhea <- subset(metadata_np_Pos, !is.na(rhinorrhea))
ids_np_rhinorrhea <- metadata_np_rhinorrhea$alias_sequencing_id
counts_np_rhinorrhea <- counts_np %>% dplyr::select(all_of(ids_np_rhinorrhea))
nrow(metadata_np_rhinorrhea)
all(rownames(metadata_np_rhinorrhea) == colnames(counts_np_rhinorrhea))
dds_np_rhinorrhea <- DESeqDataSetFromMatrix(countData = counts_np_rhinorrhea, colData = metadata_np_rhinorrhea, design = ~ batch_num + sex + age_cat + rhinorrhea + PC1 + PC2 + PC3 + PC4 + PC5)
# Require a normalized count of at least 10 in two or more samples
dds_np_rhinorrhea_norm <- estimateSizeFactors(dds_np_rhinorrhea)
dds_np_rhinorrhea_norm <- counts(dds_np_rhinorrhea_norm, normalized=TRUE)
filter <- rowSums(dds_np_rhinorrhea_norm >= 10) >= (nrow(metadata_np_rhinorrhea))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_np_rhinorrhea <- dds_np_rhinorrhea[filter,]
dds_np_rhinorrhea <- DESeq(dds_np_rhinorrhea)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_np_rhinorrhea <- dds_np_rhinorrhea[which(mcols(dds_np_rhinorrhea)$betaConv),]
resultsNames(dds_np_rhinorrhea) # lists the coefficients
dds_np_rhinorrhea <- results(dds_np_rhinorrhea, contrast=c("rhinorrhea","Y","N"))
dds_np_rhinorrhea <- data.frame(EnsemblID=row.names(dds_np_rhinorrhea), dds_np_rhinorrhea)
dds_np_rhinorrhea <- merge(dds_np_rhinorrhea, genes_np, by="row.names")
dds_np_rhinorrhea <- dds_np_rhinorrhea[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_np_rhinorrhea <- subset(dds_np_rhinorrhea, !is.na(padj))
dds_np_rhinorrhea <- dds_np_rhinorrhea[order(dds_np_rhinorrhea$padj),]
dds_np_rhinorrhea <- dds_np_rhinorrhea %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_np_rhinorrhea <- dds_np_rhinorrhea[order(dds_np_rhinorrhea$Expression, dds_np_rhinorrhea$log2FoldChange),]
write.csv(dds_np_rhinorrhea, "4_Symptoms/genes_np_rhinorrhea.csv", row.names=FALSE)

# Nasal congestion
metadata_np_congestion <- subset(metadata_np_Pos, !is.na(congestion))
ids_np_congestion <- metadata_np_congestion$alias_sequencing_id
counts_np_congestion <- counts_np %>% dplyr::select(all_of(ids_np_congestion))
nrow(metadata_np_congestion)
all(rownames(metadata_np_congestion) == colnames(counts_np_congestion))
dds_np_congestion <- DESeqDataSetFromMatrix(countData = counts_np_congestion, colData = metadata_np_congestion, design = ~ batch_num + sex + age_cat + congestion + PC1 + PC2 + PC3 + PC4 + PC5)
# Require a normalized count of at least 10 in two or more samples
dds_np_congestion_norm <- estimateSizeFactors(dds_np_congestion)
dds_np_congestion_norm <- counts(dds_np_congestion_norm, normalized=TRUE)
filter <- rowSums(dds_np_congestion_norm >= 10) >= (nrow(metadata_np_congestion))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_np_congestion <- dds_np_congestion[filter,]
dds_np_congestion <- DESeq(dds_np_congestion)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_np_congestion <- dds_np_congestion[which(mcols(dds_np_congestion)$betaConv),]
resultsNames(dds_np_congestion) # lists the coefficients
dds_np_congestion <- results(dds_np_congestion, contrast=c("congestion","Y","N"))
dds_np_congestion <- data.frame(EnsemblID=row.names(dds_np_congestion), dds_np_congestion)
dds_np_congestion <- merge(dds_np_congestion, genes_np, by="row.names")
dds_np_congestion <- dds_np_congestion[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_np_congestion <- subset(dds_np_congestion, !is.na(padj))
dds_np_congestion <- dds_np_congestion[order(dds_np_congestion$padj),]
dds_np_congestion <- dds_np_congestion %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_np_congestion <- dds_np_congestion[order(dds_np_congestion$Expression, dds_np_congestion$log2FoldChange),]
write.csv(dds_np_congestion, "4_Symptoms/genes_np_congestion.csv", row.names=FALSE)

# Viral load
metadata_np_vl_copies <- subset(metadata_np_Pos, !is.na(vl_copies))
ids_np_vl_copies <- metadata_np_vl_copies$alias_sequencing_id
counts_np_vl_copies <- counts_np %>% dplyr::select(all_of(ids_np_vl_copies))
nrow(metadata_np_vl_copies)
all(rownames(metadata_np_vl_copies) == colnames(counts_np_vl_copies))
metadata_np_vl_copies$vl_copies <- log10(metadata_np_vl_copies$vl_copies)
median(metadata_np_vl_copies$vl_copies)
metadata_np_vl_copies$vl_high[metadata_np_vl_copies$vl_copies>=6] <- "Y"
metadata_np_vl_copies$vl_high[metadata_np_vl_copies$vl_copies<6] <- "N"
dds_np_vl_copies <- DESeqDataSetFromMatrix(countData = counts_np_vl_copies, colData = metadata_np_vl_copies, design = ~ batch_num + sex + age_cat + vl_high + PC1 + PC2 + PC3 + PC4 + PC5)
# Require a normalized count of at least 10 in two or more samples
dds_np_vl_copies_norm <- estimateSizeFactors(dds_np_vl_copies)
dds_np_vl_copies_norm <- counts(dds_np_vl_copies_norm, normalized=TRUE)
filter <- rowSums(dds_np_vl_copies_norm >= 10) >= (nrow(metadata_np_vl_copies))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_np_vl_copies <- dds_np_vl_copies[filter,]
dds_np_vl_copies <- DESeq(dds_np_vl_copies)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_np_vl_copies <- dds_np_vl_copies[which(mcols(dds_np_vl_copies)$betaConv),]
resultsNames(dds_np_vl_copies) # lists the coefficients
dds_np_vl_copies <- results(dds_np_vl_copies, contrast=c("vl_high","Y","N"))
dds_np_vl_copies <- data.frame(EnsemblID=row.names(dds_np_vl_copies), dds_np_vl_copies)
dds_np_vl_copies <- merge(dds_np_vl_copies, genes_np, by="row.names")
dds_np_vl_copies <- dds_np_vl_copies[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_np_vl_copies <- subset(dds_np_vl_copies, !is.na(padj))
dds_np_vl_copies <- dds_np_vl_copies[order(dds_np_vl_copies$padj),]
dds_np_vl_copies <- dds_np_vl_copies %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_np_vl_copies <- dds_np_vl_copies[order(dds_np_vl_copies$Expression, dds_np_vl_copies$log2FoldChange),]
write.csv(dds_np_vl_copies, "4_Symptoms/genes_np_vl_copies.csv", row.names=FALSE)

################################################################################
### MODELS FOR PATIENT CHARACTERISTICS & SYMPTOMS AMONG COVID+ - PAX SAMPLES ###
################################################################################

# Exclude adults from analyses of symptom data
metadata_pax_Pos_noadults <- subset(metadata_pax_Pos, age_cat!="Adult")

# Obesity
metadata_pax_obesity <- subset(metadata_pax_Pos_noadults, !is.na(obesity))
ids_pax_obesity <- metadata_pax_obesity$alias_sequencing_id
counts_pax_obesity <- counts_pax %>% dplyr::select(all_of(ids_pax_obesity))
nrow(metadata_pax_obesity)
all(rownames(metadata_pax_obesity) == colnames(counts_pax_obesity))
dds_pax_obesity <- DESeqDataSetFromMatrix(countData = counts_pax_obesity, colData = metadata_pax_obesity, design = ~ batch_num + sex + age_cat + obesity + PC1 + PC2 + PC3 + PC4 + PC5)
# Require a normalized count of at least 10 in two or more samples
dds_pax_obesity_norm <- estimateSizeFactors(dds_pax_obesity)
dds_pax_obesity_norm <- counts(dds_pax_obesity_norm, normalized=TRUE)
filter <- rowSums(dds_pax_obesity_norm >= 10) >= (nrow(metadata_pax_obesity))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_pax_obesity <- dds_pax_obesity[filter,]
dds_pax_obesity <- DESeq(dds_pax_obesity)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_pax_obesity <- dds_pax_obesity[which(mcols(dds_pax_obesity)$betaConv),]
resultsNames(dds_pax_obesity) # lists the coefficients
dds_pax_obesity <- results(dds_pax_obesity, contrast=c("obesity","Y","N"))
dds_pax_obesity <- data.frame(EnsemblID=row.names(dds_pax_obesity), dds_pax_obesity)
dds_pax_obesity <- merge(dds_pax_obesity, genes_pax, by="row.names")
dds_pax_obesity <- dds_pax_obesity[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_pax_obesity <- subset(dds_pax_obesity, !is.na(padj))
dds_pax_obesity <- dds_pax_obesity[order(dds_pax_obesity$padj),]
dds_pax_obesity <- dds_pax_obesity %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_pax_obesity <- dds_pax_obesity[order(dds_pax_obesity$Expression, dds_pax_obesity$log2FoldChange),]
write.csv(dds_pax_obesity, "4_Symptoms/genes_pax_obesity.csv", row.names=FALSE)

# Asthma
metadata_pax_asthma <- subset(metadata_pax_Pos_noadults, !is.na(asthma))
ids_pax_asthma <- metadata_pax_asthma$alias_sequencing_id
counts_pax_asthma <- counts_pax %>% dplyr::select(all_of(ids_pax_asthma))
nrow(metadata_pax_asthma)
all(rownames(metadata_pax_asthma) == colnames(counts_pax_asthma))
dds_pax_asthma <- DESeqDataSetFromMatrix(countData = counts_pax_asthma, colData = metadata_pax_asthma, design = ~ batch_num + sex + age_cat + asthma + PC1 + PC2 + PC3 + PC4 + PC5)
# Require a normalized count of at least 10 in two or more samples
dds_pax_asthma_norm <- estimateSizeFactors(dds_pax_asthma)
dds_pax_asthma_norm <- counts(dds_pax_asthma_norm, normalized=TRUE)
filter <- rowSums(dds_pax_asthma_norm >= 10) >= (nrow(metadata_pax_asthma))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_pax_asthma <- dds_pax_asthma[filter,]
dds_pax_asthma <- DESeq(dds_pax_asthma)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_pax_asthma <- dds_pax_asthma[which(mcols(dds_pax_asthma)$betaConv),]
resultsNames(dds_pax_asthma) # lists the coefficients
dds_pax_asthma <- results(dds_pax_asthma, contrast=c("asthma","Y","N"))
dds_pax_asthma <- data.frame(EnsemblID=row.names(dds_pax_asthma), dds_pax_asthma)
dds_pax_asthma <- merge(dds_pax_asthma, genes_pax, by="row.names")
dds_pax_asthma <- dds_pax_asthma[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_pax_asthma <- subset(dds_pax_asthma, !is.na(padj))
dds_pax_asthma <- dds_pax_asthma[order(dds_pax_asthma$padj),]
dds_pax_asthma <- dds_pax_asthma %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_pax_asthma <- dds_pax_asthma[order(dds_pax_asthma$Expression, dds_pax_asthma$log2FoldChange),]
write.csv(dds_pax_asthma, "4_Symptoms/genes_pax_asthma.csv", row.names=FALSE)

# Fever
metadata_pax_fever <- subset(metadata_pax_Pos_noadults, !is.na(fever))
ids_pax_fever <- metadata_pax_fever$alias_sequencing_id
counts_pax_fever <- counts_pax %>% dplyr::select(all_of(ids_pax_fever))
nrow(metadata_pax_fever)
all(rownames(metadata_pax_fever) == colnames(counts_pax_fever))
dds_pax_fever <- DESeqDataSetFromMatrix(countData = counts_pax_fever, colData = metadata_pax_fever, design = ~ batch_num + sex + age_cat + fever + PC1 + PC2 + PC3 + PC4 + PC5)
# Require a normalized count of at least 10 in two or more samples
dds_pax_fever_norm <- estimateSizeFactors(dds_pax_fever)
dds_pax_fever_norm <- counts(dds_pax_fever_norm, normalized=TRUE)
filter <- rowSums(dds_pax_fever_norm >= 10) >= (nrow(metadata_pax_fever))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_pax_fever <- dds_pax_fever[filter,]
dds_pax_fever <- DESeq(dds_pax_fever)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_pax_fever <- dds_pax_fever[which(mcols(dds_pax_fever)$betaConv),]
resultsNames(dds_pax_fever) # lists the coefficients
dds_pax_fever <- results(dds_pax_fever, contrast=c("fever","Y","N"))
dds_pax_fever <- data.frame(EnsemblID=row.names(dds_pax_fever), dds_pax_fever)
dds_pax_fever <- merge(dds_pax_fever, genes_pax, by="row.names")
dds_pax_fever <- dds_pax_fever[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_pax_fever <- subset(dds_pax_fever, !is.na(padj))
dds_pax_fever <- dds_pax_fever[order(dds_pax_fever$padj),]
dds_pax_fever <- dds_pax_fever %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_pax_fever <- dds_pax_fever[order(dds_pax_fever$Expression, dds_pax_fever$log2FoldChange),]
write.csv(dds_pax_fever, "4_Symptoms/genes_pax_fever.csv", row.names=FALSE)

# Cough
metadata_pax_cough <- subset(metadata_pax_Pos_noadults, !is.na(cough))
ids_pax_cough <- metadata_pax_cough$alias_sequencing_id
counts_pax_cough <- counts_pax %>% dplyr::select(all_of(ids_pax_cough))
nrow(metadata_pax_cough)
all(rownames(metadata_pax_cough) == colnames(counts_pax_cough))
dds_pax_cough <- DESeqDataSetFromMatrix(countData = counts_pax_cough, colData = metadata_pax_cough, design = ~ batch_num + sex + age_cat + cough + PC1 + PC2 + PC3 + PC4 + PC5)
# Require a normalized count of at least 10 in two or more samples
dds_pax_cough_norm <- estimateSizeFactors(dds_pax_cough)
dds_pax_cough_norm <- counts(dds_pax_cough_norm, normalized=TRUE)
filter <- rowSums(dds_pax_cough_norm >= 10) >= (nrow(metadata_pax_cough))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_pax_cough <- dds_pax_cough[filter,]
dds_pax_cough <- DESeq(dds_pax_cough)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_pax_cough <- dds_pax_cough[which(mcols(dds_pax_cough)$betaConv),]
resultsNames(dds_pax_cough) # lists the coefficients
dds_pax_cough <- results(dds_pax_cough, contrast=c("cough","Y","N"))
dds_pax_cough <- data.frame(EnsemblID=row.names(dds_pax_cough), dds_pax_cough)
dds_pax_cough <- merge(dds_pax_cough, genes_pax, by="row.names")
dds_pax_cough <- dds_pax_cough[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_pax_cough <- subset(dds_pax_cough, !is.na(padj))
dds_pax_cough <- dds_pax_cough[order(dds_pax_cough$padj),]
dds_pax_cough <- dds_pax_cough %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_pax_cough <- dds_pax_cough[order(dds_pax_cough$Expression, dds_pax_cough$log2FoldChange),]
write.csv(dds_pax_cough, "4_Symptoms/genes_pax_cough.csv", row.names=FALSE)

# Rhinorrhea
metadata_pax_rhinorrhea <- subset(metadata_pax_Pos_noadults, !is.na(rhinorrhea))
ids_pax_rhinorrhea <- metadata_pax_rhinorrhea$alias_sequencing_id
counts_pax_rhinorrhea <- counts_pax %>% dplyr::select(all_of(ids_pax_rhinorrhea))
nrow(metadata_pax_rhinorrhea)
all(rownames(metadata_pax_rhinorrhea) == colnames(counts_pax_rhinorrhea))
dds_pax_rhinorrhea <- DESeqDataSetFromMatrix(countData = counts_pax_rhinorrhea, colData = metadata_pax_rhinorrhea, design = ~ batch_num + sex + age_cat + rhinorrhea + PC1 + PC2 + PC3 + PC4 + PC5)
# Require a normalized count of at least 10 in two or more samples
dds_pax_rhinorrhea_norm <- estimateSizeFactors(dds_pax_rhinorrhea)
dds_pax_rhinorrhea_norm <- counts(dds_pax_rhinorrhea_norm, normalized=TRUE)
filter <- rowSums(dds_pax_rhinorrhea_norm >= 10) >= (nrow(metadata_pax_rhinorrhea))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_pax_rhinorrhea <- dds_pax_rhinorrhea[filter,]
dds_pax_rhinorrhea <- DESeq(dds_pax_rhinorrhea)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_pax_rhinorrhea <- dds_pax_rhinorrhea[which(mcols(dds_pax_rhinorrhea)$betaConv),]
resultsNames(dds_pax_rhinorrhea) # lists the coefficients
dds_pax_rhinorrhea <- results(dds_pax_rhinorrhea, contrast=c("rhinorrhea","Y","N"))
dds_pax_rhinorrhea <- data.frame(EnsemblID=row.names(dds_pax_rhinorrhea), dds_pax_rhinorrhea)
dds_pax_rhinorrhea <- merge(dds_pax_rhinorrhea, genes_pax, by="row.names")
dds_pax_rhinorrhea <- dds_pax_rhinorrhea[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_pax_rhinorrhea <- subset(dds_pax_rhinorrhea, !is.na(padj))
dds_pax_rhinorrhea <- dds_pax_rhinorrhea[order(dds_pax_rhinorrhea$padj),]
dds_pax_rhinorrhea <- dds_pax_rhinorrhea %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_pax_rhinorrhea <- dds_pax_rhinorrhea[order(dds_pax_rhinorrhea$Expression, dds_pax_rhinorrhea$log2FoldChange),]
write.csv(dds_pax_rhinorrhea, "4_Symptoms/genes_pax_rhinorrhea.csv", row.names=FALSE)

# Nasal congestion
metadata_pax_congestion <- subset(metadata_pax_Pos_noadults, !is.na(congestion))
ids_pax_congestion <- metadata_pax_congestion$alias_sequencing_id
counts_pax_congestion <- counts_pax %>% dplyr::select(all_of(ids_pax_congestion))
nrow(metadata_pax_congestion)
all(rownames(metadata_pax_congestion) == colnames(counts_pax_congestion))
dds_pax_congestion <- DESeqDataSetFromMatrix(countData = counts_pax_congestion, colData = metadata_pax_congestion, design = ~ batch_num + sex + age_cat + congestion + PC1 + PC2 + PC3 + PC4 + PC5)
# Require a normalized count of at least 10 in two or more samples
dds_pax_congestion_norm <- estimateSizeFactors(dds_pax_congestion)
dds_pax_congestion_norm <- counts(dds_pax_congestion_norm, normalized=TRUE)
filter <- rowSums(dds_pax_congestion_norm >= 10) >= (nrow(metadata_pax_congestion))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_pax_congestion <- dds_pax_congestion[filter,]
dds_pax_congestion <- DESeq(dds_pax_congestion)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_pax_congestion <- dds_pax_congestion[which(mcols(dds_pax_congestion)$betaConv),]
resultsNames(dds_pax_congestion) # lists the coefficients
dds_pax_congestion <- results(dds_pax_congestion, contrast=c("congestion","Y","N"))
dds_pax_congestion <- data.frame(EnsemblID=row.names(dds_pax_congestion), dds_pax_congestion)
dds_pax_congestion <- merge(dds_pax_congestion, genes_pax, by="row.names")
dds_pax_congestion <- dds_pax_congestion[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_pax_congestion <- subset(dds_pax_congestion, !is.na(padj))
dds_pax_congestion <- dds_pax_congestion[order(dds_pax_congestion$padj),]
dds_pax_congestion <- dds_pax_congestion %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_pax_congestion <- dds_pax_congestion[order(dds_pax_congestion$Expression, dds_pax_congestion$log2FoldChange),]
write.csv(dds_pax_congestion, "4_Symptoms/genes_pax_congestion.csv", row.names=FALSE)

# Viral load - high (>=10^6 copies/mL vs. <10^6 copies/mL)
metadata_pax_vl_copies <- subset(metadata_pax_Pos_noadults, !is.na(vl_copies))
ids_pax_vl_copies <- metadata_pax_vl_copies$alias_sequencing_id
counts_pax_vl_copies <- counts_pax %>% dplyr::select(all_of(ids_pax_vl_copies))
nrow(metadata_pax_vl_copies)
all(rownames(metadata_pax_vl_copies) == colnames(counts_pax_vl_copies))
metadata_pax_vl_copies$vl_copies <- log10(metadata_pax_vl_copies$vl_copies)
median(metadata_pax_vl_copies$vl_copies)
metadata_pax_vl_copies$vl_high[metadata_pax_vl_copies$vl_copies>=6] <- "Y"
metadata_pax_vl_copies$vl_high[metadata_pax_vl_copies$vl_copies<6] <- "N"
dds_pax_vl_copies <- DESeqDataSetFromMatrix(countData = counts_pax_vl_copies, colData = metadata_pax_vl_copies, design = ~ batch_num + sex + age_cat + vl_high + PC1 + PC2 + PC3 + PC4 + PC5)
# Require a normalized count of at least 10 in two or more samples
dds_pax_vl_copies_norm <- estimateSizeFactors(dds_pax_vl_copies)
dds_pax_vl_copies_norm <- counts(dds_pax_vl_copies_norm, normalized=TRUE)
filter <- rowSums(dds_pax_vl_copies_norm >= 10) >= (nrow(metadata_pax_vl_copies))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_pax_vl_copies <- dds_pax_vl_copies[filter,]
dds_pax_vl_copies <- DESeq(dds_pax_vl_copies)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_pax_vl_copies <- dds_pax_vl_copies[which(mcols(dds_pax_vl_copies)$betaConv),]
resultsNames(dds_pax_vl_copies) # lists the coefficients
dds_pax_vl_copies <- results(dds_pax_vl_copies, contrast=c("vl_high","Y","N"))
dds_pax_vl_copies <- data.frame(EnsemblID=row.names(dds_pax_vl_copies), dds_pax_vl_copies)
dds_pax_vl_copies <- merge(dds_pax_vl_copies, genes_pax, by="row.names")
dds_pax_vl_copies <- dds_pax_vl_copies[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_pax_vl_copies <- subset(dds_pax_vl_copies, !is.na(padj))
dds_pax_vl_copies <- dds_pax_vl_copies[order(dds_pax_vl_copies$padj),]
dds_pax_vl_copies <- dds_pax_vl_copies %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_pax_vl_copies <- dds_pax_vl_copies[order(dds_pax_vl_copies$Expression, dds_pax_vl_copies$log2FoldChange),]
write.csv(dds_pax_vl_copies, "4_Symptoms/genes_pax_vl_copies.csv", row.names=FALSE)