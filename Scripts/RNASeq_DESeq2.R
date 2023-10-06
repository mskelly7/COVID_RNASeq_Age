# BRAVE Kids RNA Sequencing Analysis
# Aditya Mohan (MD/PhD candidate) / Matthew Kelly, MD, MPH 
# DESeq2 models for differentially expressed genes
# Models are adjusted for sequencing batch, sex, cell populations (PC1-PC5), and (as applicable) age category
# Last update: Oct. 5, 2023

remove(list=ls())
setwd("____________________________") 
options(rstudio.help.showDataPreview = FALSE)
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
phy.rnaseq.np <- readRDS("phy.rnaseq.np.rds")
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

phy.rnaseq.pax <- readRDS("phy.rnaseq.pax.rds")
metadata_pax <- data.frame(sample_data(phy.rnaseq.pax))
metadata_pax$batch_num <- as.factor(metadata_pax$batch_num)
metadata_pax$age_cat <- as.factor(metadata_pax$age_cat)
metadata_pax$age_cat2[metadata_pax$age_cat=="0-5 years" | metadata_pax$age_cat=="6-13 years" | metadata_pax$age_cat=="14-20 years"] <- "Child"
metadata_pax$age_cat2[metadata_pax$age_cat=="Adult"] <- "Adult"
metadata_pax$age_cat2 <- as.factor(metadata_pax$age_cat2)
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

##################################################################
### DIFFERENCES IN NP GENE EXPRESSION BY AGE IN COVID NEGATIVE ###
##################################################################

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
write.csv(dds_np_Neg_0to5_6to13, "Statistical_Analyses/1_COVID_Neg_by_Age/genes_np_Neg_0to5_6to13.csv", row.names=FALSE)

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
write.csv(dds_np_Neg_0to5_14to20, "Statistical_Analyses/1_COVID_Neg_by_Age/genes_np_Neg_0to5_14to20.csv", row.names=FALSE)

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
write.csv(dds_np_Neg_6to13_14to20, "Statistical_Analyses/1_COVID_Neg_by_Age/genes_np_Neg_6to13_14to20.csv", row.names=FALSE)

###################################################################
### DIFFERENCES IN PAX GENE EXPRESSION BY AGE IN COVID NEGATIVE ###
###################################################################

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
# results(dds, contrast=c("condition","C","B")) = genes with logFC > 0 are overexpressed in C
dds_pax_Neg_0to5_6to13 <- results(dds_pax_Neg, contrast=c("age_cat","0-5 years","6-13 years"))
dds_pax_Neg_0to5_6to13 <- data.frame(EnsemblID=row.names(dds_pax_Neg_0to5_6to13), dds_pax_Neg_0to5_6to13)
dds_pax_Neg_0to5_6to13 <- merge(dds_pax_Neg_0to5_6to13, genes_pax, by="row.names")
dds_pax_Neg_0to5_6to13 <- dds_pax_Neg_0to5_6to13[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_pax_Neg_0to5_6to13 <- subset(dds_pax_Neg_0to5_6to13, !is.na(padj))
dds_pax_Neg_0to5_6to13 <- dds_pax_Neg_0to5_6to13[order(dds_pax_Neg_0to5_6to13$padj),]
dds_pax_Neg_0to5_6to13 <- dds_pax_Neg_0to5_6to13 %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_pax_Neg_0to5_6to13 <- dds_pax_Neg_0to5_6to13[order(dds_pax_Neg_0to5_6to13$Expression, dds_pax_Neg_0to5_6to13$log2FoldChange),]
write.csv(dds_pax_Neg_0to5_6to13, "Statistical_Analyses/1_COVID_Neg_by_Age/genes_pax_Neg_0to5_6to13.csv", row.names=FALSE)

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
write.csv(dds_pax_Neg_0to5_14to20, "Statistical_Analyses/1_COVID_Neg_by_Age/genes_pax_Neg_0to5_14to20.csv", row.names=FALSE)

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
write.csv(dds_pax_Neg_6to13_14to20, "Statistical_Analyses/1_COVID_Neg_by_Age/genes_pax_Neg_6to13_14to20.csv", row.names=FALSE)

# Output results for 0-5 years vs. Adult = genes with logFC>0 are overexpressed in children 0-5 years
# results(dds, contrast=c("condition","C","B")) = genes with logFC > 0 are overexpressed in C
dds_pax_Neg_0to5_Adult <- results(dds_pax_Neg, contrast=c("age_cat","0-5 years","Adult"))
dds_pax_Neg_0to5_Adult <- data.frame(EnsemblID=row.names(dds_pax_Neg_0to5_Adult), dds_pax_Neg_0to5_Adult)
dds_pax_Neg_0to5_Adult <- merge(dds_pax_Neg_0to5_Adult, genes_pax, by="row.names")
dds_pax_Neg_0to5_Adult <- dds_pax_Neg_0to5_Adult[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_pax_Neg_0to5_Adult <- subset(dds_pax_Neg_0to5_Adult, !is.na(padj))
dds_pax_Neg_0to5_Adult <- dds_pax_Neg_0to5_Adult[order(dds_pax_Neg_0to5_Adult$padj),]
dds_pax_Neg_0to5_Adult <- dds_pax_Neg_0to5_Adult %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_pax_Neg_0to5_Adult <- dds_pax_Neg_0to5_Adult[order(dds_pax_Neg_0to5_Adult$Expression, dds_pax_Neg_0to5_Adult$log2FoldChange),]
write.csv(dds_pax_Neg_0to5_Adult, "Statistical_Analyses/1_COVID_Neg_by_Age/genes_pax_Neg_0to5_Adult.csv", row.names=FALSE)

# Output results for 6-13 years vs. Adult = genes with logFC>0 are overexpressed in children 6-13 years
# results(dds, contrast=c("condition","C","B")) = genes with logFC > 0 are overexpressed in C
dds_pax_Neg_6to13_Adult <- results(dds_pax_Neg, contrast=c("age_cat","6-13 years","Adult"))
dds_pax_Neg_6to13_Adult <- data.frame(EnsemblID=row.names(dds_pax_Neg_6to13_Adult), dds_pax_Neg_6to13_Adult)
dds_pax_Neg_6to13_Adult <- merge(dds_pax_Neg_6to13_Adult, genes_pax, by="row.names")
dds_pax_Neg_6to13_Adult <- dds_pax_Neg_6to13_Adult[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_pax_Neg_6to13_Adult <- subset(dds_pax_Neg_6to13_Adult, !is.na(padj))
dds_pax_Neg_6to13_Adult <- dds_pax_Neg_6to13_Adult[order(dds_pax_Neg_6to13_Adult$padj),]
dds_pax_Neg_6to13_Adult <- dds_pax_Neg_6to13_Adult %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_pax_Neg_6to13_Adult <- dds_pax_Neg_6to13_Adult[order(dds_pax_Neg_6to13_Adult$Expression, dds_pax_Neg_6to13_Adult$log2FoldChange),]
write.csv(dds_pax_Neg_6to13_Adult, "Statistical_Analyses/1_COVID_Neg_by_Age/genes_pax_Neg_6to13_Adult.csv", row.names=FALSE)

# Output results for 14-20 years vs. Adult = genes with logFC>0 are overexpressed in children 14-20 years
# results(dds, contrast=c("condition","C","B")) = genes with logFC > 0 are overexpressed in C
dds_pax_Neg_14to20_Adult <- results(dds_pax_Neg, contrast=c("age_cat","14-20 years","Adult"))
dds_pax_Neg_14to20_Adult <- data.frame(EnsemblID=row.names(dds_pax_Neg_14to20_Adult), dds_pax_Neg_14to20_Adult)
dds_pax_Neg_14to20_Adult <- merge(dds_pax_Neg_14to20_Adult, genes_pax, by="row.names")
dds_pax_Neg_14to20_Adult <- dds_pax_Neg_14to20_Adult[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_pax_Neg_14to20_Adult <- subset(dds_pax_Neg_14to20_Adult, !is.na(padj))
dds_pax_Neg_14to20_Adult <- dds_pax_Neg_14to20_Adult[order(dds_pax_Neg_14to20_Adult$padj),]
dds_pax_Neg_14to20_Adult <- dds_pax_Neg_14to20_Adult %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_pax_Neg_14to20_Adult <- dds_pax_Neg_14to20_Adult[order(dds_pax_Neg_14to20_Adult$Expression, dds_pax_Neg_14to20_Adult$log2FoldChange),]
write.csv(dds_pax_Neg_14to20_Adult, "Statistical_Analyses/1_COVID_Neg_by_Age/genes_pax_Neg_14to20_Adult.csv", row.names=FALSE)

# Compare gene expression among CHILDREN and ADULTS

nrow(metadata_pax_Neg)
all(rownames(metadata_pax_Neg) == colnames(counts_pax_Neg))
dds_pax_Neg_agecat2 <- DESeqDataSetFromMatrix(countData = counts_pax_Neg, colData = metadata_pax_Neg, design = ~ batch_num + age_cat2 + sex + PC1 + PC2 + PC3 + PC4 + PC5)

# Require a normalized count of at least 10 in two or more samples
dds_pax_Neg_agecat2_norm <- estimateSizeFactors(dds_pax_Neg_agecat2)
dds_pax_Neg_agecat2_norm <- counts(dds_pax_Neg_agecat2_norm, normalized=TRUE)
filter <- rowSums(dds_pax_Neg_agecat2_norm >= 10) >= (nrow(metadata_pax_Neg))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_pax_Neg_agecat2 <- dds_pax_Neg_agecat2[filter,]
dds_pax_Neg_agecat2 <- DESeq(dds_pax_Neg_agecat2)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_pax_Neg_agecat2 <- dds_pax_Neg_agecat2[which(mcols(dds_pax_Neg_agecat2)$betaConv),]
resultsNames(dds_pax_Neg_agecat2) # lists the coefficients

# Output results for Children vs. Adult = genes with logFC>0 are overexpressed in children 0-5 years
# results(dds, contrast=c("condition","C","B")) = genes with logFC > 0 are overexpressed in C
dds_pax_Neg_Child_Adult <- results(dds_pax_Neg_agecat2, contrast=c("age_cat2","Child","Adult"))
dds_pax_Neg_Child_Adult <- data.frame(EnsemblID=row.names(dds_pax_Neg_Child_Adult), dds_pax_Neg_Child_Adult)
dds_pax_Neg_Child_Adult <- merge(dds_pax_Neg_Child_Adult, genes_pax, by="row.names")
dds_pax_Neg_Child_Adult <- dds_pax_Neg_Child_Adult[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_pax_Neg_Child_Adult <- subset(dds_pax_Neg_Child_Adult, !is.na(padj))
dds_pax_Neg_Child_Adult <- dds_pax_Neg_Child_Adult[order(dds_pax_Neg_Child_Adult$padj),]
dds_pax_Neg_Child_Adult <- dds_pax_Neg_Child_Adult %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_pax_Neg_Child_Adult <- dds_pax_Neg_Child_Adult[order(dds_pax_Neg_Child_Adult$Expression, dds_pax_Neg_Child_Adult$log2FoldChange),]
write.csv(dds_pax_Neg_Child_Adult, "Statistical_Analyses/1_COVID_Neg_by_Age/genes_pax_Neg_Child_Adult.csv", row.names=FALSE)

####################################################
### DIFFERENCES IN NP EXPRESSION BY COVID STATUS ###
####################################################

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
write.csv(dds_np_Pos_Neg, "Statistical_Analyses/2_COVID_Pos_vs_Neg/genes_np_Pos_Neg.csv", row.names=FALSE)

##########################################################
### DIFFERENCES IN PAX GENE EXPRESSION BY COVID STATUS ###
##########################################################

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
write.csv(dds_pax_Pos_Neg, "Statistical_Analyses/2_COVID_Pos_vs_Neg/genes_pax_Pos_Neg.csv", row.names=FALSE)

##############################################################
### COMPARING NP SAMPLES BY COVID STATUS IN AGE CATEGORIES ###
##############################################################

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
write.csv(dds_np_0to5_Pos_Neg, "Statistical_Analyses/3_COVID_Pos_by_Age/genes_np_0to5_Pos_Neg.csv", row.names=FALSE)

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
write.csv(dds_np_6to13_Pos_Neg, "Statistical_Analyses/3_COVID_Pos_by_Age/genes_np_6to13_Pos_Neg.csv", row.names=FALSE)

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
write.csv(dds_np_14to20_Pos_Neg, "Statistical_Analyses/3_COVID_Pos_by_Age/genes_np_14to20_Pos_Neg.csv", row.names=FALSE)

###############################################################
### COMPARING PAX SAMPLES BY COVID STATUS IN AGE CATEGORIES ###
###############################################################

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
write.csv(dds_pax_0to5_Pos_Neg, "Statistical_Analyses/3_COVID_Pos_by_Age/genes_pax_0to5_Pos_Neg.csv", row.names=FALSE)

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
write.csv(dds_pax_6to13_Pos_Neg, "Statistical_Analyses/3_COVID_Pos_by_Age/genes_pax_6to13_Pos_Neg.csv", row.names=FALSE)

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
write.csv(dds_pax_14to20_Pos_Neg, "Statistical_Analyses/3_COVID_Pos_by_Age/genes_pax_14to20_Pos_Neg.csv", row.names=FALSE)

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
write.csv(dds_pax_adult_Pos_Neg, "Statistical_Analyses/3_COVID_Pos_by_Age/genes_pax_adult_Pos_Neg.csv", row.names=FALSE)

####################################################################
### MODELS FOR ILLNESS CHARACTERISTICS AMONG COVID+ - NP SAMPLES ###
####################################################################

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
write.csv(dds_np_fever, "Statistical_Analyses/4_Symptoms/genes_np_fever.csv", row.names=FALSE)

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
write.csv(dds_np_cough, "Statistical_Analyses/4_Symptoms/genes_np_cough.csv", row.names=FALSE)

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
write.csv(dds_np_rhinorrhea, "Statistical_Analyses/4_Symptoms/genes_np_rhinorrhea.csv", row.names=FALSE)

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
write.csv(dds_np_congestion, "Statistical_Analyses/4_Symptoms/genes_np_congestion.csv", row.names=FALSE)

# Headache
metadata_np_headache <- subset(metadata_np_Pos, !is.na(headache))
ids_np_headache <- metadata_np_headache$alias_sequencing_id
counts_np_headache <- counts_np %>% dplyr::select(all_of(ids_np_headache))
nrow(metadata_np_headache)
all(rownames(metadata_np_headache) == colnames(counts_np_headache))
dds_np_headache <- DESeqDataSetFromMatrix(countData = counts_np_headache, colData = metadata_np_headache, design = ~ batch_num + sex + age_cat + headache + PC1 + PC2 + PC3 + PC4 + PC5)
# Require a normalized count of at least 10 in two or more samples
dds_np_headache_norm <- estimateSizeFactors(dds_np_headache)
dds_np_headache_norm <- counts(dds_np_headache_norm, normalized=TRUE)
filter <- rowSums(dds_np_headache_norm >= 10) >= (nrow(metadata_np_headache))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_np_headache <- dds_np_headache[filter,]
dds_np_headache <- DESeq(dds_np_headache)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_np_headache <- dds_np_headache[which(mcols(dds_np_headache)$betaConv),]
resultsNames(dds_np_headache) # lists the coefficients
dds_np_headache <- results(dds_np_headache, contrast=c("headache","Y","N"))
dds_np_headache <- data.frame(EnsemblID=row.names(dds_np_headache), dds_np_headache)
dds_np_headache <- merge(dds_np_headache, genes_np, by="row.names")
dds_np_headache <- dds_np_headache[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_np_headache <- subset(dds_np_headache, !is.na(padj))
dds_np_headache <- dds_np_headache[order(dds_np_headache$padj),]
dds_np_headache <- dds_np_headache %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_np_headache <- dds_np_headache[order(dds_np_headache$Expression, dds_np_headache$log2FoldChange),]
write.csv(dds_np_headache, "Statistical_Analyses/4_Symptoms/genes_np_headache.csv", row.names=FALSE)

# Abdominal pain
metadata_np_abd_pain <- subset(metadata_np_Pos, !is.na(abd_pain))
ids_np_abd_pain <- metadata_np_abd_pain$alias_sequencing_id
counts_np_abd_pain <- counts_np %>% dplyr::select(all_of(ids_np_abd_pain))
nrow(metadata_np_abd_pain)
all(rownames(metadata_np_abd_pain) == colnames(counts_np_abd_pain))
dds_np_abd_pain <- DESeqDataSetFromMatrix(countData = counts_np_abd_pain, colData = metadata_np_abd_pain, design = ~ batch_num + sex + age_cat + abd_pain + PC1 + PC2 + PC3 + PC4 + PC5)
# Require a normalized count of at least 10 in two or more samples
dds_np_abd_pain_norm <- estimateSizeFactors(dds_np_abd_pain)
dds_np_abd_pain_norm <- counts(dds_np_abd_pain_norm, normalized=TRUE)
filter <- rowSums(dds_np_abd_pain_norm >= 10) >= (nrow(metadata_np_abd_pain))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_np_abd_pain <- dds_np_abd_pain[filter,]
dds_np_abd_pain <- DESeq(dds_np_abd_pain)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_np_abd_pain <- dds_np_abd_pain[which(mcols(dds_np_abd_pain)$betaConv),]
resultsNames(dds_np_abd_pain) # lists the coefficients
dds_np_abd_pain <- results(dds_np_abd_pain, contrast=c("abd_pain","Y","N"))
dds_np_abd_pain <- data.frame(EnsemblID=row.names(dds_np_abd_pain), dds_np_abd_pain)
dds_np_abd_pain <- merge(dds_np_abd_pain, genes_np, by="row.names")
dds_np_abd_pain <- dds_np_abd_pain[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_np_abd_pain <- subset(dds_np_abd_pain, !is.na(padj))
dds_np_abd_pain <- dds_np_abd_pain[order(dds_np_abd_pain$padj),]
dds_np_abd_pain <- dds_np_abd_pain %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_np_abd_pain <- dds_np_abd_pain[order(dds_np_abd_pain$Expression, dds_np_abd_pain$log2FoldChange),]
write.csv(dds_np_abd_pain, "Statistical_Analyses/4_Symptoms/genes_np_abd_pain.csv", row.names=FALSE)

# Anosmia
metadata_np_anosmia <- subset(metadata_np_Pos, !is.na(anosmia))
ids_np_anosmia <- metadata_np_anosmia$alias_sequencing_id
counts_np_anosmia <- counts_np %>% dplyr::select(all_of(ids_np_anosmia))
nrow(metadata_np_anosmia)
all(rownames(metadata_np_anosmia) == colnames(counts_np_anosmia))
dds_np_anosmia <- DESeqDataSetFromMatrix(countData = counts_np_anosmia, colData = metadata_np_anosmia, design = ~ batch_num + sex + age_cat + anosmia + PC1 + PC2 + PC3 + PC4 + PC5)
# Require a normalized count of at least 10 in two or more samples
dds_np_anosmia_norm <- estimateSizeFactors(dds_np_anosmia)
dds_np_anosmia_norm <- counts(dds_np_anosmia_norm, normalized=TRUE)
filter <- rowSums(dds_np_anosmia_norm >= 10) >= (nrow(metadata_np_anosmia))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_np_anosmia <- dds_np_anosmia[filter,]
dds_np_anosmia <- DESeq(dds_np_anosmia)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_np_anosmia <- dds_np_anosmia[which(mcols(dds_np_anosmia)$betaConv),]
resultsNames(dds_np_anosmia) # lists the coefficients
dds_np_anosmia <- results(dds_np_anosmia, contrast=c("anosmia","Y","N"))
dds_np_anosmia <- data.frame(EnsemblID=row.names(dds_np_anosmia), dds_np_anosmia)
dds_np_anosmia <- merge(dds_np_anosmia, genes_np, by="row.names")
dds_np_anosmia <- dds_np_anosmia[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_np_anosmia <- subset(dds_np_anosmia, !is.na(padj))
dds_np_anosmia <- dds_np_anosmia[order(dds_np_anosmia$padj),]
dds_np_anosmia <- dds_np_anosmia %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_np_anosmia <- dds_np_anosmia[order(dds_np_anosmia$Expression, dds_np_anosmia$log2FoldChange),]
write.csv(dds_np_anosmia, "Statistical_Analyses/4_Symptoms/genes_np_anosmia.csv", row.names=FALSE)

# Dysgeusia
metadata_np_dysgeusia <- subset(metadata_np_Pos, !is.na(dysgeusia))
ids_np_dysgeusia <- metadata_np_dysgeusia$alias_sequencing_id
counts_np_dysgeusia <- counts_np %>% dplyr::select(all_of(ids_np_dysgeusia))
nrow(metadata_np_dysgeusia)
all(rownames(metadata_np_dysgeusia) == colnames(counts_np_dysgeusia))
dds_np_dysgeusia <- DESeqDataSetFromMatrix(countData = counts_np_dysgeusia, colData = metadata_np_dysgeusia, design = ~ batch_num + sex + age_cat + dysgeusia + PC1 + PC2 + PC3 + PC4 + PC5)
# Require a normalized count of at least 10 in two or more samples
dds_np_dysgeusia_norm <- estimateSizeFactors(dds_np_dysgeusia)
dds_np_dysgeusia_norm <- counts(dds_np_dysgeusia_norm, normalized=TRUE)
filter <- rowSums(dds_np_dysgeusia_norm >= 10) >= (nrow(metadata_np_dysgeusia))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_np_dysgeusia <- dds_np_dysgeusia[filter,]
dds_np_dysgeusia <- DESeq(dds_np_dysgeusia)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_np_dysgeusia <- dds_np_dysgeusia[which(mcols(dds_np_dysgeusia)$betaConv),]
resultsNames(dds_np_dysgeusia) # lists the coefficients
dds_np_dysgeusia <- results(dds_np_dysgeusia, contrast=c("dysgeusia","Y","N"))
dds_np_dysgeusia <- data.frame(EnsemblID=row.names(dds_np_dysgeusia), dds_np_dysgeusia)
dds_np_dysgeusia <- merge(dds_np_dysgeusia, genes_np, by="row.names")
dds_np_dysgeusia <- dds_np_dysgeusia[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_np_dysgeusia <- subset(dds_np_dysgeusia, !is.na(padj))
dds_np_dysgeusia <- dds_np_dysgeusia[order(dds_np_dysgeusia$padj),]
dds_np_dysgeusia <- dds_np_dysgeusia %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_np_dysgeusia <- dds_np_dysgeusia[order(dds_np_dysgeusia$Expression, dds_np_dysgeusia$log2FoldChange),]
write.csv(dds_np_dysgeusia, "Statistical_Analyses/4_Symptoms/genes_np_dysgeusia.csv", row.names=FALSE)

# Myalgias
metadata_np_myalgias <- subset(metadata_np_Pos, !is.na(myalgias))
ids_np_myalgias <- metadata_np_myalgias$alias_sequencing_id
counts_np_myalgias <- counts_np %>% dplyr::select(all_of(ids_np_myalgias))
nrow(metadata_np_myalgias)
all(rownames(metadata_np_myalgias) == colnames(counts_np_myalgias))
dds_np_myalgias <- DESeqDataSetFromMatrix(countData = counts_np_myalgias, colData = metadata_np_myalgias, design = ~ batch_num + sex + age_cat + myalgias + PC1 + PC2 + PC3 + PC4 + PC5)
# Require a normalized count of at least 10 in two or more samples
dds_np_myalgias_norm <- estimateSizeFactors(dds_np_myalgias)
dds_np_myalgias_norm <- counts(dds_np_myalgias_norm, normalized=TRUE)
filter <- rowSums(dds_np_myalgias_norm >= 10) >= (nrow(metadata_np_myalgias))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_np_myalgias <- dds_np_myalgias[filter,]
dds_np_myalgias <- DESeq(dds_np_myalgias)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_np_myalgias <- dds_np_myalgias[which(mcols(dds_np_myalgias)$betaConv),]
resultsNames(dds_np_myalgias) # lists the coefficients
dds_np_myalgias <- results(dds_np_myalgias, contrast=c("myalgias","Y","N"))
dds_np_myalgias <- data.frame(EnsemblID=row.names(dds_np_myalgias), dds_np_myalgias)
dds_np_myalgias <- merge(dds_np_myalgias, genes_np, by="row.names")
dds_np_myalgias <- dds_np_myalgias[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_np_myalgias <- subset(dds_np_myalgias, !is.na(padj))
dds_np_myalgias <- dds_np_myalgias[order(dds_np_myalgias$padj),]
dds_np_myalgias <- dds_np_myalgias %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_np_myalgias <- dds_np_myalgias[order(dds_np_myalgias$Expression, dds_np_myalgias$log2FoldChange),]
write.csv(dds_np_myalgias, "Statistical_Analyses/4_Symptoms/genes_np_myalgias.csv", row.names=FALSE)

################################################################################
### MODELS FOR PATIENT CHARACTERISTICS & SYMPTOMS AMONG COVID+ - PAX SAMPLES ###
################################################################################

# Fever
metadata_pax_fever <- subset(metadata_pax_Pos, !is.na(fever))
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
write.csv(dds_pax_fever, "Statistical_Analyses/4_Symptoms/genes_pax_fever.csv", row.names=FALSE)

# Cough
metadata_pax_cough <- subset(metadata_pax_Pos, !is.na(cough))
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
write.csv(dds_pax_cough, "Statistical_Analyses/4_Symptoms/genes_pax_cough.csv", row.names=FALSE)

# Rhinorrhea
metadata_pax_rhinorrhea <- subset(metadata_pax_Pos, !is.na(rhinorrhea))
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
write.csv(dds_pax_rhinorrhea, "Statistical_Analyses/4_Symptoms/genes_pax_rhinorrhea.csv", row.names=FALSE)

# Nasal congestion
metadata_pax_congestion <- subset(metadata_pax_Pos, !is.na(congestion))
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
write.csv(dds_pax_congestion, "Statistical_Analyses/4_Symptoms/genes_pax_congestion.csv", row.names=FALSE)

# Headache
metadata_pax_headache <- subset(metadata_pax_Pos, !is.na(headache))
ids_pax_headache <- metadata_pax_headache$alias_sequencing_id
counts_pax_headache <- counts_pax %>% dplyr::select(all_of(ids_pax_headache))
nrow(metadata_pax_headache)
all(rownames(metadata_pax_headache) == colnames(counts_pax_headache))
dds_pax_headache <- DESeqDataSetFromMatrix(countData = counts_pax_headache, colData = metadata_pax_headache, design = ~ batch_num + sex + age_cat + headache + PC1 + PC2 + PC3 + PC4 + PC5)
# Require a normalized count of at least 10 in two or more samples
dds_pax_headache_norm <- estimateSizeFactors(dds_pax_headache)
dds_pax_headache_norm <- counts(dds_pax_headache_norm, normalized=TRUE)
filter <- rowSums(dds_pax_headache_norm >= 10) >= (nrow(metadata_pax_headache))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_pax_headache <- dds_pax_headache[filter,]
dds_pax_headache <- DESeq(dds_pax_headache)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_pax_headache <- dds_pax_headache[which(mcols(dds_pax_headache)$betaConv),]
resultsNames(dds_pax_headache) # lists the coefficients
dds_pax_headache <- results(dds_pax_headache, contrast=c("headache","Y","N"))
dds_pax_headache <- data.frame(EnsemblID=row.names(dds_pax_headache), dds_pax_headache)
dds_pax_headache <- merge(dds_pax_headache, genes_pax, by="row.names")
dds_pax_headache <- dds_pax_headache[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_pax_headache <- subset(dds_pax_headache, !is.na(padj))
dds_pax_headache <- dds_pax_headache[order(dds_pax_headache$padj),]
dds_pax_headache <- dds_pax_headache %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_pax_headache <- dds_pax_headache[order(dds_pax_headache$Expression, dds_pax_headache$log2FoldChange),]
write.csv(dds_pax_headache, "Statistical_Analyses/4_Symptoms/genes_pax_headache.csv", row.names=FALSE)

# Abdominal pain
metadata_pax_abd_pain <- subset(metadata_pax_Pos, !is.na(abd_pain))
ids_pax_abd_pain <- metadata_pax_abd_pain$alias_sequencing_id
counts_pax_abd_pain <- counts_pax %>% dplyr::select(all_of(ids_pax_abd_pain))
nrow(metadata_pax_abd_pain)
all(rownames(metadata_pax_abd_pain) == colnames(counts_pax_abd_pain))
dds_pax_abd_pain <- DESeqDataSetFromMatrix(countData = counts_pax_abd_pain, colData = metadata_pax_abd_pain, design = ~ batch_num + sex + age_cat + abd_pain + PC1 + PC2 + PC3 + PC4 + PC5)
# Require a normalized count of at least 10 in two or more samples
dds_pax_abd_pain_norm <- estimateSizeFactors(dds_pax_abd_pain)
dds_pax_abd_pain_norm <- counts(dds_pax_abd_pain_norm, normalized=TRUE)
filter <- rowSums(dds_pax_abd_pain_norm >= 10) >= (nrow(metadata_pax_abd_pain))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_pax_abd_pain <- dds_pax_abd_pain[filter,]
dds_pax_abd_pain <- DESeq(dds_pax_abd_pain)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_pax_abd_pain <- dds_pax_abd_pain[which(mcols(dds_pax_abd_pain)$betaConv),]
resultsNames(dds_pax_abd_pain) # lists the coefficients
dds_pax_abd_pain <- results(dds_pax_abd_pain, contrast=c("abd_pain","Y","N"))
dds_pax_abd_pain <- data.frame(EnsemblID=row.names(dds_pax_abd_pain), dds_pax_abd_pain)
dds_pax_abd_pain <- merge(dds_pax_abd_pain, genes_pax, by="row.names")
dds_pax_abd_pain <- dds_pax_abd_pain[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_pax_abd_pain <- subset(dds_pax_abd_pain, !is.na(padj))
dds_pax_abd_pain <- dds_pax_abd_pain[order(dds_pax_abd_pain$padj),]
dds_pax_abd_pain <- dds_pax_abd_pain %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_pax_abd_pain <- dds_pax_abd_pain[order(dds_pax_abd_pain$Expression, dds_pax_abd_pain$log2FoldChange),]
write.csv(dds_pax_abd_pain, "Statistical_Analyses/4_Symptoms/genes_pax_abd_pain.csv", row.names=FALSE)

# Anosmia
metadata_pax_anosmia <- subset(metadata_pax_Pos, !is.na(anosmia))
ids_pax_anosmia <- metadata_pax_anosmia$alias_sequencing_id
counts_pax_anosmia <- counts_pax %>% dplyr::select(all_of(ids_pax_anosmia))
nrow(metadata_pax_anosmia)
all(rownames(metadata_pax_anosmia) == colnames(counts_pax_anosmia))
dds_pax_anosmia <- DESeqDataSetFromMatrix(countData = counts_pax_anosmia, colData = metadata_pax_anosmia, design = ~ batch_num + sex + age_cat + anosmia + PC1 + PC2 + PC3 + PC4 + PC5)
# Require a normalized count of at least 10 in two or more samples
dds_pax_anosmia_norm <- estimateSizeFactors(dds_pax_anosmia)
dds_pax_anosmia_norm <- counts(dds_pax_anosmia_norm, normalized=TRUE)
filter <- rowSums(dds_pax_anosmia_norm >= 10) >= (nrow(metadata_pax_anosmia))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_pax_anosmia <- dds_pax_anosmia[filter,]
dds_pax_anosmia <- DESeq(dds_pax_anosmia)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_pax_anosmia <- dds_pax_anosmia[which(mcols(dds_pax_anosmia)$betaConv),]
resultsNames(dds_pax_anosmia) # lists the coefficients
dds_pax_anosmia <- results(dds_pax_anosmia, contrast=c("anosmia","Y","N"))
dds_pax_anosmia <- data.frame(EnsemblID=row.names(dds_pax_anosmia), dds_pax_anosmia)
dds_pax_anosmia <- merge(dds_pax_anosmia, genes_pax, by="row.names")
dds_pax_anosmia <- dds_pax_anosmia[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_pax_anosmia <- subset(dds_pax_anosmia, !is.na(padj))
dds_pax_anosmia <- dds_pax_anosmia[order(dds_pax_anosmia$padj),]
dds_pax_anosmia <- dds_pax_anosmia %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_pax_anosmia <- dds_pax_anosmia[order(dds_pax_anosmia$Expression, dds_pax_anosmia$log2FoldChange),]
write.csv(dds_pax_anosmia, "Statistical_Analyses/4_Symptoms/genes_pax_anosmia.csv", row.names=FALSE)

# Dysgeusia
metadata_pax_dysgeusia <- subset(metadata_pax_Pos, !is.na(dysgeusia))
ids_pax_dysgeusia <- metadata_pax_dysgeusia$alias_sequencing_id
counts_pax_dysgeusia <- counts_pax %>% dplyr::select(all_of(ids_pax_dysgeusia))
nrow(metadata_pax_dysgeusia)
all(rownames(metadata_pax_dysgeusia) == colnames(counts_pax_dysgeusia))
dds_pax_dysgeusia <- DESeqDataSetFromMatrix(countData = counts_pax_dysgeusia, colData = metadata_pax_dysgeusia, design = ~ batch_num + sex + age_cat + dysgeusia + PC1 + PC2 + PC3 + PC4 + PC5)
# Require a normalized count of at least 10 in two or more samples
dds_pax_dysgeusia_norm <- estimateSizeFactors(dds_pax_dysgeusia)
dds_pax_dysgeusia_norm <- counts(dds_pax_dysgeusia_norm, normalized=TRUE)
filter <- rowSums(dds_pax_dysgeusia_norm >= 10) >= (nrow(metadata_pax_dysgeusia))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_pax_dysgeusia <- dds_pax_dysgeusia[filter,]
dds_pax_dysgeusia <- DESeq(dds_pax_dysgeusia)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_pax_dysgeusia <- dds_pax_dysgeusia[which(mcols(dds_pax_dysgeusia)$betaConv),]
resultsNames(dds_pax_dysgeusia) # lists the coefficients
dds_pax_dysgeusia <- results(dds_pax_dysgeusia, contrast=c("dysgeusia","Y","N"))
dds_pax_dysgeusia <- data.frame(EnsemblID=row.names(dds_pax_dysgeusia), dds_pax_dysgeusia)
dds_pax_dysgeusia <- merge(dds_pax_dysgeusia, genes_pax, by="row.names")
dds_pax_dysgeusia <- dds_pax_dysgeusia[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_pax_dysgeusia <- subset(dds_pax_dysgeusia, !is.na(padj))
dds_pax_dysgeusia <- dds_pax_dysgeusia[order(dds_pax_dysgeusia$padj),]
dds_pax_dysgeusia <- dds_pax_dysgeusia %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_pax_dysgeusia <- dds_pax_dysgeusia[order(dds_pax_dysgeusia$Expression, dds_pax_dysgeusia$log2FoldChange),]
write.csv(dds_pax_dysgeusia, "Statistical_Analyses/4_Symptoms/genes_pax_dysgeusia.csv", row.names=FALSE)

# Myalgias
metadata_pax_myalgias <- subset(metadata_pax_Pos, !is.na(myalgias))
ids_pax_myalgias <- metadata_pax_myalgias$alias_sequencing_id
counts_pax_myalgias <- counts_pax %>% dplyr::select(all_of(ids_pax_myalgias))
nrow(metadata_pax_myalgias)
all(rownames(metadata_pax_myalgias) == colnames(counts_pax_myalgias))
dds_pax_myalgias <- DESeqDataSetFromMatrix(countData = counts_pax_myalgias, colData = metadata_pax_myalgias, design = ~ batch_num + sex + age_cat + myalgias + PC1 + PC2 + PC3 + PC4 + PC5)
# Require a normalized count of at least 10 in two or more samples
dds_pax_myalgias_norm <- estimateSizeFactors(dds_pax_myalgias)
dds_pax_myalgias_norm <- counts(dds_pax_myalgias_norm, normalized=TRUE)
filter <- rowSums(dds_pax_myalgias_norm >= 10) >= (nrow(metadata_pax_myalgias))*0.25 # require a normalized count of 10 in at least 25% of samples
dds_pax_myalgias <- dds_pax_myalgias[filter,]
dds_pax_myalgias <- DESeq(dds_pax_myalgias)
# Omit rows that did not converge (usually genes with very small counts & little power)
dds_pax_myalgias <- dds_pax_myalgias[which(mcols(dds_pax_myalgias)$betaConv),]
resultsNames(dds_pax_myalgias) # lists the coefficients
dds_pax_myalgias <- results(dds_pax_myalgias, contrast=c("myalgias","Y","N"))
dds_pax_myalgias <- data.frame(EnsemblID=row.names(dds_pax_myalgias), dds_pax_myalgias)
dds_pax_myalgias <- merge(dds_pax_myalgias, genes_pax, by="row.names")
dds_pax_myalgias <- dds_pax_myalgias[,c("EnsemblID","Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
dds_pax_myalgias <- subset(dds_pax_myalgias, !is.na(padj))
dds_pax_myalgias <- dds_pax_myalgias[order(dds_pax_myalgias$padj),]
dds_pax_myalgias <- dds_pax_myalgias %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated", log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged"))
dds_pax_myalgias <- dds_pax_myalgias[order(dds_pax_myalgias$Expression, dds_pax_myalgias$log2FoldChange),]
write.csv(dds_pax_myalgias, "Statistical_Analyses/4_Symptoms/genes_pax_myalgias.csv", row.names=FALSE)