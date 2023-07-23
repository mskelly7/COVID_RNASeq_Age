# BRAVE Kids RNA Sequencing Analysis
# Aditya Mohan (MD/PhD candidate) / Matthew Kelly, MD, MPH 
# Data preprocessing
# Last update: July 23, 2023

remove(list=ls())
setwd("_________________") 
set.seed(1234)

library(biomaRt)
library(readxl)
library(openxlsx)
library(tidyverse)
library(dplyr)
library(phyloseq)
library(DESeq2)
library(ggplot2)
library(ComplexHeatmap)
library(ggpubr)
library(data.table)
library(reshape2)
library(edgeR)
library(cowplot)
library(sva)
library(plyr)

###########################################################
# GENERATE INITIAL PHYLOSEQ OBJECT CONTAINING ALL SAMPLES #
###########################################################

# Create RNASeq count table
rnaseq <- read.csv("BRAVE_RNASeq_counts.csv", sep=",", stringsAsFactors = T, check.names=FALSE)
rnaseq$target_id <- gsub("\\..*", "", rnaseq$target_id)
nrow(rnaseq)
length(unique(rnaseq$target_id))
summary(colSums(rnaseq[,c(2:456)]))
rnaseq <- ddply(rnaseq, "target_id", numcolwise(sum))
rnaseq <- rnaseq %>% remove_rownames %>% column_to_rownames(var="target_id")
nrow(rnaseq)
summary(colSums(rnaseq))
ensembl_ids <- na.omit(as.list(rownames(rnaseq)))
ensembl_ids <- gsub("\\.*", "", ensembl_ids)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
gene_ids <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "ensembl_gene_id_version", "chromosome_name", "start_position", "end_position"), 
                  filters = "ensembl_gene_id", values = ensembl_ids, mart = mart)
gene_ids$Length <- gene_ids$end_position - gene_ids$start_position
summary(gene_ids$Length)
rnaseq$ensembl_gene_id <- row.names(rnaseq)
rnaseq_final <- merge(gene_ids, rnaseq, by="ensembl_gene_id")
rnaseq_final$ensembl_gene_id <- NULL
rnaseq_final$start_position <- NULL
rnaseq_final$end_position <- NULL
names(rnaseq_final)[names(rnaseq_final)=="ensembl_gene_id_version"] <- "EnsemblID"
names(rnaseq_final)[names(rnaseq_final)=="external_gene_name"] <- "Gene"
names(rnaseq_final)[names(rnaseq_final)=="chromosome_name"] <- "Chromosome"
remove(mart, rnaseq, ensembl_ids, gene_ids)

# Create gene count table (otu_table)
raw_counts <- rnaseq_final[,-c(1,3,4)] # dropping Gene, Chromosome, Length
raw_counts <- raw_counts %>% remove_rownames %>% column_to_rownames(var="EnsemblID")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X025CBE.M0.PAX", replacement = "025CBE.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X0719FA.M0.PAX", replacement = "0719FA.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X0BF51C.M0.PAX", replacement = "0BF51C.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X0BF51C.M6.PAX", replacement = "0BF51C.M6.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X0E1ECA.M0.PAX", replacement = "0E1ECA.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X12CE1E.M0.PAX", replacement = "12CE1E.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X150B4B.M0.PAX", replacement = "150B4B.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X17A783.M0.PAX", replacement = "17A783.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X18EF5E.M0.PAX", replacement = "18EF5E.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X1AC6B8.M0.PAX", replacement = "1AC6B8.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X1E06CA.M0.PAX", replacement = "1E06CA.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X22D142.M0.PAX", replacement = "22D142.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X234BBF.M0.PAX", replacement = "234BBF.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X322088.M0.PAX", replacement = "322088.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X3501A9.M0.PAX", replacement = "3501A9.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X3B441F.M0.PAX", replacement = "3B441F.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X496933.M0.PAX", replacement = "496933.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X5217A5.M0.PAX", replacement = "5217A5.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X58A568.M0.PAX", replacement = "58A568.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X6AA2F7.M0.PAX", replacement = "6AA2F7.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X71549F.M0.PAX", replacement = "71549F.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X72C971.M0.PAX", replacement = "72C971.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X7400BF.M0.PAX", replacement = "7400BF.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X7D880A.M0.PAX", replacement = "7D880A.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X826E5F.M0.PAX", replacement = "826E5F.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X830937.M0.PAX", replacement = "830937.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X85F0AB.M0.PAX", replacement = "85F0AB.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X8655BF.M0.PAX", replacement = "8655BF.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X8867DF.M0.PAX", replacement = "8867DF.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X88732B.M0.PAX", replacement = "88732B.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X88B192.M0.PAX", replacement = "88B192.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X8CF0CC.M0.PAX", replacement = "8CF0CC.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X8F0016.M0.PAX", replacement = "8F0016.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X8FDB69.M0.PAX", replacement = "8FDB69.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X990CA3.M0.PAX", replacement = "990CA3.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X9C0AA4.M0.PAX", replacement = "9C0AA4.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X9E259D.M0.PAX", replacement = "9E259D.M0.PAX")
names(raw_counts) <- gsub(x = names(raw_counts), pattern = "X9FAB68.M0.PAX", replacement = "9FAB68.M0.PAX")
raw_counts_otu <- otu_table(raw_counts, taxa_are_rows = TRUE)

# Create gene annotation table (tax_table)
genes <- rnaseq_final[,c(1:4)] # Keeping Gene, EnsemblID, Chromosome, Length
genes <- genes %>% remove_rownames %>% column_to_rownames(var="EnsemblID")
genes_tax <- tax_table(as.matrix(genes))

# Upload metadata file
metadata_brave <- read_excel("BRAVE_RNASeq_metadata.xlsx", sheet="brave", na="NA")
metadata_messi <- read_excel("BRAVE_RNASeq_metadata.xlsx", sheet="messi", na="NA")
rnaseq_metadata <- merge(metadata_brave, metadata_messi, by=c("alias_sequencing_id", "alias_study_id", "SampleTiming", "SampleType", "year", "study", "hospital", 
                                                              "batch_num", "RIN", "age", "age_cat", "sex", "race", "hispanic", "corona", "group", "clinical_pcr", 
                                                              "research_pcr", "vl_copies", "vaccine_doses", "timing_sx", "timing_dx",
                                                              "comorbidity", "obesity", "asthma", "pulm_oth", "htn", "cardiac_oth", "diabetes", 
                                                              "neuro", "renal", "cancer", "immuno", "symptoms_any", "fever", "cough", "rhinorrhea", "congestion"), 
                                                              all.x=TRUE, all.y=TRUE)
rnaseq_metadata <- rnaseq_metadata[,c("alias_sequencing_id", "alias_study_id", "SampleTiming", "SampleType", "year", "study", "hospital", 
                                      "batch_num", "RIN", "age", "age_cat", "sex", "race", "hispanic", "corona", "group", "clinical_pcr", 
                                      "research_pcr", "vl_copies", "vaccine_doses", "timing_sx", "timing_dx",
                                      "comorbidity", "obesity", "asthma", "pulm_oth", "htn", "cardiac_oth", "diabetes", 
                                      "neuro", "renal", "cancer", "immuno", "symptoms_any", "fever", "cough", "rhinorrhea", "congestion")]

rnaseq_metadata$vl_copies <- as.numeric(rnaseq_metadata$vl_copies)
row.names(rnaseq_metadata) <- rnaseq_metadata$alias_sequencing_id

# Assemble initial phyloseq object
phy.rnaseq.all <- phyloseq(raw_counts_otu, genes_tax)
sample_data(phy.rnaseq.all) <- rnaseq_metadata
nsamples(phy.rnaseq.all)
ntaxa(phy.rnaseq.all)
summary(sample_sums(phy.rnaseq.all))
# All samples have >5 million reads which is a good bare minimum for differential expression analysis using human samples
saveRDS(phy.rnaseq.all, "phy.rnaseq.all.rds")

########################################################################
# EVALUATE FOR MISCLASSICATION OF SAMPLES USING X & Y CHROMOSOME GENES #
########################################################################

#metadata_all <- data.frame(sample_data(phy.rnaseq.all))
#metadata_all$batch_num <- as.factor(metadata_all$batch_num)
#metadata_all$corona <- as.factor(metadata_all$corona)

# Check the gender of samples using X & Y chromosome genes
#dds.rnaseq <- DESeqDataSetFromMatrix(countData = raw_counts, colData = metadata_all, design =  ~ batch_num + corona)  
#dds.rnaseq <- DESeq(dds.rnaseq)
#vsd <- vst(dds.rnaseq) # variance stabilizing transformation
#normCount <- counts(dds.rnaseq, normalized=TRUE) # normalized count
#genes <- data.frame(tax_table(phy.rnaseq.all))
#genes$EnsemblID <- row.names(genes)
#y_gene <- genes %>% filter(Chromosome == "Y" & EnsemblID %in% row.names(genes)) %>% pull(EnsemblID)
#x_gene <- genes %>% filter(Gene == "XIST") %>% pull(EnsemblID)
#y_counts <- colSums(normCount[y_gene,])
#x_counts <- normCount[x_gene,]
#gender_df <- cbind(metadata_all, x_counts, y_counts)
#gender_df <- gender_df[order(gender_df$sex),]

# Inspect plots to confirm that all samples have RNA-Seq data that matches patient gender

#y_plot <- ggplot(gender_df, aes(x=alias_study_id, y=y_counts, fill=sex)) + geom_bar(stat="identity") +
#  labs(title="Genes on the Y chromosome", x="Sample", y="Total Counts") + 
#  theme(plot.title = element_text(hjust = 0.5, size = 10), axis.text.x=element_text(angle=45, hjust = 1, size=6))
#png(file="R_Plots/Y_Genes.png", 
#    width = 30, height = 10, units = 'in', res = 600)
#print(y_plot)
#dev.off()

#x_plot <- ggplot(gender_df, aes(x=alias_study_id, y=((x_counts)/(y_counts+1)), fill=sex)) + geom_bar(stat="identity") + 
#  labs(title="XIST gene (female marker)", x="Sample", y="Counts Relative to Y X-Some Counts") + 
#  theme(plot.title = element_text(hjust = 0.5, size = 10), axis.text.x=element_text(angle=45, hjust = 1, size=6))
#png(file="R_Plots/X_Genes.png", 
#    width = 30, height = 10, units = 'in', res = 600)
#print(x_plot)
#dev.off()

##################
# SAMPLE PRUNING #
##################

phy.rnaseq.all <- readRDS("phy.rnaseq.all.rds")
nsamples(phy.rnaseq.all)
ntaxa(phy.rnaseq.all)
# Exclude samples from subjects with SARS-CoV-2 infections that required hospitalization (n=3)
phy.rnaseq.pruned <- subset_samples(phy.rnaseq.all, hospital!="Y")
nsamples(phy.rnaseq.pruned)
# Exclude nasal swabs (n=14)
phy.rnaseq.pruned <- subset_samples(phy.rnaseq.pruned, SampleType!="nasal")
nsamples(phy.rnaseq.pruned)
# Exclude 1-month convalescent samples (n=14)
phy.rnaseq.pruned <- subset_samples(phy.rnaseq.pruned, SampleTiming!="convalescent_1mo")
nsamples(phy.rnaseq.pruned)
# Remove samples from symptomatic SARS-CoV-2-uninfected individuals (n=15) 
phy.rnaseq.pruned <- subset_samples(phy.rnaseq.pruned, group!="NEG_SX")
nsamples(phy.rnaseq.pruned)
# Remove SARS-CoV-2-uninfected subjects in BRAVE study with single-replicate positive PCR testing in Denny Lab (n=6) 
phy.rnaseq.pruned <- subset_samples(phy.rnaseq.pruned, corona=="Positive" | vl_copies==0 | is.na(vl_copies))
nsamples(phy.rnaseq.pruned)
# Remove individuals with any history of SARS-CoV-2 vaccination, all of whom are SARS-CoV-2+ (n=21)
phy.rnaseq.pruned <- subset_samples(phy.rnaseq.pruned, vaccine_doses==0)
nsamples(phy.rnaseq.pruned)
# Focus only on COVID+ individuals enrolled in 2020 with likely wild-type virus infection (n=47)
phy.rnaseq.pruned <- subset_samples(phy.rnaseq.pruned, corona=="Negative" | year=="2020")
nsamples(phy.rnaseq.pruned)
# Remove samples from COVID+ collected prior to or >14 days from SARS-CoV-2 diagnosis (n=0)
phy.rnaseq.pruned <- subset_samples(phy.rnaseq.pruned, is.na(timing_dx) | (timing_dx>=0 & timing_dx<=14))
nsamples(phy.rnaseq.pruned)
# Remove samples from COVID+ collected >14 days from symptom onset (n=3)
phy.rnaseq.pruned <- subset_samples(phy.rnaseq.pruned, is.na(timing_sx) | (timing_sx>=-3 & timing_sx<=14))
nsamples(phy.rnaseq.pruned)
# Remove LE3362 who had abdominal pain but prolonged and preceded COVID diagnosis (n=2)
phy.rnaseq.pruned <- subset_samples(phy.rnaseq.pruned, alias_study_id!="LE3362")
nsamples(phy.rnaseq.pruned)
# Remove single NP sample in this dataset from an adult (n=1)
phy.rnaseq.pruned <- subset_samples(phy.rnaseq.pruned, SampleType!="np" | age_cat!="Adult")
nsamples(phy.rnaseq.pruned)
# Remove genes that are present in no samples
ntaxa(phy.rnaseq.pruned)
phy.rnaseq.pruned <- prune_taxa(taxa_sums(phy.rnaseq.pruned) > 0, phy.rnaseq.pruned)
ntaxa(phy.rnaseq.pruned)

# Generate tables summarizing final study population by age category and COVID status
metadata_pruned <- data.frame(sample_data(phy.rnaseq.pruned))
metadata_pruned_np <- subset(metadata_pruned, SampleType=="np")
table(metadata_pruned_np$age_cat, metadata_pruned_np$corona)
metadata_pruned_pax <- subset(metadata_pruned, SampleType=="pax")
table(metadata_pruned_pax$age_cat, metadata_pruned_pax$corona)
saveRDS(phy.rnaseq.pruned, "phy.rnaseq.pruned.rds")

############################################
# ADD CELL POPULATION DATA FROM CIBERSORTX #
############################################

# CIBERSORT is a deconvolution-based method that uses the v-support vector regression method (v-SVR) to estimate each immune cell type’s relative proportion from 
# a gene expression profile (GEP). 

# The “mixture file” is a single tab-delimited text file containing 1 or more gene expression profiles (GEPs) of biological mixture samples. The first column 
# contains gene names and should have “Gene” (or similar) as a column header (i.e., in the space occupying column 1, row 1). Multiple samples may be analyzed 
# in parallel, with the remaining columns (2, 3, etc.) dedicated to mixture GEPs, where each row represents the expression value for a given gene and the column 
# header is the name of the mixture sample. Note that the mixture file and the signature matrix must share the same naming scheme for gene identifiers.

# The “signature matrix” is a tab-delimited text file consisting of sets of “barcode genes” whose expression values collectively define unique gene expression 
# signatures for each cell subset of interest. The file format is similar to the mixture file, with gene names in column 1, however the remaining columns consist 
# of signature GEPs from individual cell subsets. A validated leukocyte gene signature matrix (LM22) is available for the deconvolution of 22 functionally defined 
# human hematopoietic subsets. LM22 was generated using Affymetrix HGU133A microarray data and has been rigorously tested on Affymetrix HGU133 and Illumina 
# Beadchip platforms.

phy.rnaseq.pruned <- readRDS("phy.rnaseq.pruned.rds")
rnaseq_metadata <- data.frame(sample_data(phy.rnaseq.pruned))
raw_counts_pruned <- data.frame(otu_table(phy.rnaseq.pruned))
genes_pruned <- data.frame(tax_table(phy.rnaseq.pruned))
genes_pruned$Length <- as.numeric(genes_pruned$Length)

# Upload LM22 signature matrix for inspection
lm22 <- read.delim("Cibersortx/LM22.txt", header = TRUE, sep = "\t", dec = ".")
lm22_genes <- lm22$Gene.symbol

# Convert read counts to transcripts per kilobase million (TPMs)
# Function to calculate TPMs from a table of raw gene counts (https://github.com/davidrequena/drfun/blob/main/R/tpm.R)
# The input table is numeric:
# - The row names are the gene identifiers (ensemblID).
# - The column names represent the samples.
# The gene lengths are in a column of a dataframe with the same row order.
tpm <- function(raw_counts, gene_lengths) {
  x <- raw_counts*1e3 / gene_lengths
  return(t(t(x)*1e6 / colSums(x)))
}

tpm_counts <- data.frame(tpm(raw_counts_pruned, genes_pruned$Length))
tpm_counts <- merge(genes_pruned, tpm_counts, by='row.names', all=TRUE)

# Confirm that all genes in signature matrix are only present once in dataframe
tpm_test <- subset(tpm_counts, Gene %in% lm22_genes)
nrow(tpm_test)
length(unique(tpm_test$Gene))
row.names(tpm_counts) <- tpm_counts$Row.names # Set EnsemblIDs as rownames again
tpm_counts <- tpm_counts[,-c(1,3,4)] # Drop Row.names, Chromosome, Length
write.table(tpm_counts, file = "Cibersortx/mixture_file.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

# Run this count matrix using the "Impute Cell Fractions" analysis module on the CIBERSORTx website
# Run with B-mode batch correction, disable quantile normalization, run in relative mode, and set permutations to 500 (recommended is >=100)
cibersort <- read.delim("Cibersortx/Relative_mode/CIBERSORTx_Results.txt", header = TRUE, sep = "\t", dec = ".")
names(cibersort) <- gsub(x = names(cibersort), pattern = "\\.", replacement = "_") 
cibersort$Mixture[cibersort$Mixture=="X025CBE.M0.PAX"] <- "025CBE.M0.PAX"
cibersort$Mixture[cibersort$Mixture=="X0719FA.M0.PAX"] <- "0719FA.M0.PAX"
cibersort$Mixture[cibersort$Mixture=="X0BF51C.M0.PAX"] <- "0BF51C.M0.PAX"
cibersort$Mixture[cibersort$Mixture=="X0BF51C.M6.PAX"] <- "0BF51C.M6.PAX"
cibersort$Mixture[cibersort$Mixture=="X0E1ECA.M0.PAX"] <- "0E1ECA.M0.PAX"
cibersort$Mixture[cibersort$Mixture=="X18EF5E.M0.PAX"] <- "18EF5E.M0.PAX"
cibersort$Mixture[cibersort$Mixture=="X1AC6B8.M0.PAX"] <- "1AC6B8.M0.PAX"
cibersort$Mixture[cibersort$Mixture=="X1E06CA.M0.PAX"] <- "1E06CA.M0.PAX"
cibersort$Mixture[cibersort$Mixture=="X234BBF.M0.PAX"] <- "234BBF.M0.PAX"
cibersort$Mixture[cibersort$Mixture=="X322088.M0.PAX"] <- "322088.M0.PAX"
cibersort$Mixture[cibersort$Mixture=="X3B441F.M0.PAX"] <- "3B441F.M0.PAX"
cibersort$Mixture[cibersort$Mixture=="X496933.M0.PAX"] <- "496933.M0.PAX"
cibersort$Mixture[cibersort$Mixture=="X58A568.M0.PAX"] <- "58A568.M0.PAX"
cibersort$Mixture[cibersort$Mixture=="X6AA2F7.M0.PAX"] <- "6AA2F7.M0.PAX"
cibersort$Mixture[cibersort$Mixture=="X71549F.M0.PAX"] <- "71549F.M0.PAX"
cibersort$Mixture[cibersort$Mixture=="X72C971.M0.PAX"] <- "72C971.M0.PAX"
cibersort$Mixture[cibersort$Mixture=="X7D880A.M0.PAX"] <- "7D880A.M0.PAX"
cibersort$Mixture[cibersort$Mixture=="X826E5F.M0.PAX"] <- "826E5F.M0.PAX"
cibersort$Mixture[cibersort$Mixture=="X830937.M0.PAX"] <- "830937.M0.PAX"
cibersort$Mixture[cibersort$Mixture=="X8867DF.M0.PAX"] <- "8867DF.M0.PAX"
cibersort$Mixture[cibersort$Mixture=="X88732B.M0.PAX"] <- "88732B.M0.PAX"
cibersort$Mixture[cibersort$Mixture=="X88B192.M0.PAX"] <- "88B192.M0.PAX"
row.names(cibersort) <- cibersort$Mixture
cibersort <- cibersort[,c(2:23)]
cibersort$B_cells <- cibersort$B_cells_naive + cibersort$B_cells_memory
cibersort$T_cells_CD4 <- cibersort$T_cells_CD4_naive + cibersort$T_cells_CD4_memory_resting + cibersort$T_cells_CD4_memory_activated + cibersort$T_cells_follicular_helper +
                       cibersort$T_cells_follicular_helper
cibersort$NK_cells <- cibersort$NK_cells_resting + cibersort$NK_cells_activated
cibersort$Mono_Macrophages <- cibersort$Monocytes + cibersort$Macrophages_M0 + cibersort$Macrophages_M1 + cibersort$Macrophages_M2
cibersort$Dendritic_cells <- cibersort$Dendritic_cells_resting + cibersort$Dendritic_cells_activated
cibersort$Mast_cells <- cibersort$Mast_cells_resting + cibersort$Mast_cells_activated
rnaseq_metadata <- merge(rnaseq_metadata, cibersort, by='row.names')
row.names(rnaseq_metadata) <- rnaseq_metadata$Row.names
rnaseq_metadata$Row.names <- NULL
sample_data(phy.rnaseq.pruned) <- data.frame(rnaseq_metadata)
phy.rnaseq.cibersort <- phy.rnaseq.pruned
saveRDS(phy.rnaseq.cibersort, "phy.rnaseq.cibersort.rds")
write.csv(rnaseq_metadata, "rnaseq_metadata.csv", row.names=FALSE)

##################
# GENE FILTERING #
##################

phy.rnaseq.cibersort <- readRDS("phy.rnaseq.cibersort.rds")
rnaseq_genes <- data.frame(tax_table(phy.rnaseq.cibersort))
nsamples(phy.rnaseq.cibersort)
ntaxa(phy.rnaseq.cibersort)

# Remove transcripts without gene identifier
phy.rnaseq.cibersort <- subset_taxa(phy.rnaseq.cibersort, Gene!="")
ntaxa(phy.rnaseq.cibersort)

# Remove mitochondrial genes
MT_genes <- rnaseq_genes[grep("^MT", rnaseq_genes$Gene, perl=T),]
MT_genes <- MT_genes$Gene
phy.rnaseq.cibersort <- subset_taxa(phy.rnaseq.cibersort, !(Gene %in% MT_genes))
ntaxa(phy.rnaseq.cibersort)
# Remove L ribosomal genes
RPL_genes <- rnaseq_genes[grep("^RPL", rnaseq_genes$Gene, perl=T),]
RPL_genes <- RPL_genes$Gene
phy.rnaseq.cibersort <- subset_taxa(phy.rnaseq.cibersort, !(Gene %in% RPL_genes))
ntaxa(phy.rnaseq.cibersort)
# Remove S ribosomal genes
RPS_genes <- rnaseq_genes[grep("^RPS", rnaseq_genes$Gene, perl=T),]
RPS_genes <- RPS_genes$Gene
phy.rnaseq.cibersort <- subset_taxa(phy.rnaseq.cibersort, !(Gene %in% RPS_genes))
ntaxa(phy.rnaseq.cibersort)
# Remove large subunit mitochondrial ribosomal proteins
MRPL_genes <- rnaseq_genes[grep("^MRPL", rnaseq_genes$Gene, perl=T),]
MRPL_genes <- MRPL_genes$Gene
phy.rnaseq.cibersort <- subset_taxa(phy.rnaseq.cibersort, !(Gene %in% MRPL_genes))
ntaxa(phy.rnaseq.cibersort)
# Remove large subunit mitochondrial ribosomal proteins
MRPS_genes <- rnaseq_genes[grep("^MRPS", rnaseq_genes$Gene, perl=T),]
MRPS_genes <- MRPS_genes$Gene
phy.rnaseq.cibersort <- subset_taxa(phy.rnaseq.cibersort, !(Gene %in% MRPS_genes))
ntaxa(phy.rnaseq.cibersort)
# Remove additional rRNA genes
rRNA_genes <- rnaseq_genes[grep("rRNA", rnaseq_genes$Gene, perl=T),]
rRNA_genes <- rRNA_genes$Gene
phy.rnaseq.cibersort <- subset_taxa(phy.rnaseq.cibersort, !(Gene %in% rRNA_genes))
ntaxa(phy.rnaseq.cibersort)
RNA5_genes <- rnaseq_genes[grep("RNA5", rnaseq_genes$Gene, perl=T),]
RNA5_genes <- RNA5_genes$Gene
phy.rnaseq.cibersort <- subset_taxa(phy.rnaseq.cibersort, !(Gene %in% RNA5_genes))
# Remove non-coding RNA genes
Y_RNA_genes <- rnaseq_genes[grep("^Y_RNA", rnaseq_genes$Gene, perl=T),]
Y_RNA_genes <- Y_RNA_genes$Gene
phy.rnaseq.cibersort <- subset_taxa(phy.rnaseq.cibersort, !(Gene %in% Y_RNA_genes))
ntaxa(phy.rnaseq.cibersort)
# Remove small nuclear RNA genes
RNU_genes <- rnaseq_genes[grep("RNU", rnaseq_genes$Gene, perl=T),]
RNU_genes <- RNU_genes$Gene
phy.rnaseq.cibersort <- subset_taxa(phy.rnaseq.cibersort, !(Gene %in% RNU_genes))
ntaxa(phy.rnaseq.cibersort)
rnaseq_genes_filtered <- data.frame(tax_table(phy.rnaseq.cibersort))
write.csv(rnaseq_genes_filtered, "rnaseq_genes_filtered.csv")

### NP SAMPLES ###

# Remove genes not expressed at level of at least 1 per million reads in at least 50% of samples
phy.rnaseq.np <- subset_samples(phy.rnaseq.cibersort, SampleType=="np")
nsamples(phy.rnaseq.np)
summary(sample_sums(phy.rnaseq.np))
counts_np <- data.frame(otu_table(phy.rnaseq.np))
summary(rowSums(counts_np))
counts_np_cpm <- cpm(counts_np)
phy.rnaseq.np.cpm <- phy.rnaseq.np
counts_np_cpm <- otu_table(counts_np_cpm, taxa_are_rows = TRUE)
otu_table(phy.rnaseq.np.cpm) <- counts_np_cpm 
phy.rnaseq.np.cpm <- filter_taxa(phy.rnaseq.np.cpm, function(x) (sum(x>0.000001)) > (0.5*length(x)), TRUE)
genes_to_keep <- taxa_names(phy.rnaseq.np.cpm)
phy.rnaseq.np <- subset_taxa(phy.rnaseq.np, taxa_names(phy.rnaseq.np) %in% genes_to_keep)
ntaxa(phy.rnaseq.np)
summary(sample_sums(phy.rnaseq.np))

# Generate PCA plot to evaluate for batch effect (will adjust for this in DESeq2)
metadata_np <- data.frame(sample_data(phy.rnaseq.np))
counts_raw_np <- data.frame(otu_table(phy.rnaseq.np))
DGE_raw_np <- DGEList(counts = counts_raw_np, samples = metadata_np)
DGE_raw_np_cpm <- cpm(DGE_raw_np, log=TRUE)
# Create PCA plot by batch
DGE_raw_np_pca <- prcomp(t(DGE_raw_np_cpm), center=TRUE, scale=TRUE)
batch_num <- as.factor(DGE_raw_np$samples$batch_num)
counts_raw_np_pca <- data.frame(pc1=DGE_raw_np_pca$x[, 1], pc2=DGE_raw_np_pca$x[, 2], group=batch_num)
batch_np_plot <- ggplot(counts_raw_np_pca, aes(x=pc1, y=pc2)) + geom_point(aes(colour=group), size=4) + labs(title="NP samples by batch") + 
  theme_bw(18) + theme(plot.title = element_text(size = 15, face="bold", hjust=0.5))
batch_np_plot

### PAX SAMPLES ###

# Remove genes not expressed at level of at least 1 per million reads in at least 50% of samples
phy.rnaseq.pax <- subset_samples(phy.rnaseq.cibersort, SampleType=="pax")
nsamples(phy.rnaseq.pax)
summary(sample_sums(phy.rnaseq.pax))
counts_pax <- data.frame(otu_table(phy.rnaseq.pax))
summary(rowSums(counts_pax))
counts_pax_cpm <- cpm(counts_pax)
phy.rnaseq.pax.cpm <- phy.rnaseq.pax
counts_pax_cpm <- otu_table(counts_pax_cpm, taxa_are_rows = TRUE)
otu_table(phy.rnaseq.pax.cpm) <- counts_pax_cpm 
phy.rnaseq.pax.cpm <- filter_taxa(phy.rnaseq.pax.cpm, function(x) (sum(x>0.000001)) > (0.5*length(x)), TRUE)
genes_to_keep <- taxa_names(phy.rnaseq.pax.cpm)
phy.rnaseq.pax <- subset_taxa(phy.rnaseq.pax, taxa_names(phy.rnaseq.pax) %in% genes_to_keep)
ntaxa(phy.rnaseq.pax)
summary(sample_sums(phy.rnaseq.pax))

# Generate PCA plot to evaluate for batch effect (will adjust for this in DESeq2)
metadata_pax <- data.frame(sample_data(phy.rnaseq.pax))
counts_raw_pax <- data.frame(otu_table(phy.rnaseq.pax))
DGE_raw_pax <- DGEList(counts = counts_raw_pax, samples = metadata_pax)
DGE_raw_pax_cpm <- cpm(DGE_raw_pax, log=TRUE)
# Create PCA plot by batch
DGE_raw_pax_pca <- prcomp(t(DGE_raw_pax_cpm), center=TRUE, scale=TRUE)
batch_num <- as.factor(DGE_raw_pax$samples$batch_num)
counts_raw_pax_pca <- data.frame(pc1=DGE_raw_pax_pca$x[, 1], pc2=DGE_raw_pax_pca$x[, 2], group=batch_num)
batch_pax_plot <- ggplot(counts_raw_pax_pca, aes(x=pc1, y=pc2)) + geom_point(aes(colour=group), size=4) + labs(title="PAX samples by batch") + 
  theme_bw(18) + theme(plot.title = element_text(size = 15, face="bold", hjust=0.5))
batch_pax_plot

########################################################
### CREATE PRINCIPAL COMPONENTS FOR CELL POPULATIONS ###
########################################################

metadata_np_pca <- metadata_np[,c("B_cells", "Plasma_cells", "T_cells_CD4", "T_cells_CD8", "T_cells_gamma_delta", "NK_cells",
                                  "Mono_Macrophages", "Dendritic_cells", "Mast_cells", "Eosinophils", "Neutrophils")]
np_pca <- prcomp(metadata_np_pca, center = TRUE, scale = TRUE)
print(np_pca)
summary(np_pca)
# PC1:PC5 explain 76.6% of variability in cell populations in upper respiratory samples
str(np_pca)
np_PC1toPC5 <- np_pca$x[,1:5]
cibersort_np_pca <- merge(metadata_np, np_PC1toPC5, by="row.names")
cibersort_np_pca$Row.names <- NULL
row.names(cibersort_np_pca) <- cibersort_np_pca$alias_sequencing_id
sample_data(phy.rnaseq.np) <- cibersort_np_pca
saveRDS(phy.rnaseq.np, "phy.rnaseq.np.rds")

metadata_pax_pca <- metadata_pax[,c("B_cells", "Plasma_cells", "T_cells_CD4", "T_cells_CD8", "T_cells_gamma_delta", "NK_cells",
                                    "Mono_Macrophages", "Dendritic_cells", "Mast_cells", "Eosinophils", "Neutrophils")]
summary(metadata_pax_pca$T_cells_gamma_delta)
metadata_pax_pca$T_cells_gamma_delta <- NULL
summary(metadata_pax_pca$Eosinophils)
metadata_pax_pca$Eosinophils <- NULL
pax_pca <- prcomp(metadata_pax_pca, center = TRUE, scale = TRUE)
print(pax_pca)
summary(pax_pca)
# PC1:PC5 explain 83.3% of variability in cell populations in peripheral blood samples
str(pax_pca)
pax_PC1toPC5 <- pax_pca$x[,1:5]
cibersort_pax_pca <- merge(metadata_pax, pax_PC1toPC5, by="row.names")
cibersort_pax_pca$Row.names <- NULL
row.names(cibersort_pax_pca) <- cibersort_pax_pca$alias_sequencing_id
sample_data(phy.rnaseq.pax) <- cibersort_pax_pca
saveRDS(phy.rnaseq.pax, "phy.rnaseq.pax.rds")