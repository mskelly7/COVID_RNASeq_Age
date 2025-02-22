# BRAVE Kids RNA Sequencing Analysis
# Aditya Mohan (MD/PhD candidate)  / Matthew Kelly, MD, MPH 
# Gene Set Enrichment Analysis using FGSEA
# Last update: February 21, 2025

remove(list=ls())
setwd("___________________") 
set.seed(1234)
getRversion()
options(rstudio.help.showDataPreview = FALSE)

if(any(grepl("package:plyr", search()))) detach("package:plyr") else message("plyr not loaded")
library(readr)
library(dplyr)
library(readxl)
library(writexl)
library(tidyr) 
library(ggpubr)
library(data.table)
library(tibble)
library(fgsea)
packageVersion("fgsea")

# Upload module information 
modules_61 <- gmtPathways("Statistical_Analyses/modules_61.gmt")

# Upload DESeq2 output files

dds_np_neg_0to5_vs_6to13_nocibersort <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_np_neg_0to5_vs_6to13_nocibersort.csv")
dds_np_neg_0to5_vs_14to20_nocibersort <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_np_neg_0to5_vs_14to20_nocibersort.csv")
dds_np_neg_6to13_vs_14to20_nocibersort <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_np_neg_6to13_vs_14to20_nocibersort.csv")
dds_pax_neg_0to5_vs_6to13_cibersort <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_pax_neg_0to5_vs_6to13_cibersort.csv")
dds_pax_neg_0to5_vs_14to20_cibersort <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_pax_neg_0to5_vs_14to20_cibersort.csv")
dds_pax_neg_6to13_vs_14to20_cibersort <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_pax_neg_6to13_vs_14to20_cibersort.csv")
dds_pax_neg_0to5_vs_adult_cibersort <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_pax_neg_0to5_vs_adult_cibersort.csv")
dds_pax_neg_6to13_vs_adult_cibersort <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_pax_neg_6to13_vs_adult_cibersort.csv")
dds_pax_neg_14to20_vs_adult_cibersort <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_pax_neg_14to20_vs_adult_cibersort.csv")

dds_np_pos_neg_nocibersort <- read.csv("Statistical_Analyses/2_COVID_Pos_vs_Neg/genes_np_pos_neg_nocibersort.csv")
dds_pax_pos_neg_cibersort <- read.csv("Statistical_Analyses/2_COVID_Pos_vs_Neg/genes_pax_pos_neg_cibersort.csv")

dds_np_0to5_pos_neg_nocibersort <- read.csv("Statistical_Analyses/3_COVID_Pos_by_Age/genes_np_0to5_pos_neg_nocibersort.csv")
dds_np_6to13_pos_neg_nocibersort <- read.csv("Statistical_Analyses/3_COVID_Pos_by_Age/genes_np_6to13_pos_neg_nocibersort.csv")
dds_np_14to20_pos_neg_nocibersort <- read.csv("Statistical_Analyses/3_COVID_Pos_by_Age/genes_np_14to20_pos_neg_nocibersort.csv")
dds_pax_0to5_pos_neg_cibersort <- read.csv("Statistical_Analyses/3_COVID_Pos_by_Age/genes_pax_0to5_pos_neg_cibersort.csv")
dds_pax_6to13_pos_neg_cibersort <- read.csv("Statistical_Analyses/3_COVID_Pos_by_Age/genes_pax_6to13_pos_neg_cibersort.csv")
dds_pax_14to20_pos_neg_cibersort <- read.csv("Statistical_Analyses/3_COVID_Pos_by_Age/genes_pax_14to20_pos_neg_cibersort.csv")
dds_pax_adult_pos_neg_cibersort <- read.csv("Statistical_Analyses/3_COVID_Pos_by_Age/genes_pax_adult_pos_neg_cibersort.csv")

dds_np_fever_nocibersort <- read.csv("Statistical_Analyses/4_Symptoms/genes_np_fever_nocibersort.csv")
dds_np_cough_nocibersort <- read.csv("Statistical_Analyses/4_Symptoms/genes_np_cough_nocibersort.csv")
dds_np_rhinorrhea_nocibersort <- read.csv("Statistical_Analyses/4_Symptoms/genes_np_rhinorrhea_nocibersort.csv")
dds_np_congestion_nocibersort <- read.csv("Statistical_Analyses/4_Symptoms/genes_np_congestion_nocibersort.csv")
dds_np_headache_nocibersort <- read.csv("Statistical_Analyses/4_Symptoms/genes_np_headache_nocibersort.csv")
dds_np_abd_pain_nocibersort <- read.csv("Statistical_Analyses/4_Symptoms/genes_np_abd_pain_nocibersort.csv")
dds_np_anosmia_nocibersort <- read.csv("Statistical_Analyses/4_Symptoms/genes_np_anosmia_nocibersort.csv")
dds_np_dysgeusia_nocibersort <- read.csv("Statistical_Analyses/4_Symptoms/genes_np_dysgeusia_nocibersort.csv")
dds_np_myalgias_nocibersort <- read.csv("Statistical_Analyses/4_Symptoms/genes_np_myalgias_nocibersort.csv")
dds_pax_fever_cibersort <- read.csv("Statistical_Analyses/4_Symptoms/genes_pax_fever_cibersort.csv")
dds_pax_cough_cibersort <- read.csv("Statistical_Analyses/4_Symptoms/genes_pax_cough_cibersort.csv")
dds_pax_rhinorrhea_cibersort <- read.csv("Statistical_Analyses/4_Symptoms/genes_pax_rhinorrhea_cibersort.csv")
dds_pax_congestion_cibersort <- read.csv("Statistical_Analyses/4_Symptoms/genes_pax_congestion_cibersort.csv")
dds_pax_headache_cibersort <- read.csv("Statistical_Analyses/4_Symptoms/genes_pax_headache_cibersort.csv")
dds_pax_abd_pain_cibersort <- read.csv("Statistical_Analyses/4_Symptoms/genes_pax_abd_pain_cibersort.csv")
dds_pax_anosmia_cibersort <- read.csv("Statistical_Analyses/4_Symptoms/genes_pax_anosmia_cibersort.csv")
dds_pax_dysgeusia_cibersort <- read.csv("Statistical_Analyses/4_Symptoms/genes_pax_dysgeusia_cibersort.csv")
dds_pax_myalgias_cibersort <- read.csv("Statistical_Analyses/4_Symptoms/genes_pax_myalgias_cibersort.csv")

#####################################################################
###  DIFFERENCES IN NP MODULE EXPRESSION BY AGE IN COVID-NEGATIVE ###
#####################################################################

res_np_neg_0to5_vs_6to13_nocibersort <- as_tibble(dds_np_neg_0to5_vs_6to13_nocibersort)
res2_np_neg_0to5_vs_6to13_nocibersort <- res_np_neg_0to5_vs_6to13_nocibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_neg_0to5_vs_6to13_nocibersort <- deframe(res2_np_neg_0to5_vs_6to13_nocibersort)
fgsea_np_neg_0to5_vs_6to13_nocibersort <- fgsea(pathways=modules_61, stats=ranks_np_neg_0to5_vs_6to13_nocibersort)
fgsea_np_neg_0to5_vs_6to13_nocibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_neg_0to5_vs_6to13_nocibersort))) {fgsea_np_neg_0to5_vs_6to13_nocibersort[i,9]<-(length(unlist(fgsea_np_neg_0to5_vs_6to13_nocibersort[i,8])))/fgsea_np_neg_0to5_vs_6to13_nocibersort[i,7]}
fgsea_np_neg_0to5_vs_6to13_nocibersort <- data.frame(fgsea_np_neg_0to5_vs_6to13_nocibersort %>% 
                                        mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_neg_0to5_vs_6to13_nocibersort$leadingEdge <- NULL
fgsea_np_neg_0to5_vs_6to13_nocibersort$Sample <- "0-5 vs. 6-13 yr"
write_xlsx(fgsea_np_neg_0to5_vs_6to13_nocibersort, "Statistical_Analyses/1_COVID_Neg_by_Age/modules_np_neg_0to5_vs_6to13_nocibersort.xlsx")

res_np_neg_0to5_vs_14to20_nocibersort <- as_tibble(dds_np_neg_0to5_vs_14to20_nocibersort)
res2_np_neg_0to5_vs_14to20_nocibersort <- res_np_neg_0to5_vs_14to20_nocibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_neg_0to5_vs_14to20_nocibersort <- deframe(res2_np_neg_0to5_vs_14to20_nocibersort)
fgsea_np_neg_0to5_vs_14to20_nocibersort <- fgsea(pathways=modules_61, stats=ranks_np_neg_0to5_vs_14to20_nocibersort)
fgsea_np_neg_0to5_vs_14to20_nocibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_neg_0to5_vs_14to20_nocibersort))) {fgsea_np_neg_0to5_vs_14to20_nocibersort[i,9]<-(length(unlist(fgsea_np_neg_0to5_vs_14to20_nocibersort[i,8])))/fgsea_np_neg_0to5_vs_14to20_nocibersort[i,7]}
fgsea_np_neg_0to5_vs_14to20_nocibersort <- data.frame(fgsea_np_neg_0to5_vs_14to20_nocibersort %>% 
                                         mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_neg_0to5_vs_14to20_nocibersort$leadingEdge <- NULL
fgsea_np_neg_0to5_vs_14to20_nocibersort$Sample <- "0-5 vs. 14-20 yr"
write_xlsx(fgsea_np_neg_0to5_vs_14to20_nocibersort, "Statistical_Analyses/1_COVID_Neg_by_Age/modules_np_neg_0to5_vs_14to20_nocibersort.xlsx")

res_np_neg_6to13_vs_14to20_nocibersort <- as_tibble(dds_np_neg_6to13_vs_14to20_nocibersort)
res2_np_neg_6to13_vs_14to20_nocibersort <- res_np_neg_6to13_vs_14to20_nocibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_neg_6to13_vs_14to20_nocibersort <- deframe(res2_np_neg_6to13_vs_14to20_nocibersort)
fgsea_np_neg_6to13_vs_14to20_nocibersort <- fgsea(pathways=modules_61, stats=ranks_np_neg_6to13_vs_14to20_nocibersort)
fgsea_np_neg_6to13_vs_14to20_nocibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_neg_6to13_vs_14to20_nocibersort))) {fgsea_np_neg_6to13_vs_14to20_nocibersort[i,9]<-(length(unlist(fgsea_np_neg_6to13_vs_14to20_nocibersort[i,8])))/fgsea_np_neg_6to13_vs_14to20_nocibersort[i,7]}
fgsea_np_neg_6to13_vs_14to20_nocibersort <- data.frame(fgsea_np_neg_6to13_vs_14to20_nocibersort %>% 
                                          mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_neg_6to13_vs_14to20_nocibersort$Sample <- "6-13 vs. 14-20 yr"
fgsea_np_neg_6to13_vs_14to20_nocibersort$leadingEdge <- NULL
write_xlsx(fgsea_np_neg_6to13_vs_14to20_nocibersort, "Statistical_Analyses/1_COVID_Neg_by_Age/modules_np_neg_6to13_vs_14to20_nocibersort.xlsx")

#####################################################################
### DIFFERENCES IN PAX MODULE EXPRESSION BY AGE IN COVID-NEGATIVE ###
####### ANALYSES WITH ADJUSTMENT FOR IMPUTED CELL PROPORTIONS #######
#####################################################################

res_pax_neg_0to5_vs_6to13_cibersort <- as_tibble(dds_pax_neg_0to5_vs_6to13_cibersort)
res2_pax_neg_0to5_vs_6to13_cibersort <- res_pax_neg_0to5_vs_6to13_cibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_neg_0to5_vs_6to13_cibersort <- deframe(res2_pax_neg_0to5_vs_6to13_cibersort)
fgsea_pax_neg_0to5_vs_6to13_cibersort <- fgsea(pathways=modules_61, stats=ranks_pax_neg_0to5_vs_6to13_cibersort)
fgsea_pax_neg_0to5_vs_6to13_cibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_neg_0to5_vs_6to13_cibersort))) {fgsea_pax_neg_0to5_vs_6to13_cibersort[i,9]<-(length(unlist(fgsea_pax_neg_0to5_vs_6to13_cibersort[i,8])))/fgsea_pax_neg_0to5_vs_6to13_cibersort[i,7]}
fgsea_pax_neg_0to5_vs_6to13_cibersort <- data.frame(fgsea_pax_neg_0to5_vs_6to13_cibersort %>% 
                                         mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_neg_0to5_vs_6to13_cibersort$leadingEdge <- NULL
fgsea_pax_neg_0to5_vs_6to13_cibersort$Sample <- "0-5 vs. 6-13 yr"
write_xlsx(fgsea_pax_neg_0to5_vs_6to13_cibersort, "Statistical_Analyses/1_COVID_Neg_by_Age/modules_pax_neg_0to5_vs_6to13_cibersort.xlsx")

res_pax_neg_0to5_vs_14to20_cibersort <- as_tibble(dds_pax_neg_0to5_vs_14to20_cibersort)
res2_pax_neg_0to5_vs_14to20_cibersort <- res_pax_neg_0to5_vs_14to20_cibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_neg_0to5_vs_14to20_cibersort <- deframe(res2_pax_neg_0to5_vs_14to20_cibersort)
fgsea_pax_neg_0to5_vs_14to20_cibersort <- fgsea(pathways=modules_61, stats=ranks_pax_neg_0to5_vs_14to20_cibersort)
fgsea_pax_neg_0to5_vs_14to20_cibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_neg_0to5_vs_14to20_cibersort))) {fgsea_pax_neg_0to5_vs_14to20_cibersort[i,9]<-(length(unlist(fgsea_pax_neg_0to5_vs_14to20_cibersort[i,8])))/fgsea_pax_neg_0to5_vs_14to20_cibersort[i,7]}
fgsea_pax_neg_0to5_vs_14to20_cibersort <- data.frame(fgsea_pax_neg_0to5_vs_14to20_cibersort %>% 
                                          mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_neg_0to5_vs_14to20_cibersort$leadingEdge <- NULL
fgsea_pax_neg_0to5_vs_14to20_cibersort$Sample <- "0-5 vs. 14-20 yr"
write_xlsx(fgsea_pax_neg_0to5_vs_14to20_cibersort, "Statistical_Analyses/1_COVID_Neg_by_Age/modules_pax_neg_0to5_vs_14to20_cibersort.xlsx")

res_pax_neg_6to13_vs_14to20_cibersort <- as_tibble(dds_pax_neg_6to13_vs_14to20_cibersort)
res2_pax_neg_6to13_vs_14to20_cibersort <- res_pax_neg_6to13_vs_14to20_cibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_neg_6to13_vs_14to20_cibersort <- deframe(res2_pax_neg_6to13_vs_14to20_cibersort)
fgsea_pax_neg_6to13_vs_14to20_cibersort <- fgsea(pathways=modules_61, stats=ranks_pax_neg_6to13_vs_14to20_cibersort)
fgsea_pax_neg_6to13_vs_14to20_cibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_neg_6to13_vs_14to20_cibersort))) {fgsea_pax_neg_6to13_vs_14to20_cibersort[i,9]<-(length(unlist(fgsea_pax_neg_6to13_vs_14to20_cibersort[i,8])))/fgsea_pax_neg_6to13_vs_14to20_cibersort[i,7]}
fgsea_pax_neg_6to13_vs_14to20_cibersort <- data.frame(fgsea_pax_neg_6to13_vs_14to20_cibersort %>% 
                                           mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_neg_6to13_vs_14to20_cibersort$leadingEdge <- NULL
fgsea_pax_neg_6to13_vs_14to20_cibersort$Sample <- "6-13 vs. 14-20 yr"
write_xlsx(fgsea_pax_neg_6to13_vs_14to20_cibersort, "Statistical_Analyses/1_COVID_Neg_by_Age/modules_pax_neg_6to13_vs_14to20_cibersort.xlsx")

res_pax_neg_0to5_vs_adult_cibersort <- as_tibble(dds_pax_neg_0to5_vs_adult_cibersort)
res2_pax_neg_0to5_vs_adult_cibersort <- res_pax_neg_0to5_vs_adult_cibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_neg_0to5_vs_adult_cibersort <- deframe(res2_pax_neg_0to5_vs_adult_cibersort)
fgsea_pax_neg_0to5_vs_adult_cibersort <- fgsea(pathways=modules_61, stats=ranks_pax_neg_0to5_vs_adult_cibersort)
fgsea_pax_neg_0to5_vs_adult_cibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_neg_0to5_vs_adult_cibersort))) {fgsea_pax_neg_0to5_vs_adult_cibersort[i,9]<-(length(unlist(fgsea_pax_neg_0to5_vs_adult_cibersort[i,8])))/fgsea_pax_neg_0to5_vs_adult_cibersort[i,7]}
fgsea_pax_neg_0to5_vs_adult_cibersort <- data.frame(fgsea_pax_neg_0to5_vs_adult_cibersort %>% 
                                          mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_neg_0to5_vs_adult_cibersort$leadingEdge <- NULL
fgsea_pax_neg_0to5_vs_adult_cibersort$Sample <- "0-5 yr vs. Adult"
write_xlsx(fgsea_pax_neg_0to5_vs_adult_cibersort, "Statistical_Analyses/1_COVID_Neg_by_Age/modules_pax_neg_0to5_vs_adult_cibersort.xlsx")

res_pax_neg_6to13_vs_adult_cibersort <- as_tibble(dds_pax_neg_6to13_vs_adult_cibersort)
res2_pax_neg_6to13_vs_adult_cibersort <- res_pax_neg_6to13_vs_adult_cibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_neg_6to13_vs_adult_cibersort <- deframe(res2_pax_neg_6to13_vs_adult_cibersort)
fgsea_pax_neg_6to13_vs_adult_cibersort <- fgsea(pathways=modules_61, stats=ranks_pax_neg_6to13_vs_adult_cibersort)
fgsea_pax_neg_6to13_vs_adult_cibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_neg_6to13_vs_adult_cibersort))) {fgsea_pax_neg_6to13_vs_adult_cibersort[i,9]<-(length(unlist(fgsea_pax_neg_6to13_vs_adult_cibersort[i,8])))/fgsea_pax_neg_6to13_vs_adult_cibersort[i,7]}
fgsea_pax_neg_6to13_vs_adult_cibersort <- data.frame(fgsea_pax_neg_6to13_vs_adult_cibersort %>% 
                                           mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_neg_6to13_vs_adult_cibersort$leadingEdge <- NULL
fgsea_pax_neg_6to13_vs_adult_cibersort$Sample <- "6-13 yr vs. Adult"
write_xlsx(fgsea_pax_neg_6to13_vs_adult_cibersort, "Statistical_Analyses/1_COVID_Neg_by_Age/modules_pax_neg_6to13_vs_adult_cibersort.xlsx")

res_pax_neg_14to20_vs_adult_cibersort <- as_tibble(dds_pax_neg_14to20_vs_adult_cibersort)
res2_pax_neg_14to20_vs_adult_cibersort <- res_pax_neg_14to20_vs_adult_cibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_neg_14to20_vs_adult_cibersort <- deframe(res2_pax_neg_14to20_vs_adult_cibersort)
fgsea_pax_neg_14to20_vs_adult_cibersort <- fgsea(pathways=modules_61, stats=ranks_pax_neg_14to20_vs_adult_cibersort)
fgsea_pax_neg_14to20_vs_adult_cibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_neg_14to20_vs_adult_cibersort))) {fgsea_pax_neg_14to20_vs_adult_cibersort[i,9]<-(length(unlist(fgsea_pax_neg_14to20_vs_adult_cibersort[i,8])))/fgsea_pax_neg_14to20_vs_adult_cibersort[i,7]}
fgsea_pax_neg_14to20_vs_adult_cibersort <- data.frame(fgsea_pax_neg_14to20_vs_adult_cibersort %>% 
                                          mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_neg_14to20_vs_adult_cibersort$leadingEdge <- NULL
fgsea_pax_neg_14to20_vs_adult_cibersort$Sample <- "14-20 yr vs. Adult"
write_xlsx(fgsea_pax_neg_14to20_vs_adult_cibersort, "Statistical_Analyses/1_COVID_Neg_by_Age/modules_pax_neg_14to20_vs_adult_cibersort.xlsx")

####################################################
### DIFFERENCES IN NP EXPRESSION BY COVID STATUS ###
####################################################

res_np_pos_neg_nocibersort <- as_tibble(dds_np_pos_neg_nocibersort)
res2_np_pos_neg_nocibersort <- res_np_pos_neg_nocibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_pos_neg_nocibersort <- deframe(res2_np_pos_neg_nocibersort)
fgsea_np_pos_neg_nocibersort <- fgsea(pathways=modules_61, stats=ranks_np_pos_neg_nocibersort)
fgsea_np_pos_neg_nocibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_pos_neg_nocibersort))) {fgsea_np_pos_neg_nocibersort[i,9]<-(length(unlist(fgsea_np_pos_neg_nocibersort[i,8])))/fgsea_np_pos_neg_nocibersort[i,7]}
fgsea_np_pos_neg_nocibersort <- data.frame(fgsea_np_pos_neg_nocibersort %>% 
                                         mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_pos_neg_nocibersort$leadingEdge <- NULL
write_xlsx(fgsea_np_pos_neg_nocibersort, "Statistical_Analyses/2_COVID_Pos_vs_Neg/modules_np_pos_neg_nocibersort.xlsx")

##########################################################
### DIFFERENCES IN PAX GENE EXPRESSION BY COVID STATUS ###
# ANALYSES WITH ADJUSTMENT FOR IMPUTED CELL PROPORTIONS ##
##########################################################

res_pax_pos_neg_cibersort <- as_tibble(dds_pax_pos_neg_cibersort)
res2_pax_pos_neg_cibersort <- res_pax_pos_neg_cibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_pos_neg_cibersort <- deframe(res2_pax_pos_neg_cibersort)
fgsea_pax_pos_neg_cibersort <- fgsea(pathways=modules_61, stats=ranks_pax_pos_neg_cibersort)
fgsea_pax_pos_neg_cibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_pos_neg_cibersort))) {fgsea_pax_pos_neg_cibersort[i,9]<-(length(unlist(fgsea_pax_pos_neg_cibersort[i,8])))/fgsea_pax_pos_neg_cibersort[i,7]}
fgsea_pax_pos_neg_cibersort <- data.frame(fgsea_pax_pos_neg_cibersort %>% 
                                 mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_pos_neg_cibersort$leadingEdge <- NULL
write_xlsx(fgsea_pax_pos_neg_cibersort, "Statistical_Analyses/2_COVID_Pos_vs_Neg/modules_pax_pos_neg_cibersort.xlsx")

##############################################################
### COMPARING NP SAMPLES BY COVID STATUS IN AGE CATEGORIES ###
##############################################################

res_np_0to5_pos_neg_nocibersort <- as_tibble(dds_np_0to5_pos_neg_nocibersort)
res2_np_0to5_pos_neg_nocibersort <- res_np_0to5_pos_neg_nocibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_0to5_pos_neg_nocibersort <- deframe(res2_np_0to5_pos_neg_nocibersort)
fgsea_np_0to5_pos_neg_nocibersort <- fgsea(pathways=modules_61, stats=ranks_np_0to5_pos_neg_nocibersort)
fgsea_np_0to5_pos_neg_nocibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_0to5_pos_neg_nocibersort))) {fgsea_np_0to5_pos_neg_nocibersort[i,9]<-(length(unlist(fgsea_np_0to5_pos_neg_nocibersort[i,8])))/fgsea_np_0to5_pos_neg_nocibersort[i,7]}
fgsea_np_0to5_pos_neg_nocibersort <- data.frame(fgsea_np_0to5_pos_neg_nocibersort %>% 
                                        mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_0to5_pos_neg_nocibersort$leadingEdge <- NULL
fgsea_np_0to5_pos_neg_nocibersort$Sample <- "0-5 yr"
write_xlsx(fgsea_np_0to5_pos_neg_nocibersort, "Statistical_Analyses/3_COVID_Pos_by_Age/modules_np_0to5_pos_neg_nocibersort.xlsx")

res_np_6to13_pos_neg_nocibersort <- as_tibble(dds_np_6to13_pos_neg_nocibersort)
res2_np_6to13_pos_neg_nocibersort <- res_np_6to13_pos_neg_nocibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_6to13_pos_neg_nocibersort <- deframe(res2_np_6to13_pos_neg_nocibersort)
fgsea_np_6to13_pos_neg_nocibersort <- fgsea(pathways=modules_61, stats=ranks_np_6to13_pos_neg_nocibersort)
fgsea_np_6to13_pos_neg_nocibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_6to13_pos_neg_nocibersort))) {fgsea_np_6to13_pos_neg_nocibersort[i,9]<-(length(unlist(fgsea_np_6to13_pos_neg_nocibersort[i,8])))/fgsea_np_6to13_pos_neg_nocibersort[i,7]}
fgsea_np_6to13_pos_neg_nocibersort <- data.frame(fgsea_np_6to13_pos_neg_nocibersort %>% 
                                      mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_6to13_pos_neg_nocibersort$leadingEdge <- NULL
fgsea_np_6to13_pos_neg_nocibersort$Sample <- "6-13 yr"
write_xlsx(fgsea_np_6to13_pos_neg_nocibersort, "Statistical_Analyses/3_COVID_Pos_by_Age/modules_np_6to13_pos_neg_nocibersort.xlsx")

res_np_14to20_pos_neg_nocibersort <- as_tibble(dds_np_14to20_pos_neg_nocibersort)
res2_np_14to20_pos_neg_nocibersort <- res_np_14to20_pos_neg_nocibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_14to20_pos_neg_nocibersort <- deframe(res2_np_14to20_pos_neg_nocibersort)
fgsea_np_14to20_pos_neg_nocibersort <- fgsea(pathways=modules_61, stats=ranks_np_14to20_pos_neg_nocibersort)
fgsea_np_14to20_pos_neg_nocibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_14to20_pos_neg_nocibersort))) {fgsea_np_14to20_pos_neg_nocibersort[i,9]<-(length(unlist(fgsea_np_14to20_pos_neg_nocibersort[i,8])))/fgsea_np_14to20_pos_neg_nocibersort[i,7]}
fgsea_np_14to20_pos_neg_nocibersort <- data.frame(fgsea_np_14to20_pos_neg_nocibersort %>% 
                                       mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_14to20_pos_neg_nocibersort$leadingEdge <- NULL
fgsea_np_14to20_pos_neg_nocibersort$Sample <- "14-20 yr"
write_xlsx(fgsea_np_14to20_pos_neg_nocibersort, "Statistical_Analyses/3_COVID_Pos_by_Age/modules_np_14to20_pos_neg_nocibersort.xlsx")

###############################################################
### COMPARING PAX SAMPLES BY COVID STATUS IN AGE CATEGORIES ###
#### ANALYSES WITH ADJUSTMENT FOR IMPUTED CELL PROPORTIONS ####
###############################################################

res_pax_0to5_pos_neg_cibersort <- as_tibble(dds_pax_0to5_pos_neg_cibersort)
res2_pax_0to5_pos_neg_cibersort <- res_pax_0to5_pos_neg_cibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_0to5_pos_neg_cibersort <- deframe(res2_pax_0to5_pos_neg_cibersort)
fgsea_pax_0to5_pos_neg_cibersort <- fgsea(pathways=modules_61, stats=ranks_pax_0to5_pos_neg_cibersort)
fgsea_pax_0to5_pos_neg_cibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_0to5_pos_neg_cibersort))) {fgsea_pax_0to5_pos_neg_cibersort[i,9]<-(length(unlist(fgsea_pax_0to5_pos_neg_cibersort[i,8])))/fgsea_pax_0to5_pos_neg_cibersort[i,7]}
fgsea_pax_0to5_pos_neg_cibersort <- data.frame(fgsea_pax_0to5_pos_neg_cibersort %>% 
                                      mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_0to5_pos_neg_cibersort$leadingEdge <- NULL
fgsea_pax_0to5_pos_neg_cibersort$Sample <- "0-5 yr"
write_xlsx(fgsea_pax_0to5_pos_neg_cibersort, "Statistical_Analyses/3_COVID_Pos_by_Age/modules_pax_0to5_pos_neg_cibersort.xlsx")

res_pax_6to13_pos_neg_cibersort <- as_tibble(dds_pax_6to13_pos_neg_cibersort)
res2_pax_6to13_pos_neg_cibersort <- res_pax_6to13_pos_neg_cibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_6to13_pos_neg_cibersort <- deframe(res2_pax_6to13_pos_neg_cibersort)
fgsea_pax_6to13_pos_neg_cibersort <- fgsea(pathways=modules_61, stats=ranks_pax_6to13_pos_neg_cibersort)
fgsea_pax_6to13_pos_neg_cibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_6to13_pos_neg_cibersort))) {fgsea_pax_6to13_pos_neg_cibersort[i,9]<-(length(unlist(fgsea_pax_6to13_pos_neg_cibersort[i,8])))/fgsea_pax_6to13_pos_neg_cibersort[i,7]}
fgsea_pax_6to13_pos_neg_cibersort <- data.frame(fgsea_pax_6to13_pos_neg_cibersort %>% 
                                       mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_6to13_pos_neg_cibersort$leadingEdge <- NULL
fgsea_pax_6to13_pos_neg_cibersort$Sample <- "6-13 yr"
write_xlsx(fgsea_pax_6to13_pos_neg_cibersort, "Statistical_Analyses/3_COVID_Pos_by_Age/modules_pax_6to13_pos_neg_cibersort.xlsx")

res_pax_14to20_pos_neg_cibersort <- as_tibble(dds_pax_14to20_pos_neg_cibersort)
res2_pax_14to20_pos_neg_cibersort <- res_pax_14to20_pos_neg_cibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_14to20_pos_neg_cibersort <- deframe(res2_pax_14to20_pos_neg_cibersort)
fgsea_pax_14to20_pos_neg_cibersort <- fgsea(pathways=modules_61, stats=ranks_pax_14to20_pos_neg_cibersort)
fgsea_pax_14to20_pos_neg_cibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_14to20_pos_neg_cibersort))) {fgsea_pax_14to20_pos_neg_cibersort[i,9]<-(length(unlist(fgsea_pax_14to20_pos_neg_cibersort[i,8])))/fgsea_pax_14to20_pos_neg_cibersort[i,7]}
fgsea_pax_14to20_pos_neg_cibersort <- data.frame(fgsea_pax_14to20_pos_neg_cibersort %>% 
                                        mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_14to20_pos_neg_cibersort$leadingEdge <- NULL
fgsea_pax_14to20_pos_neg_cibersort$Sample <- "14-20 yr"
write_xlsx(fgsea_pax_14to20_pos_neg_cibersort, "Statistical_Analyses/3_COVID_Pos_by_Age/modules_pax_14to20_pos_neg_cibersort.xlsx")

res_pax_adult_pos_neg_cibersort <- as_tibble(dds_pax_adult_pos_neg_cibersort)
res2_pax_adult_pos_neg_cibersort <- res_pax_adult_pos_neg_cibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_adult_pos_neg_cibersort <- deframe(res2_pax_adult_pos_neg_cibersort)
fgsea_pax_adult_pos_neg_cibersort <- fgsea(pathways=modules_61, stats=ranks_pax_adult_pos_neg_cibersort)
fgsea_pax_adult_pos_neg_cibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_adult_pos_neg_cibersort))) {fgsea_pax_adult_pos_neg_cibersort[i,9]<-(length(unlist(fgsea_pax_adult_pos_neg_cibersort[i,8])))/fgsea_pax_adult_pos_neg_cibersort[i,7]}
fgsea_pax_adult_pos_neg_cibersort <- data.frame(fgsea_pax_adult_pos_neg_cibersort %>% 
                                         mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_adult_pos_neg_cibersort$leadingEdge <- NULL
fgsea_pax_adult_pos_neg_cibersort$Sample <- "Adult"
write_xlsx(fgsea_pax_adult_pos_neg_cibersort, "Statistical_Analyses/3_COVID_Pos_by_Age/modules_pax_adult_pos_neg_cibersort.xlsx")

###############################################################################
### MODELS FOR PATIENT CHARACTERISTICS & SYMPTOMS AMONG COVID+ - NP SAMPLES ###
###############################################################################

# Fever
res_np_fever_nocibersort <- as_tibble(dds_np_fever_nocibersort)
res2_np_fever_nocibersort <- res_np_fever_nocibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_fever_nocibersort <- deframe(res2_np_fever_nocibersort)
fgsea_np_fever_nocibersort <- fgsea(pathways=modules_61, stats=ranks_np_fever_nocibersort)
fgsea_np_fever_nocibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_fever_nocibersort))) {fgsea_np_fever_nocibersort[i,9]<-(length(unlist(fgsea_np_fever_nocibersort[i,8])))/fgsea_np_fever_nocibersort[i,7]}
fgsea_np_fever_nocibersort <- data.frame(fgsea_np_fever_nocibersort %>% 
                                mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_fever_nocibersort$leadingEdge <- NULL
write_xlsx(fgsea_np_fever_nocibersort, "Statistical_Analyses/4_Symptoms/modules_np_fever_nocibersort.xlsx")

# Cough
res_np_cough_nocibersort <- as_tibble(dds_np_cough_nocibersort)
res2_np_cough_nocibersort <- res_np_cough_nocibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_cough_nocibersort <- deframe(res2_np_cough_nocibersort)
fgsea_np_cough_nocibersort <- fgsea(pathways=modules_61, stats=ranks_np_cough_nocibersort)
fgsea_np_cough_nocibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_cough_nocibersort))) {fgsea_np_cough_nocibersort[i,9]<-(length(unlist(fgsea_np_cough_nocibersort[i,8])))/fgsea_np_cough_nocibersort[i,7]}
fgsea_np_cough_nocibersort <- data.frame(fgsea_np_cough_nocibersort %>% 
                                mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_cough_nocibersort$leadingEdge <- NULL
write_xlsx(fgsea_np_cough_nocibersort, "Statistical_Analyses/4_Symptoms/modules_np_cough_nocibersort.xlsx")

# Rhinorrhea
res_np_rhinorrhea_nocibersort <- as_tibble(dds_np_rhinorrhea_nocibersort)
res2_np_rhinorrhea_nocibersort <- res_np_rhinorrhea_nocibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_rhinorrhea_nocibersort <- deframe(res2_np_rhinorrhea_nocibersort)
fgsea_np_rhinorrhea_nocibersort <- fgsea(pathways=modules_61, stats=ranks_np_rhinorrhea_nocibersort)
fgsea_np_rhinorrhea_nocibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_rhinorrhea_nocibersort))) {fgsea_np_rhinorrhea_nocibersort[i,9]<-(length(unlist(fgsea_np_rhinorrhea_nocibersort[i,8])))/fgsea_np_rhinorrhea_nocibersort[i,7]}
fgsea_np_rhinorrhea_nocibersort <- data.frame(fgsea_np_rhinorrhea_nocibersort %>% 
                                     mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_rhinorrhea_nocibersort$leadingEdge <- NULL
write_xlsx(fgsea_np_rhinorrhea_nocibersort, "Statistical_Analyses/4_Symptoms/modules_np_rhinorrhea_nocibersort.xlsx")

# Nasal congestion
res_np_congestion_nocibersort <- as_tibble(dds_np_congestion_nocibersort)
res2_np_congestion_nocibersort <- res_np_congestion_nocibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_congestion_nocibersort <- deframe(res2_np_congestion_nocibersort)
fgsea_np_congestion_nocibersort <- fgsea(pathways=modules_61, stats=ranks_np_congestion_nocibersort)
fgsea_np_congestion_nocibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_congestion_nocibersort))) {fgsea_np_congestion_nocibersort[i,9]<-(length(unlist(fgsea_np_congestion_nocibersort[i,8])))/fgsea_np_congestion_nocibersort[i,7]}
fgsea_np_congestion_nocibersort <- data.frame(fgsea_np_congestion_nocibersort %>% 
                                     mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_congestion_nocibersort$leadingEdge <- NULL
write_xlsx(fgsea_np_congestion_nocibersort, "Statistical_Analyses/4_Symptoms/modules_np_congestion_nocibersort.xlsx")

# Headache
res_np_headache_nocibersort <- as_tibble(dds_np_headache_nocibersort)
res2_np_headache_nocibersort <- res_np_headache_nocibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_headache_nocibersort <- deframe(res2_np_headache_nocibersort)
fgsea_np_headache_nocibersort <- fgsea(pathways=modules_61, stats=ranks_np_headache_nocibersort)
fgsea_np_headache_nocibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_headache_nocibersort))) {fgsea_np_headache_nocibersort[i,9]<-(length(unlist(fgsea_np_headache_nocibersort[i,8])))/fgsea_np_headache_nocibersort[i,7]}
fgsea_np_headache_nocibersort <- data.frame(fgsea_np_headache_nocibersort %>% 
                               mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_headache_nocibersort$leadingEdge <- NULL
write_xlsx(fgsea_np_headache_nocibersort, "Statistical_Analyses/4_Symptoms/modules_np_headache_nocibersort.xlsx")

# Abdominal pain
res_np_abd_pain_nocibersort <- as_tibble(dds_np_abd_pain_nocibersort)
res2_np_abd_pain_nocibersort <- res_np_abd_pain_nocibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_abd_pain_nocibersort <- deframe(res2_np_abd_pain_nocibersort)
fgsea_np_abd_pain_nocibersort <- fgsea(pathways=modules_61, stats=ranks_np_abd_pain_nocibersort)
fgsea_np_abd_pain_nocibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_abd_pain_nocibersort))) {fgsea_np_abd_pain_nocibersort[i,9]<-(length(unlist(fgsea_np_abd_pain_nocibersort[i,8])))/fgsea_np_abd_pain_nocibersort[i,7]}
fgsea_np_abd_pain_nocibersort <- data.frame(fgsea_np_abd_pain_nocibersort %>% 
                               mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_abd_pain_nocibersort$leadingEdge <- NULL
write_xlsx(fgsea_np_abd_pain_nocibersort, "Statistical_Analyses/4_Symptoms/modules_np_abd_pain_nocibersort.xlsx")

# Anosmia
res_np_anosmia_nocibersort <- as_tibble(dds_np_anosmia_nocibersort)
res2_np_anosmia_nocibersort <- res_np_anosmia_nocibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_anosmia_nocibersort <- deframe(res2_np_anosmia_nocibersort)
fgsea_np_anosmia_nocibersort <- fgsea(pathways=modules_61, stats=ranks_np_anosmia_nocibersort)
fgsea_np_anosmia_nocibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_anosmia_nocibersort))) {fgsea_np_anosmia_nocibersort[i,9]<-(length(unlist(fgsea_np_anosmia_nocibersort[i,8])))/fgsea_np_anosmia_nocibersort[i,7]}
fgsea_np_anosmia_nocibersort <- data.frame(fgsea_np_anosmia_nocibersort %>% 
                               mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_anosmia_nocibersort$leadingEdge <- NULL
write_xlsx(fgsea_np_anosmia_nocibersort, "Statistical_Analyses/4_Symptoms/modules_np_anosmia_nocibersort.xlsx")

# Dysgeusia
res_np_dysgeusia_nocibersort <- as_tibble(dds_np_dysgeusia_nocibersort)
res2_np_dysgeusia_nocibersort <- res_np_dysgeusia_nocibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_dysgeusia_nocibersort <- deframe(res2_np_dysgeusia_nocibersort)
fgsea_np_dysgeusia_nocibersort <- fgsea(pathways=modules_61, stats=ranks_np_dysgeusia_nocibersort)
fgsea_np_dysgeusia_nocibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_dysgeusia_nocibersort))) {fgsea_np_dysgeusia_nocibersort[i,9]<-(length(unlist(fgsea_np_dysgeusia_nocibersort[i,8])))/fgsea_np_dysgeusia_nocibersort[i,7]}
fgsea_np_dysgeusia_nocibersort <- data.frame(fgsea_np_dysgeusia_nocibersort %>% 
                               mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_dysgeusia_nocibersort$leadingEdge <- NULL
write_xlsx(fgsea_np_dysgeusia_nocibersort, "Statistical_Analyses/4_Symptoms/modules_np_dysgeusia_nocibersort.xlsx")

# Myalgias
res_np_myalgias_nocibersort <- as_tibble(dds_np_myalgias_nocibersort)
res2_np_myalgias_nocibersort <- res_np_myalgias_nocibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_myalgias_nocibersort <- deframe(res2_np_myalgias_nocibersort)
fgsea_np_myalgias_nocibersort <- fgsea(pathways=modules_61, stats=ranks_np_myalgias_nocibersort)
fgsea_np_myalgias_nocibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_myalgias_nocibersort))) {fgsea_np_myalgias_nocibersort[i,9]<-(length(unlist(fgsea_np_myalgias_nocibersort[i,8])))/fgsea_np_myalgias_nocibersort[i,7]}
fgsea_np_myalgias_nocibersort <- data.frame(fgsea_np_myalgias_nocibersort %>% 
                               mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_myalgias_nocibersort$leadingEdge <- NULL
write_xlsx(fgsea_np_myalgias_nocibersort, "Statistical_Analyses/4_Symptoms/modules_np_myalgias_nocibersort.xlsx")

##############################################################################
## MODELS FOR PATIENT CHARACTERISTICS & SYMPTOMS AMONG COVID+ - PAX SAMPLES ##
########## ANALYSES WITH ADJUSTMENT FOR IMPUTED CELL PROPORTIONS #############
##############################################################################

# Fever
res_pax_fever_cibersort <- as_tibble(dds_pax_fever_cibersort)
res2_pax_fever_cibersort <- res_pax_fever_cibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_fever_cibersort <- deframe(res2_pax_fever_cibersort)
fgsea_pax_fever_cibersort <- fgsea(pathways=modules_61, stats=ranks_pax_fever_cibersort)
fgsea_pax_fever_cibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_fever_cibersort))) {fgsea_pax_fever_cibersort[i,9]<-(length(unlist(fgsea_pax_fever_cibersort[i,8])))/fgsea_pax_fever_cibersort[i,7]}
fgsea_pax_fever_cibersort <- data.frame(fgsea_pax_fever_cibersort %>% 
                                            mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_fever_cibersort$leadingEdge <- NULL
write_xlsx(fgsea_pax_fever_cibersort, "Statistical_Analyses/4_Symptoms/modules_pax_fever_cibersort.xlsx")

# Cough
res_pax_cough_cibersort <- as_tibble(dds_pax_cough_cibersort)
res2_pax_cough_cibersort <- res_pax_cough_cibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_cough_cibersort <- deframe(res2_pax_cough_cibersort)
fgsea_pax_cough_cibersort <- fgsea(pathways=modules_61, stats=ranks_pax_cough_cibersort)
fgsea_pax_cough_cibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_cough_cibersort))) {fgsea_pax_cough_cibersort[i,9]<-(length(unlist(fgsea_pax_cough_cibersort[i,8])))/fgsea_pax_cough_cibersort[i,7]}
fgsea_pax_cough_cibersort <- data.frame(fgsea_pax_cough_cibersort %>% 
                                            mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_cough_cibersort$leadingEdge <- NULL
write_xlsx(fgsea_pax_cough_cibersort, "Statistical_Analyses/4_Symptoms/modules_pax_cough_cibersort.xlsx")

# Rhinorrhea
res_pax_rhinorrhea_cibersort <- as_tibble(dds_pax_rhinorrhea_cibersort)
res2_pax_rhinorrhea_cibersort <- res_pax_rhinorrhea_cibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_rhinorrhea_cibersort <- deframe(res2_pax_rhinorrhea_cibersort)
fgsea_pax_rhinorrhea_cibersort <- fgsea(pathways=modules_61, stats=ranks_pax_rhinorrhea_cibersort)
fgsea_pax_rhinorrhea_cibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_rhinorrhea_cibersort))) {fgsea_pax_rhinorrhea_cibersort[i,9]<-(length(unlist(fgsea_pax_rhinorrhea_cibersort[i,8])))/fgsea_pax_rhinorrhea_cibersort[i,7]}
fgsea_pax_rhinorrhea_cibersort <- data.frame(fgsea_pax_rhinorrhea_cibersort %>% 
                                                 mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_rhinorrhea_cibersort$leadingEdge <- NULL
write_xlsx(fgsea_pax_rhinorrhea_cibersort, "Statistical_Analyses/4_Symptoms/modules_pax_rhinorrhea_cibersort.xlsx")

# Nasal congestion
res_pax_congestion_cibersort <- as_tibble(dds_pax_congestion_cibersort)
res2_pax_congestion_cibersort <- res_pax_congestion_cibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_congestion_cibersort <- deframe(res2_pax_congestion_cibersort)
fgsea_pax_congestion_cibersort <- fgsea(pathways=modules_61, stats=ranks_pax_congestion_cibersort)
fgsea_pax_congestion_cibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_congestion_cibersort))) {fgsea_pax_congestion_cibersort[i,9]<-(length(unlist(fgsea_pax_congestion_cibersort[i,8])))/fgsea_pax_congestion_cibersort[i,7]}
fgsea_pax_congestion_cibersort <- data.frame(fgsea_pax_congestion_cibersort %>% 
                                                 mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_congestion_cibersort$leadingEdge <- NULL
write_xlsx(fgsea_pax_congestion_cibersort, "Statistical_Analyses/4_Symptoms/modules_pax_congestion_cibersort.xlsx")

# Headache
res_pax_headache_cibersort <- as_tibble(dds_pax_headache_cibersort)
res2_pax_headache_cibersort <- res_pax_headache_cibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_headache_cibersort <- deframe(res2_pax_headache_cibersort)
fgsea_pax_headache_cibersort <- fgsea(pathways=modules_61, stats=ranks_pax_headache_cibersort)
fgsea_pax_headache_cibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_headache_cibersort))) {fgsea_pax_headache_cibersort[i,9]<-(length(unlist(fgsea_pax_headache_cibersort[i,8])))/fgsea_pax_headache_cibersort[i,7]}
fgsea_pax_headache_cibersort <- data.frame(fgsea_pax_headache_cibersort %>% 
                                               mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_headache_cibersort$leadingEdge <- NULL
write_xlsx(fgsea_pax_headache_cibersort, "Statistical_Analyses/4_Symptoms/modules_pax_headache_cibersort.xlsx")

# Abdominal pain
res_pax_abd_pain_cibersort <- as_tibble(dds_pax_abd_pain_cibersort)
res2_pax_abd_pain_cibersort <- res_pax_abd_pain_cibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_abd_pain_cibersort <- deframe(res2_pax_abd_pain_cibersort)
fgsea_pax_abd_pain_cibersort <- fgsea(pathways=modules_61, stats=ranks_pax_abd_pain_cibersort)
fgsea_pax_abd_pain_cibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_abd_pain_cibersort))) {fgsea_pax_abd_pain_cibersort[i,9]<-(length(unlist(fgsea_pax_abd_pain_cibersort[i,8])))/fgsea_pax_abd_pain_cibersort[i,7]}
fgsea_pax_abd_pain_cibersort <- data.frame(fgsea_pax_abd_pain_cibersort %>% 
                                               mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_abd_pain_cibersort$leadingEdge <- NULL
write_xlsx(fgsea_pax_abd_pain_cibersort, "Statistical_Analyses/4_Symptoms/modules_pax_abd_pain_cibersort.xlsx")

# Anosmia
res_pax_anosmia_cibersort <- as_tibble(dds_pax_anosmia_cibersort)
res2_pax_anosmia_cibersort <- res_pax_anosmia_cibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_anosmia_cibersort <- deframe(res2_pax_anosmia_cibersort)
fgsea_pax_anosmia_cibersort <- fgsea(pathways=modules_61, stats=ranks_pax_anosmia_cibersort)
fgsea_pax_anosmia_cibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_anosmia_cibersort))) {fgsea_pax_anosmia_cibersort[i,9]<-(length(unlist(fgsea_pax_anosmia_cibersort[i,8])))/fgsea_pax_anosmia_cibersort[i,7]}
fgsea_pax_anosmia_cibersort <- data.frame(fgsea_pax_anosmia_cibersort %>% 
                                              mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_anosmia_cibersort$leadingEdge <- NULL
write_xlsx(fgsea_pax_anosmia_cibersort, "Statistical_Analyses/4_Symptoms/modules_pax_anosmia_cibersort.xlsx")

# Dysgeusia
res_pax_dysgeusia_cibersort <- as_tibble(dds_pax_dysgeusia_cibersort)
res2_pax_dysgeusia_cibersort <- res_pax_dysgeusia_cibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_dysgeusia_cibersort <- deframe(res2_pax_dysgeusia_cibersort)
fgsea_pax_dysgeusia_cibersort <- fgsea(pathways=modules_61, stats=ranks_pax_dysgeusia_cibersort)
fgsea_pax_dysgeusia_cibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_dysgeusia_cibersort))) {fgsea_pax_dysgeusia_cibersort[i,9]<-(length(unlist(fgsea_pax_dysgeusia_cibersort[i,8])))/fgsea_pax_dysgeusia_cibersort[i,7]}
fgsea_pax_dysgeusia_cibersort <- data.frame(fgsea_pax_dysgeusia_cibersort %>% 
                                                mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_dysgeusia_cibersort$leadingEdge <- NULL
write_xlsx(fgsea_pax_dysgeusia_cibersort, "Statistical_Analyses/4_Symptoms/modules_pax_dysgeusia_cibersort.xlsx")

# Myalgias
res_pax_myalgias_cibersort <- as_tibble(dds_pax_myalgias_cibersort)
res2_pax_myalgias_cibersort <- res_pax_myalgias_cibersort %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_myalgias_cibersort <- deframe(res2_pax_myalgias_cibersort)
fgsea_pax_myalgias_cibersort <- fgsea(pathways=modules_61, stats=ranks_pax_myalgias_cibersort)
fgsea_pax_myalgias_cibersort$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_myalgias_cibersort))) {fgsea_pax_myalgias_cibersort[i,9]<-(length(unlist(fgsea_pax_myalgias_cibersort[i,8])))/fgsea_pax_myalgias_cibersort[i,7]}
fgsea_pax_myalgias_cibersort <- data.frame(fgsea_pax_myalgias_cibersort %>% 
                                               mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_myalgias_cibersort$leadingEdge <- NULL
write_xlsx(fgsea_pax_myalgias_cibersort, "Statistical_Analyses/4_Symptoms/modules_pax_myalgias_cibersort.xlsx")