# BRAVE Kids RNA Sequencing Analysis
# Aditya Mohan (MD/PhD candidate)  / Matthew Kelly, MD, MPH 
# Gene Set Enrichment Analysis using FGSEA
# Last update: Oct. 5, 2023

remove(list=ls())
setwd("__________________________________") 
options(rstudio.help.showDataPreview = FALSE)
set.seed(1234)
version

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

# Upload module information 
modules_61 <- gmtPathways("Statistical_Analyses/modules_61.gmt")

# Upload DESeq2 output files
dds_np_Neg_0to5_6to13 <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_np_Neg_0to5_6to13.csv")
dds_np_Neg_0to5_14to20 <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_np_Neg_0to5_14to20.csv")
dds_np_Neg_6to13_14to20 <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_np_Neg_6to13_14to20.csv")
dds_pax_Neg_0to5_6to13 <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_pax_Neg_0to5_6to13.csv")
dds_pax_Neg_0to5_14to20 <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_pax_Neg_0to5_14to20.csv")
dds_pax_Neg_6to13_14to20 <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_pax_Neg_6to13_14to20.csv")
dds_pax_Neg_0to5_Adult <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_pax_Neg_0to5_Adult.csv")
dds_pax_Neg_6to13_Adult <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_pax_Neg_6to13_Adult.csv")
dds_pax_Neg_14to20_Adult <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_pax_Neg_14to20_Adult.csv")
dds_pax_Neg_Child_Adult <- read.csv("Statistical_Analyses/1_COVID_Neg_by_Age/genes_pax_Neg_Child_Adult.csv")
dds_np_Pos_Neg <- read.csv("Statistical_Analyses/2_COVID_Pos_vs_Neg/genes_np_Pos_Neg.csv")
dds_pax_Pos_Neg <- read.csv("Statistical_Analyses/2_COVID_Pos_vs_Neg/genes_pax_Pos_Neg.csv")
dds_np_0to5_Pos_Neg <- read.csv("Statistical_Analyses/3_COVID_Pos_by_Age/genes_np_0to5_Pos_Neg.csv")
dds_np_6to13_Pos_Neg <- read.csv("Statistical_Analyses/3_COVID_Pos_by_Age/genes_np_6to13_Pos_Neg.csv")
dds_np_14to20_Pos_Neg <- read.csv("Statistical_Analyses/3_COVID_Pos_by_Age/genes_np_14to20_Pos_Neg.csv")
dds_pax_0to5_Pos_Neg <- read.csv("Statistical_Analyses/3_COVID_Pos_by_Age/genes_pax_0to5_Pos_Neg.csv")
dds_pax_6to13_Pos_Neg <- read.csv("Statistical_Analyses/3_COVID_Pos_by_Age/genes_pax_6to13_Pos_Neg.csv")
dds_pax_14to20_Pos_Neg <- read.csv("Statistical_Analyses/3_COVID_Pos_by_Age/genes_pax_14to20_Pos_Neg.csv")
dds_pax_adult_Pos_Neg <- read.csv("Statistical_Analyses/3_COVID_Pos_by_Age/genes_pax_adult_Pos_Neg.csv")
dds_np_fever <- read.csv("Statistical_Analyses/4_Symptoms/genes_np_fever.csv")
dds_np_cough <- read.csv("Statistical_Analyses/4_Symptoms/genes_np_cough.csv")
dds_np_rhinorrhea <- read.csv("Statistical_Analyses/4_Symptoms/genes_np_rhinorrhea.csv")
dds_np_congestion <- read.csv("Statistical_Analyses/4_Symptoms/genes_np_congestion.csv")
dds_np_headache <- read.csv("Statistical_Analyses/4_Symptoms/genes_np_headache.csv")
dds_np_abd_pain <- read.csv("Statistical_Analyses/4_Symptoms/genes_np_abd_pain.csv")
dds_np_anosmia <- read.csv("Statistical_Analyses/4_Symptoms/genes_np_anosmia.csv")
dds_np_dysgeusia <- read.csv("Statistical_Analyses/4_Symptoms/genes_np_dysgeusia.csv")
dds_np_myalgias <- read.csv("Statistical_Analyses/4_Symptoms/genes_np_myalgias.csv")
dds_pax_fever <- read.csv("Statistical_Analyses/4_Symptoms/genes_pax_fever.csv")
dds_pax_cough <- read.csv("Statistical_Analyses/4_Symptoms/genes_pax_cough.csv")
dds_pax_rhinorrhea <- read.csv("Statistical_Analyses/4_Symptoms/genes_pax_rhinorrhea.csv")
dds_pax_congestion <- read.csv("Statistical_Analyses/4_Symptoms/genes_pax_congestion.csv")
dds_pax_headache <- read.csv("Statistical_Analyses/4_Symptoms/genes_pax_headache.csv")
dds_pax_abd_pain <- read.csv("Statistical_Analyses/4_Symptoms/genes_pax_abd_pain.csv")
dds_pax_anosmia <- read.csv("Statistical_Analyses/4_Symptoms/genes_pax_anosmia.csv")
dds_pax_dysgeusia <- read.csv("Statistical_Analyses/4_Symptoms/genes_pax_dysgeusia.csv")
dds_pax_myalgias <- read.csv("Statistical_Analyses/4_Symptoms/genes_pax_myalgias.csv")

#####################################################################
###  DIFFERENCES IN NP MODULE EXPRESSION BY AGE IN COVID NEGATIVE ###
#####################################################################

res_np_Neg_0to5_6to13 <- as_tibble(dds_np_Neg_0to5_6to13)
res2_np_Neg_0to5_6to13 <- res_np_Neg_0to5_6to13 %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_Neg_0to5_6to13 <- deframe(res2_np_Neg_0to5_6to13)
fgsea_np_Neg_0to5_6to13 <- fgsea(pathways=modules_61, stats=ranks_np_Neg_0to5_6to13)
fgsea_np_Neg_0to5_6to13$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_Neg_0to5_6to13))) {fgsea_np_Neg_0to5_6to13[i,9]<-(length(unlist(fgsea_np_Neg_0to5_6to13[i,8])))/fgsea_np_Neg_0to5_6to13[i,7]}
fgsea_np_Neg_0to5_6to13 <- data.frame(fgsea_np_Neg_0to5_6to13 %>% 
                                        mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_Neg_0to5_6to13$leadingEdge <- NULL
fgsea_np_Neg_0to5_6to13$Sample <- "0-5 vs. 6-13 yr"
write_xlsx(fgsea_np_Neg_0to5_6to13, "Statistical_Analyses/1_COVID_Neg_by_Age/modules_np_Neg_0to5_6to13.xlsx")

res_np_Neg_0to5_14to20 <- as_tibble(dds_np_Neg_0to5_14to20)
res2_np_Neg_0to5_14to20 <- res_np_Neg_0to5_14to20 %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_Neg_0to5_14to20 <- deframe(res2_np_Neg_0to5_14to20)
fgsea_np_Neg_0to5_14to20 <- fgsea(pathways=modules_61, stats=ranks_np_Neg_0to5_14to20)
fgsea_np_Neg_0to5_14to20$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_Neg_0to5_14to20))) {fgsea_np_Neg_0to5_14to20[i,9]<-(length(unlist(fgsea_np_Neg_0to5_14to20[i,8])))/fgsea_np_Neg_0to5_14to20[i,7]}
fgsea_np_Neg_0to5_14to20 <- data.frame(fgsea_np_Neg_0to5_14to20 %>% 
                                         mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_Neg_0to5_14to20$leadingEdge <- NULL
fgsea_np_Neg_0to5_14to20$Sample <- "0-5 vs. 14-20 yr"
write_xlsx(fgsea_np_Neg_0to5_14to20, "Statistical_Analyses/1_COVID_Neg_by_Age/modules_np_Neg_0to5_14to20.xlsx")

res_np_Neg_6to13_14to20 <- as_tibble(dds_np_Neg_6to13_14to20)
res2_np_Neg_6to13_14to20 <- res_np_Neg_6to13_14to20 %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_Neg_6to13_14to20 <- deframe(res2_np_Neg_6to13_14to20)
fgsea_np_Neg_6to13_14to20 <- fgsea(pathways=modules_61, stats=ranks_np_Neg_6to13_14to20)
fgsea_np_Neg_6to13_14to20$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_Neg_6to13_14to20))) {fgsea_np_Neg_6to13_14to20[i,9]<-(length(unlist(fgsea_np_Neg_6to13_14to20[i,8])))/fgsea_np_Neg_6to13_14to20[i,7]}
fgsea_np_Neg_6to13_14to20 <- data.frame(fgsea_np_Neg_6to13_14to20 %>% 
                                          mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_Neg_6to13_14to20$Sample <- "6-13 vs. 14-20 yr"
fgsea_np_Neg_6to13_14to20$leadingEdge <- NULL
write_xlsx(fgsea_np_Neg_6to13_14to20, "Statistical_Analyses/1_COVID_Neg_by_Age/modules_np_Neg_6to13_14to20.xlsx")

#####################################################################
### DIFFERENCES IN PAX MODULE EXPRESSION BY AGE IN COVID NEGATIVE ###
#####################################################################

res_pax_Neg_0to5_6to13 <- as_tibble(dds_pax_Neg_0to5_6to13)
res2_pax_Neg_0to5_6to13 <- res_pax_Neg_0to5_6to13 %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_Neg_0to5_6to13 <- deframe(res2_pax_Neg_0to5_6to13)
fgsea_pax_Neg_0to5_6to13 <- fgsea(pathways=modules_61, stats=ranks_pax_Neg_0to5_6to13)
fgsea_pax_Neg_0to5_6to13$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_Neg_0to5_6to13))) {fgsea_pax_Neg_0to5_6to13[i,9]<-(length(unlist(fgsea_pax_Neg_0to5_6to13[i,8])))/fgsea_pax_Neg_0to5_6to13[i,7]}
fgsea_pax_Neg_0to5_6to13 <- data.frame(fgsea_pax_Neg_0to5_6to13 %>% 
                                         mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_Neg_0to5_6to13$leadingEdge <- NULL
fgsea_pax_Neg_0to5_6to13$Sample <- "0-5 vs. 6-13 yr"
write_xlsx(fgsea_pax_Neg_0to5_6to13, "Statistical_Analyses/1_COVID_Neg_by_Age/modules_pax_Neg_0to5_6to13.xlsx")

res_pax_Neg_0to5_14to20 <- as_tibble(dds_pax_Neg_0to5_14to20)
res2_pax_Neg_0to5_14to20 <- res_pax_Neg_0to5_14to20 %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_Neg_0to5_14to20 <- deframe(res2_pax_Neg_0to5_14to20)
fgsea_pax_Neg_0to5_14to20 <- fgsea(pathways=modules_61, stats=ranks_pax_Neg_0to5_14to20)
fgsea_pax_Neg_0to5_14to20$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_Neg_0to5_14to20))) {fgsea_pax_Neg_0to5_14to20[i,9]<-(length(unlist(fgsea_pax_Neg_0to5_14to20[i,8])))/fgsea_pax_Neg_0to5_14to20[i,7]}
fgsea_pax_Neg_0to5_14to20 <- data.frame(fgsea_pax_Neg_0to5_14to20 %>% 
                                          mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_Neg_0to5_14to20$leadingEdge <- NULL
fgsea_pax_Neg_0to5_14to20$Sample <- "0-5 vs. 14-20 yr"
write_xlsx(fgsea_pax_Neg_0to5_14to20, "Statistical_Analyses/1_COVID_Neg_by_Age/modules_pax_Neg_0to5_14to20.xlsx")

res_pax_Neg_6to13_14to20 <- as_tibble(dds_pax_Neg_6to13_14to20)
res2_pax_Neg_6to13_14to20 <- res_pax_Neg_6to13_14to20 %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_Neg_6to13_14to20 <- deframe(res2_pax_Neg_6to13_14to20)
fgsea_pax_Neg_6to13_14to20 <- fgsea(pathways=modules_61, stats=ranks_pax_Neg_6to13_14to20)
fgsea_pax_Neg_6to13_14to20$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_Neg_6to13_14to20))) {fgsea_pax_Neg_6to13_14to20[i,9]<-(length(unlist(fgsea_pax_Neg_6to13_14to20[i,8])))/fgsea_pax_Neg_6to13_14to20[i,7]}
fgsea_pax_Neg_6to13_14to20 <- data.frame(fgsea_pax_Neg_6to13_14to20 %>% 
                                           mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_Neg_6to13_14to20$leadingEdge <- NULL
fgsea_pax_Neg_6to13_14to20$Sample <- "6-13 vs. 14-20 yr"
write_xlsx(fgsea_pax_Neg_6to13_14to20, "Statistical_Analyses/1_COVID_Neg_by_Age/modules_pax_Neg_6to13_14to20.xlsx")

res_pax_Neg_0to5_Adult <- as_tibble(dds_pax_Neg_0to5_Adult)
res2_pax_Neg_0to5_Adult <- res_pax_Neg_0to5_Adult %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_Neg_0to5_Adult <- deframe(res2_pax_Neg_0to5_Adult)
fgsea_pax_Neg_0to5_Adult <- fgsea(pathways=modules_61, stats=ranks_pax_Neg_0to5_Adult)
fgsea_pax_Neg_0to5_Adult$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_Neg_0to5_Adult))) {fgsea_pax_Neg_0to5_Adult[i,9]<-(length(unlist(fgsea_pax_Neg_0to5_Adult[i,8])))/fgsea_pax_Neg_0to5_Adult[i,7]}
fgsea_pax_Neg_0to5_Adult <- data.frame(fgsea_pax_Neg_0to5_Adult %>% 
                                          mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_Neg_0to5_Adult$leadingEdge <- NULL
fgsea_pax_Neg_0to5_Adult$Sample <- "0-5 yr vs. Adult"
write_xlsx(fgsea_pax_Neg_0to5_Adult, "Statistical_Analyses/1_COVID_Neg_by_Age/modules_pax_Neg_0to5_Adult.xlsx")

res_pax_Neg_6to13_Adult <- as_tibble(dds_pax_Neg_6to13_Adult)
res2_pax_Neg_6to13_Adult <- res_pax_Neg_6to13_Adult %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_Neg_6to13_Adult <- deframe(res2_pax_Neg_6to13_Adult)
fgsea_pax_Neg_6to13_Adult <- fgsea(pathways=modules_61, stats=ranks_pax_Neg_6to13_Adult)
fgsea_pax_Neg_6to13_Adult$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_Neg_6to13_Adult))) {fgsea_pax_Neg_6to13_Adult[i,9]<-(length(unlist(fgsea_pax_Neg_6to13_Adult[i,8])))/fgsea_pax_Neg_6to13_Adult[i,7]}
fgsea_pax_Neg_6to13_Adult <- data.frame(fgsea_pax_Neg_6to13_Adult %>% 
                                           mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_Neg_6to13_Adult$leadingEdge <- NULL
fgsea_pax_Neg_6to13_Adult$Sample <- "6-13 yr vs. Adult"
write_xlsx(fgsea_pax_Neg_6to13_Adult, "Statistical_Analyses/1_COVID_Neg_by_Age/modules_pax_Neg_6to13_Adult.xlsx")

res_pax_Neg_14to20_Adult <- as_tibble(dds_pax_Neg_14to20_Adult)
res2_pax_Neg_14to20_Adult <- res_pax_Neg_14to20_Adult %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_Neg_14to20_Adult <- deframe(res2_pax_Neg_14to20_Adult)
fgsea_pax_Neg_14to20_Adult <- fgsea(pathways=modules_61, stats=ranks_pax_Neg_14to20_Adult)
fgsea_pax_Neg_14to20_Adult$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_Neg_14to20_Adult))) {fgsea_pax_Neg_14to20_Adult[i,9]<-(length(unlist(fgsea_pax_Neg_14to20_Adult[i,8])))/fgsea_pax_Neg_14to20_Adult[i,7]}
fgsea_pax_Neg_14to20_Adult <- data.frame(fgsea_pax_Neg_14to20_Adult %>% 
                                          mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_Neg_14to20_Adult$leadingEdge <- NULL
fgsea_pax_Neg_14to20_Adult$Sample <- "14-20 yr vs. Adult"
write_xlsx(fgsea_pax_Neg_14to20_Adult, "Statistical_Analyses/1_COVID_Neg_by_Age/modules_pax_Neg_14to20_Adult.xlsx")

res_pax_Neg_Child_Adult <- as_tibble(dds_pax_Neg_Child_Adult)
res2_pax_Neg_Child_Adult <- res_pax_Neg_Child_Adult %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_Neg_Child_Adult <- deframe(res2_pax_Neg_Child_Adult)
fgsea_pax_Neg_Child_Adult <- fgsea(pathways=modules_61, stats=ranks_pax_Neg_Child_Adult)
fgsea_pax_Neg_Child_Adult$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_Neg_Child_Adult))) {fgsea_pax_Neg_Child_Adult[i,9]<-(length(unlist(fgsea_pax_Neg_Child_Adult[i,8])))/fgsea_pax_Neg_Child_Adult[i,7]}
fgsea_pax_Neg_Child_Adult <- data.frame(fgsea_pax_Neg_Child_Adult %>% 
                                         mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_Neg_Child_Adult$leadingEdge <- NULL
fgsea_pax_Neg_Child_Adult$Sample <- "Child vs. Adult"
write_xlsx(fgsea_pax_Neg_Child_Adult, "Statistical_Analyses/1_COVID_Neg_by_Age/modules_pax_Neg_Child_Adult.xlsx")

####################################################
### DIFFERENCES IN NP EXPRESSION BY COVID STATUS ###
####################################################

res_np_Pos_Neg <- as_tibble(dds_np_Pos_Neg)
res2_np_Pos_Neg <- res_np_Pos_Neg %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_Pos_Neg <- deframe(res2_np_Pos_Neg)
fgsea_np_Pos_Neg <- fgsea(pathways=modules_61, stats=ranks_np_Pos_Neg)
fgsea_np_Pos_Neg$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_Pos_Neg))) {fgsea_np_Pos_Neg[i,9]<-(length(unlist(fgsea_np_Pos_Neg[i,8])))/fgsea_np_Pos_Neg[i,7]}
fgsea_np_Pos_Neg <- data.frame(fgsea_np_Pos_Neg %>% 
                                         mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_Pos_Neg$leadingEdge <- NULL
write_xlsx(fgsea_np_Pos_Neg, "Statistical_Analyses/2_COVID_Pos_vs_Neg/modules_np_Pos_Neg.xlsx")

##########################################################
### DIFFERENCES IN PAX GENE EXPRESSION BY COVID STATUS ###
##########################################################

res_pax_Pos_Neg <- as_tibble(dds_pax_Pos_Neg)
res2_pax_Pos_Neg <- res_pax_Pos_Neg %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_Pos_Neg <- deframe(res2_pax_Pos_Neg)
fgsea_pax_Pos_Neg <- fgsea(pathways=modules_61, stats=ranks_pax_Pos_Neg)
fgsea_pax_Pos_Neg$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_Pos_Neg))) {fgsea_pax_Pos_Neg[i,9]<-(length(unlist(fgsea_pax_Pos_Neg[i,8])))/fgsea_pax_Pos_Neg[i,7]}
fgsea_pax_Pos_Neg <- data.frame(fgsea_pax_Pos_Neg %>% 
                                 mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_Pos_Neg$leadingEdge <- NULL
write_xlsx(fgsea_pax_Pos_Neg, "Statistical_Analyses/2_COVID_Pos_vs_Neg/modules_pax_Pos_Neg.xlsx")

##############################################################
### COMPARING NP SAMPLES BY COVID STATUS IN AGE CATEGORIES ###
##############################################################

res_np_0to5_Pos_Neg <- as_tibble(dds_np_0to5_Pos_Neg)
res2_np_0to5_Pos_Neg <- res_np_0to5_Pos_Neg %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_0to5_Pos_Neg <- deframe(res2_np_0to5_Pos_Neg)
fgsea_np_0to5_Pos_Neg <- fgsea(pathways=modules_61, stats=ranks_np_0to5_Pos_Neg)
fgsea_np_0to5_Pos_Neg$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_0to5_Pos_Neg))) {fgsea_np_0to5_Pos_Neg[i,9]<-(length(unlist(fgsea_np_0to5_Pos_Neg[i,8])))/fgsea_np_0to5_Pos_Neg[i,7]}
fgsea_np_0to5_Pos_Neg <- data.frame(fgsea_np_0to5_Pos_Neg %>% 
                                        mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_0to5_Pos_Neg$leadingEdge <- NULL
fgsea_np_0to5_Pos_Neg$Sample <- "0-5 yr"
write_xlsx(fgsea_np_0to5_Pos_Neg, "Statistical_Analyses/3_COVID_Pos_by_Age/modules_np_0to5_Pos_Neg.xlsx")

res_np_6to13_Pos_Neg <- as_tibble(dds_np_6to13_Pos_Neg)
res2_np_6to13_Pos_Neg <- res_np_6to13_Pos_Neg %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_6to13_Pos_Neg <- deframe(res2_np_6to13_Pos_Neg)
fgsea_np_6to13_Pos_Neg <- fgsea(pathways=modules_61, stats=ranks_np_6to13_Pos_Neg)
fgsea_np_6to13_Pos_Neg$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_6to13_Pos_Neg))) {fgsea_np_6to13_Pos_Neg[i,9]<-(length(unlist(fgsea_np_6to13_Pos_Neg[i,8])))/fgsea_np_6to13_Pos_Neg[i,7]}
fgsea_np_6to13_Pos_Neg <- data.frame(fgsea_np_6to13_Pos_Neg %>% 
                                      mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_6to13_Pos_Neg$leadingEdge <- NULL
fgsea_np_6to13_Pos_Neg$Sample <- "6-13 yr"
write_xlsx(fgsea_np_6to13_Pos_Neg, "Statistical_Analyses/3_COVID_Pos_by_Age/modules_np_6to13_Pos_Neg.xlsx")

res_np_14to20_Pos_Neg <- as_tibble(dds_np_14to20_Pos_Neg)
res2_np_14to20_Pos_Neg <- res_np_14to20_Pos_Neg %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_14to20_Pos_Neg <- deframe(res2_np_14to20_Pos_Neg)
fgsea_np_14to20_Pos_Neg <- fgsea(pathways=modules_61, stats=ranks_np_14to20_Pos_Neg)
fgsea_np_14to20_Pos_Neg$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_14to20_Pos_Neg))) {fgsea_np_14to20_Pos_Neg[i,9]<-(length(unlist(fgsea_np_14to20_Pos_Neg[i,8])))/fgsea_np_14to20_Pos_Neg[i,7]}
fgsea_np_14to20_Pos_Neg <- data.frame(fgsea_np_14to20_Pos_Neg %>% 
                                       mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_14to20_Pos_Neg$leadingEdge <- NULL
fgsea_np_14to20_Pos_Neg$Sample <- "14-20 yr"
write_xlsx(fgsea_np_14to20_Pos_Neg, "Statistical_Analyses/3_COVID_Pos_by_Age/modules_np_14to20_Pos_Neg.xlsx")

###############################################################
### COMPARING PAX SAMPLES BY COVID STATUS IN AGE CATEGORIES ###
###############################################################

res_pax_0to5_Pos_Neg <- as_tibble(dds_pax_0to5_Pos_Neg)
res2_pax_0to5_Pos_Neg <- res_pax_0to5_Pos_Neg %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_0to5_Pos_Neg <- deframe(res2_pax_0to5_Pos_Neg)
fgsea_pax_0to5_Pos_Neg <- fgsea(pathways=modules_61, stats=ranks_pax_0to5_Pos_Neg)
fgsea_pax_0to5_Pos_Neg$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_0to5_Pos_Neg))) {fgsea_pax_0to5_Pos_Neg[i,9]<-(length(unlist(fgsea_pax_0to5_Pos_Neg[i,8])))/fgsea_pax_0to5_Pos_Neg[i,7]}
fgsea_pax_0to5_Pos_Neg <- data.frame(fgsea_pax_0to5_Pos_Neg %>% 
                                      mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_0to5_Pos_Neg$leadingEdge <- NULL
fgsea_pax_0to5_Pos_Neg$Sample <- "0-5 yr"
write_xlsx(fgsea_pax_0to5_Pos_Neg, "Statistical_Analyses/3_COVID_Pos_by_Age/modules_pax_0to5_Pos_Neg.xlsx")

res_pax_6to13_Pos_Neg <- as_tibble(dds_pax_6to13_Pos_Neg)
res2_pax_6to13_Pos_Neg <- res_pax_6to13_Pos_Neg %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_6to13_Pos_Neg <- deframe(res2_pax_6to13_Pos_Neg)
fgsea_pax_6to13_Pos_Neg <- fgsea(pathways=modules_61, stats=ranks_pax_6to13_Pos_Neg)
fgsea_pax_6to13_Pos_Neg$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_6to13_Pos_Neg))) {fgsea_pax_6to13_Pos_Neg[i,9]<-(length(unlist(fgsea_pax_6to13_Pos_Neg[i,8])))/fgsea_pax_6to13_Pos_Neg[i,7]}
fgsea_pax_6to13_Pos_Neg <- data.frame(fgsea_pax_6to13_Pos_Neg %>% 
                                       mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_6to13_Pos_Neg$leadingEdge <- NULL
fgsea_pax_6to13_Pos_Neg$Sample <- "6-13 yr"
write_xlsx(fgsea_pax_6to13_Pos_Neg, "Statistical_Analyses/3_COVID_Pos_by_Age/modules_pax_6to13_Pos_Neg.xlsx")

res_pax_14to20_Pos_Neg <- as_tibble(dds_pax_14to20_Pos_Neg)
res2_pax_14to20_Pos_Neg <- res_pax_14to20_Pos_Neg %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_14to20_Pos_Neg <- deframe(res2_pax_14to20_Pos_Neg)
fgsea_pax_14to20_Pos_Neg <- fgsea(pathways=modules_61, stats=ranks_pax_14to20_Pos_Neg)
fgsea_pax_14to20_Pos_Neg$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_14to20_Pos_Neg))) {fgsea_pax_14to20_Pos_Neg[i,9]<-(length(unlist(fgsea_pax_14to20_Pos_Neg[i,8])))/fgsea_pax_14to20_Pos_Neg[i,7]}
fgsea_pax_14to20_Pos_Neg <- data.frame(fgsea_pax_14to20_Pos_Neg %>% 
                                        mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_14to20_Pos_Neg$leadingEdge <- NULL
fgsea_pax_14to20_Pos_Neg$Sample <- "14-20 yr"
write_xlsx(fgsea_pax_14to20_Pos_Neg, "Statistical_Analyses/3_COVID_Pos_by_Age/modules_pax_14to20_Pos_Neg.xlsx")

res_pax_adult_Pos_Neg <- as_tibble(dds_pax_adult_Pos_Neg)
res2_pax_adult_Pos_Neg <- res_pax_adult_Pos_Neg %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_adult_Pos_Neg <- deframe(res2_pax_adult_Pos_Neg)
fgsea_pax_adult_Pos_Neg <- fgsea(pathways=modules_61, stats=ranks_pax_adult_Pos_Neg)
fgsea_pax_adult_Pos_Neg$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_adult_Pos_Neg))) {fgsea_pax_adult_Pos_Neg[i,9]<-(length(unlist(fgsea_pax_adult_Pos_Neg[i,8])))/fgsea_pax_adult_Pos_Neg[i,7]}
fgsea_pax_adult_Pos_Neg <- data.frame(fgsea_pax_adult_Pos_Neg %>% 
                                         mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_adult_Pos_Neg$leadingEdge <- NULL
fgsea_pax_adult_Pos_Neg$Sample <- "Adult"
write_xlsx(fgsea_pax_adult_Pos_Neg, "Statistical_Analyses/3_COVID_Pos_by_Age/modules_pax_adult_Pos_Neg.xlsx")

###############################################################################
### MODELS FOR PATIENT CHARACTERISTICS & SYMPTOMS AMONG COVID+ - NP SAMPLES ###
###############################################################################

# Fever
res_np_fever <- as_tibble(dds_np_fever)
res2_np_fever <- res_np_fever %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_fever <- deframe(res2_np_fever)
fgsea_np_fever <- fgsea(pathways=modules_61, stats=ranks_np_fever)
fgsea_np_fever$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_fever))) {fgsea_np_fever[i,9]<-(length(unlist(fgsea_np_fever[i,8])))/fgsea_np_fever[i,7]}
fgsea_np_fever <- data.frame(fgsea_np_fever %>% 
                                mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_fever$leadingEdge <- NULL
write_xlsx(fgsea_np_fever, "Statistical_Analyses/4_Symptoms/modules_np_fever.xlsx")

# Cough
res_np_cough <- as_tibble(dds_np_cough)
res2_np_cough <- res_np_cough %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_cough <- deframe(res2_np_cough)
fgsea_np_cough <- fgsea(pathways=modules_61, stats=ranks_np_cough)
fgsea_np_cough$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_cough))) {fgsea_np_cough[i,9]<-(length(unlist(fgsea_np_cough[i,8])))/fgsea_np_cough[i,7]}
fgsea_np_cough <- data.frame(fgsea_np_cough %>% 
                                mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_cough$leadingEdge <- NULL
write_xlsx(fgsea_np_cough, "Statistical_Analyses/4_Symptoms/modules_np_cough.xlsx")

# Rhinorrhea
res_np_rhinorrhea <- as_tibble(dds_np_rhinorrhea)
res2_np_rhinorrhea <- res_np_rhinorrhea %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_rhinorrhea <- deframe(res2_np_rhinorrhea)
fgsea_np_rhinorrhea <- fgsea(pathways=modules_61, stats=ranks_np_rhinorrhea)
fgsea_np_rhinorrhea$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_rhinorrhea))) {fgsea_np_rhinorrhea[i,9]<-(length(unlist(fgsea_np_rhinorrhea[i,8])))/fgsea_np_rhinorrhea[i,7]}
fgsea_np_rhinorrhea <- data.frame(fgsea_np_rhinorrhea %>% 
                                     mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_rhinorrhea$leadingEdge <- NULL
write_xlsx(fgsea_np_rhinorrhea, "Statistical_Analyses/4_Symptoms/modules_np_rhinorrhea.xlsx")

# Nasal congestion
res_np_congestion <- as_tibble(dds_np_congestion)
res2_np_congestion <- res_np_congestion %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_congestion <- deframe(res2_np_congestion)
fgsea_np_congestion <- fgsea(pathways=modules_61, stats=ranks_np_congestion)
fgsea_np_congestion$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_congestion))) {fgsea_np_congestion[i,9]<-(length(unlist(fgsea_np_congestion[i,8])))/fgsea_np_congestion[i,7]}
fgsea_np_congestion <- data.frame(fgsea_np_congestion %>% 
                                     mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_congestion$leadingEdge <- NULL
write_xlsx(fgsea_np_congestion, "Statistical_Analyses/4_Symptoms/modules_np_congestion.xlsx")

# Headache
res_np_headache <- as_tibble(dds_np_headache)
res2_np_headache <- res_np_headache %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_headache <- deframe(res2_np_headache)
fgsea_np_headache <- fgsea(pathways=modules_61, stats=ranks_np_headache)
fgsea_np_headache$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_headache))) {fgsea_np_headache[i,9]<-(length(unlist(fgsea_np_headache[i,8])))/fgsea_np_headache[i,7]}
fgsea_np_headache <- data.frame(fgsea_np_headache %>% 
                               mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_headache$leadingEdge <- NULL
write_xlsx(fgsea_np_headache, "Statistical_Analyses/4_Symptoms/modules_np_headache.xlsx")

# Abdominal pain
res_np_abd_pain <- as_tibble(dds_np_abd_pain)
res2_np_abd_pain <- res_np_abd_pain %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_abd_pain <- deframe(res2_np_abd_pain)
fgsea_np_abd_pain <- fgsea(pathways=modules_61, stats=ranks_np_abd_pain)
fgsea_np_abd_pain$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_abd_pain))) {fgsea_np_abd_pain[i,9]<-(length(unlist(fgsea_np_abd_pain[i,8])))/fgsea_np_abd_pain[i,7]}
fgsea_np_abd_pain <- data.frame(fgsea_np_abd_pain %>% 
                               mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_abd_pain$leadingEdge <- NULL
write_xlsx(fgsea_np_abd_pain, "Statistical_Analyses/4_Symptoms/modules_np_abd_pain.xlsx")

# Anosmia
res_np_anosmia <- as_tibble(dds_np_anosmia)
res2_np_anosmia <- res_np_anosmia %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_anosmia <- deframe(res2_np_anosmia)
fgsea_np_anosmia <- fgsea(pathways=modules_61, stats=ranks_np_anosmia)
fgsea_np_anosmia$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_anosmia))) {fgsea_np_anosmia[i,9]<-(length(unlist(fgsea_np_anosmia[i,8])))/fgsea_np_anosmia[i,7]}
fgsea_np_anosmia <- data.frame(fgsea_np_anosmia %>% 
                               mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_anosmia$leadingEdge <- NULL
write_xlsx(fgsea_np_anosmia, "Statistical_Analyses/4_Symptoms/modules_np_anosmia.xlsx")

# Dysgeusia
res_np_dysgeusia <- as_tibble(dds_np_dysgeusia)
res2_np_dysgeusia <- res_np_dysgeusia %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_dysgeusia <- deframe(res2_np_dysgeusia)
fgsea_np_dysgeusia <- fgsea(pathways=modules_61, stats=ranks_np_dysgeusia)
fgsea_np_dysgeusia$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_dysgeusia))) {fgsea_np_dysgeusia[i,9]<-(length(unlist(fgsea_np_dysgeusia[i,8])))/fgsea_np_dysgeusia[i,7]}
fgsea_np_dysgeusia <- data.frame(fgsea_np_dysgeusia %>% 
                               mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_dysgeusia$leadingEdge <- NULL
write_xlsx(fgsea_np_dysgeusia, "Statistical_Analyses/4_Symptoms/modules_np_dysgeusia.xlsx")

# Myalgias
res_np_myalgias <- as_tibble(dds_np_myalgias)
res2_np_myalgias <- res_np_myalgias %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_myalgias <- deframe(res2_np_myalgias)
fgsea_np_myalgias <- fgsea(pathways=modules_61, stats=ranks_np_myalgias)
fgsea_np_myalgias$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_myalgias))) {fgsea_np_myalgias[i,9]<-(length(unlist(fgsea_np_myalgias[i,8])))/fgsea_np_myalgias[i,7]}
fgsea_np_myalgias <- data.frame(fgsea_np_myalgias %>% 
                               mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_myalgias$leadingEdge <- NULL
write_xlsx(fgsea_np_myalgias, "Statistical_Analyses/4_Symptoms/modules_np_myalgias.xlsx")

################################################################################
### MODELS FOR PATIENT CHARACTERISTICS & SYMPTOMS AMONG COVID+ - PAX SAMPLES ###
################################################################################

# Fever
res_pax_fever <- as_tibble(dds_pax_fever)
res2_pax_fever <- res_pax_fever %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_fever <- deframe(res2_pax_fever)
fgsea_pax_fever <- fgsea(pathways=modules_61, stats=ranks_pax_fever)
fgsea_pax_fever$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_fever))) {fgsea_pax_fever[i,9]<-(length(unlist(fgsea_pax_fever[i,8])))/fgsea_pax_fever[i,7]}
fgsea_pax_fever <- data.frame(fgsea_pax_fever %>% 
                               mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_fever$leadingEdge <- NULL
write_xlsx(fgsea_pax_fever, "Statistical_Analyses/4_Symptoms/modules_pax_fever.xlsx")

# Cough
res_pax_cough <- as_tibble(dds_pax_cough)
res2_pax_cough <- res_pax_cough %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_cough <- deframe(res2_pax_cough)
fgsea_pax_cough <- fgsea(pathways=modules_61, stats=ranks_pax_cough)
fgsea_pax_cough$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_cough))) {fgsea_pax_cough[i,9]<-(length(unlist(fgsea_pax_cough[i,8])))/fgsea_pax_cough[i,7]}
fgsea_pax_cough <- data.frame(fgsea_pax_cough %>% 
                               mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_cough$leadingEdge <- NULL
write_xlsx(fgsea_pax_cough, "Statistical_Analyses/4_Symptoms/modules_pax_cough.xlsx")

# Rhinorrhea
res_pax_rhinorrhea <- as_tibble(dds_pax_rhinorrhea)
res2_pax_rhinorrhea <- res_pax_rhinorrhea %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_rhinorrhea <- deframe(res2_pax_rhinorrhea)
fgsea_pax_rhinorrhea <- fgsea(pathways=modules_61, stats=ranks_pax_rhinorrhea)
fgsea_pax_rhinorrhea$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_rhinorrhea))) {fgsea_pax_rhinorrhea[i,9]<-(length(unlist(fgsea_pax_rhinorrhea[i,8])))/fgsea_pax_rhinorrhea[i,7]}
fgsea_pax_rhinorrhea <- data.frame(fgsea_pax_rhinorrhea %>% 
                                    mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_rhinorrhea$leadingEdge <- NULL
write_xlsx(fgsea_pax_rhinorrhea, "Statistical_Analyses/4_Symptoms/modules_pax_rhinorrhea.xlsx")

# Nasal congestion
res_pax_congestion <- as_tibble(dds_pax_congestion)
res2_pax_congestion <- res_pax_congestion %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_congestion <- deframe(res2_pax_congestion)
fgsea_pax_congestion <- fgsea(pathways=modules_61, stats=ranks_pax_congestion)
fgsea_pax_congestion$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_congestion))) {fgsea_pax_congestion[i,9]<-(length(unlist(fgsea_pax_congestion[i,8])))/fgsea_pax_congestion[i,7]}
fgsea_pax_congestion <- data.frame(fgsea_pax_congestion %>% 
                                    mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_congestion$leadingEdge <- NULL
write_xlsx(fgsea_pax_congestion, "Statistical_Analyses/4_Symptoms/modules_pax_congestion.xlsx")

# Headache
res_pax_headache <- as_tibble(dds_pax_headache)
res2_pax_headache <- res_pax_headache %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_headache <- deframe(res2_pax_headache)
fgsea_pax_headache <- fgsea(pathways=modules_61, stats=ranks_pax_headache)
fgsea_pax_headache$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_headache))) {fgsea_pax_headache[i,9]<-(length(unlist(fgsea_pax_headache[i,8])))/fgsea_pax_headache[i,7]}
fgsea_pax_headache <- data.frame(fgsea_pax_headache %>% 
                                  mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_headache$leadingEdge <- NULL
write_xlsx(fgsea_pax_headache, "Statistical_Analyses/4_Symptoms/modules_pax_headache.xlsx")

# Abdominal pain
res_pax_abd_pain <- as_tibble(dds_pax_abd_pain)
res2_pax_abd_pain <- res_pax_abd_pain %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_abd_pain <- deframe(res2_pax_abd_pain)
fgsea_pax_abd_pain <- fgsea(pathways=modules_61, stats=ranks_pax_abd_pain)
fgsea_pax_abd_pain$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_abd_pain))) {fgsea_pax_abd_pain[i,9]<-(length(unlist(fgsea_pax_abd_pain[i,8])))/fgsea_pax_abd_pain[i,7]}
fgsea_pax_abd_pain <- data.frame(fgsea_pax_abd_pain %>% 
                                  mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_abd_pain$leadingEdge <- NULL
write_xlsx(fgsea_pax_abd_pain, "Statistical_Analyses/4_Symptoms/modules_pax_abd_pain.xlsx")

# Anosmia
res_pax_anosmia <- as_tibble(dds_pax_anosmia)
res2_pax_anosmia <- res_pax_anosmia %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_anosmia <- deframe(res2_pax_anosmia)
fgsea_pax_anosmia <- fgsea(pathways=modules_61, stats=ranks_pax_anosmia)
fgsea_pax_anosmia$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_anosmia))) {fgsea_pax_anosmia[i,9]<-(length(unlist(fgsea_pax_anosmia[i,8])))/fgsea_pax_anosmia[i,7]}
fgsea_pax_anosmia <- data.frame(fgsea_pax_anosmia %>% 
                                 mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_anosmia$leadingEdge <- NULL
write_xlsx(fgsea_pax_anosmia, "Statistical_Analyses/4_Symptoms/modules_pax_anosmia.xlsx")

# Dysgeusia
res_pax_dysgeusia <- as_tibble(dds_pax_dysgeusia)
res2_pax_dysgeusia <- res_pax_dysgeusia %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_dysgeusia <- deframe(res2_pax_dysgeusia)
fgsea_pax_dysgeusia <- fgsea(pathways=modules_61, stats=ranks_pax_dysgeusia)
fgsea_pax_dysgeusia$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_dysgeusia))) {fgsea_pax_dysgeusia[i,9]<-(length(unlist(fgsea_pax_dysgeusia[i,8])))/fgsea_pax_dysgeusia[i,7]}
fgsea_pax_dysgeusia <- data.frame(fgsea_pax_dysgeusia %>% 
                                   mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_dysgeusia$leadingEdge <- NULL
write_xlsx(fgsea_pax_dysgeusia, "Statistical_Analyses/4_Symptoms/modules_pax_dysgeusia.xlsx")

# Myalgias
res_pax_myalgias <- as_tibble(dds_pax_myalgias)
res2_pax_myalgias <- res_pax_myalgias %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_myalgias <- deframe(res2_pax_myalgias)
fgsea_pax_myalgias <- fgsea(pathways=modules_61, stats=ranks_pax_myalgias)
fgsea_pax_myalgias$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_myalgias))) {fgsea_pax_myalgias[i,9]<-(length(unlist(fgsea_pax_myalgias[i,8])))/fgsea_pax_myalgias[i,7]}
fgsea_pax_myalgias <- data.frame(fgsea_pax_myalgias %>% 
                                  mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_myalgias$leadingEdge <- NULL
write_xlsx(fgsea_pax_myalgias, "Statistical_Analyses/4_Symptoms/modules_pax_myalgias.xlsx")