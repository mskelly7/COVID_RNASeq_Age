# BRAVE Kids RNA Sequencing Analysis
# Aditya Mohan (MD/PhD candidate)  / Matthew Kelly, MD, MPH 
# Gene Set Enrichment Analysis using FGSEA
# Last update: July 23, 2023

remove(list=ls())
setwd("___________________") 
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
modules_61 <- gmtPathways("modules_61.gmt")

# Upload DESeq2 output files
dds_np_Neg_0to5_6to13 <- read.csv("1_COVID_Neg_by_Age/genes_np_Neg_0to5_6to13.csv")
dds_np_Neg_0to5_14to20 <- read.csv("1_COVID_Neg_by_Age/genes_np_Neg_0to5_14to20.csv")
dds_np_Neg_6to13_14to20 <- read.csv("1_COVID_Neg_by_Age/genes_np_Neg_6to13_14to20.csv")
dds_pax_Neg_0to5_6to13 <- read.csv("1_COVID_Neg_by_Age/genes_pax_Neg_0to5_6to13.csv")
dds_pax_Neg_0to5_14to20 <- read.csv("1_COVID_Neg_by_Age/genes_pax_Neg_0to5_14to20.csv")
dds_pax_Neg_6to13_14to20 <- read.csv("1_COVID_Neg_by_Age/genes_pax_Neg_6to13_14to20.csv")
dds_pax_Neg_0to5_Adult <- read.csv("1_COVID_Neg_by_Age/genes_pax_Neg_0to5_Adult.csv")
dds_pax_Neg_6to13_Adult <- read.csv("1_COVID_Neg_by_Age/genes_pax_Neg_6to13_Adult.csv")
dds_pax_Neg_14to20_Adult <- read.csv("1_COVID_Neg_by_Age/genes_pax_Neg_14to20_Adult.csv")
dds_pax_Neg_Child_Adult <- read.csv("1_COVID_Neg_by_Age/genes_pax_Neg_Child_Adult.csv")
dds_np_Pos_Neg <- read.csv("2_COVID_Pos_vs_Neg/genes_np_Pos_Neg.csv")
dds_pax_Pos_Neg <- read.csv("2_COVID_Pos_vs_Neg/genes_pax_Pos_Neg.csv")
dds_np_0to5_Pos_Neg <- read.csv("3_COVID_Pos_by_Age/genes_np_0to5_Pos_Neg.csv")
dds_np_6to13_Pos_Neg <- read.csv("3_COVID_Pos_by_Age/genes_np_6to13_Pos_Neg.csv")
dds_np_14to20_Pos_Neg <- read.csv("3_COVID_Pos_by_Age/genes_np_14to20_Pos_Neg.csv")
dds_pax_0to5_Pos_Neg <- read.csv("3_COVID_Pos_by_Age/genes_pax_0to5_Pos_Neg.csv")
dds_pax_6to13_Pos_Neg <- read.csv("3_COVID_Pos_by_Age/genes_pax_6to13_Pos_Neg.csv")
dds_pax_14to20_Pos_Neg <- read.csv("3_COVID_Pos_by_Age/genes_pax_14to20_Pos_Neg.csv")
dds_pax_adult_Pos_Neg <- read.csv("3_COVID_Pos_by_Age/genes_pax_adult_Pos_Neg.csv")
dds_np_obesity <- read.csv("4_Symptoms/genes_np_obesity.csv")
dds_np_asthma <- read.csv("4_Symptoms/genes_np_asthma.csv")
dds_np_fever <- read.csv("4_Symptoms/genes_np_fever.csv")
dds_np_cough <- read.csv("4_Symptoms/genes_np_cough.csv")
dds_np_rhinorrhea <- read.csv("4_Symptoms/genes_np_rhinorrhea.csv")
dds_np_congestion <- read.csv("4_Symptoms/genes_np_congestion.csv")
dds_np_vl_copies <- read.csv("4_Symptoms/genes_np_vl_copies.csv")
dds_pax_obesity <- read.csv("4_Symptoms/genes_pax_obesity.csv")
dds_pax_asthma <- read.csv("4_Symptoms/genes_pax_asthma.csv")
dds_pax_fever <- read.csv("4_Symptoms/genes_pax_fever.csv")
dds_pax_cough <- read.csv("4_Symptoms/genes_pax_cough.csv")
dds_pax_rhinorrhea <- read.csv("4_Symptoms/genes_pax_rhinorrhea.csv")
dds_pax_congestion <- read.csv("4_Symptoms/genes_pax_congestion.csv")
dds_pax_vl_copies <- read.csv("4_Symptoms/genes_pax_vl_copies.csv")

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
write_xlsx(fgsea_np_Neg_0to5_6to13, "1_COVID_Neg_by_Age/modules_np_Neg_0to5_6to13.xlsx")

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
write_xlsx(fgsea_np_Neg_0to5_14to20, "1_COVID_Neg_by_Age/modules_np_Neg_0to5_14to20.xlsx")

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
write_xlsx(fgsea_np_Neg_6to13_14to20, "1_COVID_Neg_by_Age/modules_np_Neg_6to13_14to20.xlsx")

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
write_xlsx(fgsea_pax_Neg_0to5_6to13, "1_COVID_Neg_by_Age/modules_pax_Neg_0to5_6to13.xlsx")

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
write_xlsx(fgsea_pax_Neg_0to5_14to20, "1_COVID_Neg_by_Age/modules_pax_Neg_0to5_14to20.xlsx")

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
write_xlsx(fgsea_pax_Neg_6to13_14to20, "1_COVID_Neg_by_Age/modules_pax_Neg_6to13_14to20.xlsx")

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
write_xlsx(fgsea_pax_Neg_0to5_Adult, "1_COVID_Neg_by_Age/modules_pax_Neg_0to5_Adult.xlsx")

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
write_xlsx(fgsea_pax_Neg_6to13_Adult, "1_COVID_Neg_by_Age/modules_pax_Neg_6to13_Adult.xlsx")

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
write_xlsx(fgsea_pax_Neg_14to20_Adult, "1_COVID_Neg_by_Age/modules_pax_Neg_14to20_Adult.xlsx")

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
write_xlsx(fgsea_pax_Neg_Child_Adult, "1_COVID_Neg_by_Age/modules_pax_Neg_Child_Adult.xlsx")

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
write_xlsx(fgsea_np_Pos_Neg, "2_COVID_Pos_vs_Neg/modules_np_Pos_Neg.xlsx")

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
write_xlsx(fgsea_pax_Pos_Neg, "2_COVID_Pos_vs_Neg/modules_pax_Pos_Neg.xlsx")

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
write_xlsx(fgsea_np_0to5_Pos_Neg, "3_COVID_Pos_by_Age/modules_np_0to5_Pos_Neg.xlsx")

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
write_xlsx(fgsea_np_6to13_Pos_Neg, "3_COVID_Pos_by_Age/modules_np_6to13_Pos_Neg.xlsx")

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
write_xlsx(fgsea_np_14to20_Pos_Neg, "3_COVID_Pos_by_Age/modules_np_14to20_Pos_Neg.xlsx")

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
write_xlsx(fgsea_pax_0to5_Pos_Neg, "3_COVID_Pos_by_Age/modules_pax_0to5_Pos_Neg.xlsx")

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
write_xlsx(fgsea_pax_6to13_Pos_Neg, "3_COVID_Pos_by_Age/modules_pax_6to13_Pos_Neg.xlsx")

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
write_xlsx(fgsea_pax_14to20_Pos_Neg, "3_COVID_Pos_by_Age/modules_pax_14to20_Pos_Neg.xlsx")

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
write_xlsx(fgsea_pax_adult_Pos_Neg, "3_COVID_Pos_by_Age/modules_pax_adult_Pos_Neg.xlsx")

###############################################################################
### MODELS FOR PATIENT CHARACTERISTICS & SYMPTOMS AMONG COVID+ - NP SAMPLES ###
###############################################################################

# Obesity
res_np_obesity <- as_tibble(dds_np_obesity)
res2_np_obesity <- res_np_obesity %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_obesity <- deframe(res2_np_obesity)
fgsea_np_obesity <- fgsea(pathways=modules_61, stats=ranks_np_obesity)
fgsea_np_obesity$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_obesity))) {fgsea_np_obesity[i,9]<-(length(unlist(fgsea_np_obesity[i,8])))/fgsea_np_obesity[i,7]}
fgsea_np_obesity <- data.frame(fgsea_np_obesity %>% 
                                  mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_obesity$leadingEdge <- NULL
write_xlsx(fgsea_np_obesity, "4_Symptoms/modules_np_obesity.xlsx")

# Asthma
res_np_asthma <- as_tibble(dds_np_asthma)
res2_np_asthma <- res_np_asthma %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_asthma <- deframe(res2_np_asthma)
fgsea_np_asthma <- fgsea(pathways=modules_61, stats=ranks_np_asthma)
fgsea_np_asthma$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_asthma))) {fgsea_np_asthma[i,9]<-(length(unlist(fgsea_np_asthma[i,8])))/fgsea_np_asthma[i,7]}
fgsea_np_asthma <- data.frame(fgsea_np_asthma %>% 
                                 mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_asthma$leadingEdge <- NULL
write_xlsx(fgsea_np_asthma, "4_Symptoms/modules_np_asthma.xlsx")

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
write_xlsx(fgsea_np_fever, "4_Symptoms/modules_np_fever.xlsx")

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
write_xlsx(fgsea_np_cough, "4_Symptoms/modules_np_cough.xlsx")

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
write_xlsx(fgsea_np_rhinorrhea, "4_Symptoms/modules_np_rhinorrhea.xlsx")

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
write_xlsx(fgsea_np_congestion, "4_Symptoms/modules_np_congestion.xlsx")

# Viral load
res_np_vl_copies <- as_tibble(dds_np_vl_copies)
res2_np_vl_copies <- res_np_vl_copies %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_np_vl_copies <- deframe(res2_np_vl_copies)
fgsea_np_vl_copies <- fgsea(pathways=modules_61, stats=ranks_np_vl_copies)
fgsea_np_vl_copies$GeneRatio <- 1
for(i in 1:(nrow(fgsea_np_vl_copies))) {fgsea_np_vl_copies[i,9]<-(length(unlist(fgsea_np_vl_copies[i,8])))/fgsea_np_vl_copies[i,7]}
fgsea_np_vl_copies <- data.frame(fgsea_np_vl_copies %>% 
                                    mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_np_vl_copies$leadingEdge <- NULL
write_xlsx(fgsea_np_vl_copies, "4_Symptoms/modules_np_vl_copies.xlsx")

################################################################################
### MODELS FOR PATIENT CHARACTERISTICS & SYMPTOMS AMONG COVID+ - PAX SAMPLES ###
################################################################################

# Obesity
res_pax_obesity <- as_tibble(dds_pax_obesity)
res2_pax_obesity <- res_pax_obesity %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_obesity <- deframe(res2_pax_obesity)
fgsea_pax_obesity <- fgsea(pathways=modules_61, stats=ranks_pax_obesity)
fgsea_pax_obesity$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_obesity))) {fgsea_pax_obesity[i,9]<-(length(unlist(fgsea_pax_obesity[i,8])))/fgsea_pax_obesity[i,7]}
fgsea_pax_obesity <- data.frame(fgsea_pax_obesity %>% 
                                mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_obesity$leadingEdge <- NULL
write_xlsx(fgsea_pax_obesity, "4_Symptoms/modules_pax_obesity.xlsx")

# Asthma
res_pax_asthma <- as_tibble(dds_pax_asthma)
res2_pax_asthma <- res_pax_asthma %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_asthma <- deframe(res2_pax_asthma)
fgsea_pax_asthma <- fgsea(pathways=modules_61, stats=ranks_pax_asthma)
fgsea_pax_asthma$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_asthma))) {fgsea_pax_asthma[i,9]<-(length(unlist(fgsea_pax_asthma[i,8])))/fgsea_pax_asthma[i,7]}
fgsea_pax_asthma <- data.frame(fgsea_pax_asthma %>% 
                                mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_asthma$leadingEdge <- NULL
write_xlsx(fgsea_pax_asthma, "4_Symptoms/modules_pax_asthma.xlsx")

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
write_xlsx(fgsea_pax_fever, "4_Symptoms/modules_pax_fever.xlsx")

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
write_xlsx(fgsea_pax_cough, "4_Symptoms/modules_pax_cough.xlsx")

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
write_xlsx(fgsea_pax_rhinorrhea, "4_Symptoms/modules_pax_rhinorrhea.xlsx")

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
write_xlsx(fgsea_pax_congestion, "4_Symptoms/modules_pax_congestion.xlsx")

# Viral load - high (>=10^6 copies/mL vs. <10^6 copies/mL)
res_pax_vl_copies <- as_tibble(dds_pax_vl_copies)
res2_pax_vl_copies <- res_pax_vl_copies %>% dplyr::select(Gene, stat) %>% na.omit() %>% distinct() %>% group_by(Gene) %>% dplyr::summarize(stat=mean(stat))
ranks_pax_vl_copies <- deframe(res2_pax_vl_copies)
fgsea_pax_vl_copies <- fgsea(pathways=modules_61, stats=ranks_pax_vl_copies)
fgsea_pax_vl_copies$GeneRatio <- 1
for(i in 1:(nrow(fgsea_pax_vl_copies))) {fgsea_pax_vl_copies[i,9]<-(length(unlist(fgsea_pax_vl_copies[i,8])))/fgsea_pax_vl_copies[i,7]}
fgsea_pax_vl_copies <- data.frame(fgsea_pax_vl_copies %>% 
                                     mutate(Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated", NES <= -1 & padj <= 0.05 ~ "Down-regulated", TRUE ~ "Unchanged")))
fgsea_pax_vl_copies$leadingEdge <- NULL
write_xlsx(fgsea_pax_vl_copies, "4_Symptoms/modules_pax_vl_copies.xlsx")