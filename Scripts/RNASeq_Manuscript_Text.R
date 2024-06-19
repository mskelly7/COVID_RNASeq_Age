# BRAVE Kids RNA Sequencing Analysis
# Aditya Mohan (MD/PhD candidate) / Matthew Kelly, MD, MPH 
# Results presented in manuscript text and Tables 1, S1, and S2
# Last update: June 18, 2024

remove(list=ls())
setwd("____________________________") 
set.seed(1234)

library(readxl)
library(openxlsx)
library(tidyverse)
library(dplyr)
library(plyr)
library(phyloseq)
library(data.table)
library(reshape2)

###############################################################
# UPLOAD METADATA FILE FOR SUMMARIES OF SUBJECT CHARACTERISTICS
###############################################################

# Load the required datasets
phy.rnaseq.np <- readRDS("phy.rnaseq.np.rds")
metadata_np <- data.frame(sample_data(phy.rnaseq.np))
summary(metadata_np$age)
phy.rnaseq.pax <- readRDS("phy.rnaseq.pax.rds")
metadata_pax <- data.frame(sample_data(phy.rnaseq.pax))
metadata_pax_brave <- subset(metadata_pax, age<21)
summary(metadata_pax_brave$age)
metadata_pax_messi <- subset(metadata_pax, age>=21)
summary(metadata_pax_messi$age)
metadata_all <- rbind(metadata_np, metadata_pax)
metadata_all <- metadata_all[!duplicated(cbind(metadata_all$alias_study_id,metadata_all$SampleTiming)),]
metadata_brave <- rbind(metadata_np, metadata_pax_brave)
metadata_brave <- metadata_brave[!duplicated(metadata_brave$alias_study_id),]
# Note that there is one adult who provided both acute infection and healthy samples
metadata_covid <- subset(metadata_all, corona=="Positive")
metadata_covid_vl <- subset(metadata_all, corona=="Positive" & vl_copies>0)

#######################
# ABSTRACT/INTRODUCTION
#######################

nrow(metadata_all)
summary(metadata_all$age)
table(metadata_all$corona)

#########
# RESULTS
#########

#########################################
# Characteristics of the study population
#########################################

# BRAVE Kids subjects
nrow(metadata_np)
nrow(metadata_pax_brave)
nrow(metadata_brave)
summary(metadata_brave$age)
table(metadata_brave$sex)
prop.table(table(metadata_brave$sex))
table(metadata_brave$age_cat)
table(metadata_brave$corona)
prop.table(table(metadata_brave$corona))
table(metadata_brave$group)
tapply(metadata_brave$timing_sx, metadata_brave$group, summary)
tapply(metadata_brave$timing_dx, metadata_brave$group, summary)

# MESSI participants
nrow(metadata_pax_messi)
table(metadata_pax_messi$corona)
summary(metadata_pax_messi$age)
table(metadata_pax_messi$sex)
prop.table(table(metadata_pax_messi$sex))
table(metadata_pax_messi$corona)
prop.table(table(metadata_pax_messi$corona))
table(metadata_pax_messi$group)
tapply(metadata_pax_messi$timing_sx, metadata_pax_messi$group, summary)
tapply(metadata_pax_messi$timing_dx, metadata_pax_messi$group, summary)

# Symptoms among SARS-CoV-2-infected by age group
pos_sx_metadata <- subset(metadata_all, group=="POS_SX")
tapply(pos_sx_metadata$timing_sx, pos_sx_metadata$age_cat, summary)
kruskal.test(timing_sx ~ age_cat, data=pos_sx_metadata)
# No difference in timing of sample collection relative to symptom onset in POS_SX
pos_asx_metadata <- subset(metadata_all, group=="POS_ASX")
tapply(pos_asx_metadata$timing_dx, pos_asx_metadata$age_cat, summary)
kruskal.test(timing_dx ~ age_cat, data=pos_asx_metadata)
# No difference in timing of sample collection relative to SARS-CoV-2 diagnosis in POS_ASX

# Comparing adults/adolescents to young/school-age children
metadata_covid$age_cat2[metadata_covid$age_cat=="0-5 years" | metadata_covid$age_cat=="6-13 years"] <- "Children"
metadata_covid$age_cat2[metadata_covid$age_cat=="14-20 years" | metadata_covid$age_cat=="Adult"] <- "Adolescents/Adults"
prop.table(table(metadata_covid$age_cat2, metadata_covid$cough),1)
chisq.test(metadata_covid$age_cat2, metadata_covid$cough)
table(metadata_covid$age_cat2, metadata_covid$headache)
prop.table(table(metadata_covid$age_cat2, metadata_covid$headache),1)
chisq.test(metadata_covid$age_cat2, metadata_covid$headache)
table(metadata_covid$age_cat2, metadata_covid$myalgias)
prop.table(table(metadata_covid$age_cat2, metadata_covid$myalgias),1)
chisq.test(metadata_covid$age_cat2, metadata_covid$myalgias)
table(metadata_covid$age_cat2, metadata_covid$congestion)
prop.table(table(metadata_covid$age_cat2, metadata_covid$congestion),1)
chisq.test(metadata_covid$age_cat2, metadata_covid$congestion)
table(metadata_covid$age_cat2, metadata_covid$rhinorrhea)
prop.table(table(metadata_covid$age_cat2, metadata_covid$rhinorrhea),1)
chisq.test(metadata_covid$age_cat2, metadata_covid$rhinorrhea)
table(metadata_covid$age_cat2, metadata_covid$anosmia)
prop.table(table(metadata_covid$age_cat2, metadata_covid$anosmia),1)
chisq.test(metadata_covid$age_cat2, metadata_covid$anosmia)
table(metadata_covid$age_cat2, metadata_covid$dysgeusia)
prop.table(table(metadata_covid$age_cat2, metadata_covid$dysgeusia),1)
chisq.test(metadata_covid$age_cat2, metadata_covid$dysgeusia)
table(metadata_pax_messi$corona, metadata_pax_messi$comorbidity_oth)
prop.table(table(metadata_pax_messi$corona, metadata_pax_messi$comorbidity_oth),1)
fisher.test(metadata_pax_messi$corona, metadata_pax_messi$comorbidity_oth)

metadata_covid_vl$vl_copies <- log10(metadata_covid_vl$vl_copies)
table(metadata_covid_vl$age_cat)
tapply(metadata_covid_vl$vl_copies, metadata_covid_vl$age_cat, summary)
kruskal.test(vl_copies ~ age_cat, data=metadata_covid_vl)
# No difference in viral load by age category 
table(metadata_covid_vl$group)
wilcox.test(vl_copies ~ group, data=metadata_covid_vl)
# No difference in VL in symptomatic and asymptomatic SARS-CoV-2-infected subjects 

table(metadata_covid$sars2_variant)
# Only ancestral strains identified by genomic sequencing

##############################################################################################
# The mucosal and systemic transcriptional profiles of healthy individuals differ based on age
##############################################################################################

metadata_no_covid_brave <- subset(metadata_all, corona=="Negative" & age_cat!="Adult")
table(metadata_no_covid_brave$age_cat)
metadata_pax_messi_no_covid <- subset(metadata_pax_messi, corona=="Negative")
nrow(metadata_pax_messi_no_covid)

# Use beta regression to evaluate for changes in cellular proportions associated with age
# Limit to analyses of cell populations that are present in >=25% of samples by sample type

library(betareg)
library(lmtest)
library(rcompanion)
library(broom)

# UPPER RESPIRATORY SAMPLES (COVID-NEGATIVE BY AGE)
# Adjust comparisons for multiple testing using BH (n=10)
# B cells: 0.74
# Plasma cells: >0.99
# CD8+ T cells: >0.99
# CD4+ T cells: >0.99
# NK cells: >0.99
# Mono/macrophages: >0.99
# Dendritic cells: 0.27
# Mast cells: >0.99
# Eosinophils: >0.99
# Neutrophils: >0.99

phy.rnaseq.np.Neg <- subset_samples(phy.rnaseq.np, corona=="Negative" & age_cat!="Adult")
metadata_np_neg <- data.frame(sample_data(phy.rnaseq.np.Neg))
metadata_np_neg_pops <- metadata_np_neg[,c("age","age_cat","B_cells","Plasma_cells", "T_cells_CD8","T_cells_CD4","T_cells_gamma_delta","NK_cells",
                                             "Mono_Macrophages","Dendritic_cells","Mast_cells","Eosinophils","Neutrophils")]

# Exclude cell populations with <25% of samples containing this immune cell population
summary(metadata_np_neg_pops$T_cells_gamma_delta)

metadata_np_neg_pops$B_cells[metadata_np_neg_pops$B_cells==0] <- 0.000001
B_cells_model_np <- betareg(B_cells ~ age, data=metadata_np_neg_pops)
summary(B_cells_model_np)
B_cells_model_np <- data.frame(tidy(coeftest(B_cells_model_np)))
p.adjust(B_cells_model_np[2,5], method="fdr", n=10)

metadata_np_neg_pops$Plasma_cells[metadata_np_neg_pops$Plasma_cells==0] <- 0.000001
Plasma_cells_model_np <- betareg(Plasma_cells ~ age, data=metadata_np_neg_pops)
summary(Plasma_cells_model_np)
Plasma_cells_model_np <- data.frame(tidy(coeftest(Plasma_cells_model_np)))
p.adjust(Plasma_cells_model_np[2,5], method="fdr", n=10)

metadata_np_neg_pops$T_cells_CD8[metadata_np_neg_pops$T_cells_CD8==0] <- 0.000001
T_cells_CD8_model_np <- betareg(T_cells_CD8 ~ age, data=metadata_np_neg_pops)
summary(T_cells_CD8_model_np)
T_cells_CD8_model_np <- data.frame(tidy(coeftest(T_cells_CD8_model_np)))
p.adjust(T_cells_CD8_model_np[2,5], method="fdr", n=10)

metadata_np_neg_pops$T_cells_CD4[metadata_np_neg_pops$T_cells_CD4==0] <- 0.000001
T_cells_CD4_model_np <- betareg(T_cells_CD4 ~ age, data=metadata_np_neg_pops)
summary(T_cells_CD4_model_np)
T_cells_CD4_model_np <- data.frame(tidy(coeftest(T_cells_CD4_model_np)))
p.adjust(T_cells_CD4_model_np[2,5], method="fdr", n=10)

metadata_np_neg_pops$NK_cells[metadata_np_neg_pops$NK_cells==0] <- 0.000001
NK_cells_model_np <- betareg(NK_cells ~ age, data=metadata_np_neg_pops)
summary(NK_cells_model_np)
NK_cells_model_np <- data.frame(tidy(coeftest(NK_cells_model_np)))
p.adjust(NK_cells_model_np[2,5], method="fdr", n=10)

metadata_np_neg_pops$Mono_Macrophages[metadata_np_neg_pops$Mono_Macrophages==0] <- 0.000001
Mono_Macrophages_model_np <- betareg(Mono_Macrophages ~ age, data=metadata_np_neg_pops)
summary(Mono_Macrophages_model_np)
Mono_Macrophages_model_np <- data.frame(tidy(coeftest(Mono_Macrophages_model_np)))
p.adjust(Mono_Macrophages_model_np[2,5], method="fdr", n=10)

metadata_np_neg_pops$Dendritic_cells[metadata_np_neg_pops$Dendritic_cells==0] <- 0.000001
Dendritic_cells_model_np <- betareg(Dendritic_cells ~ age, data=metadata_np_neg_pops)
summary(Dendritic_cells_model_np)
Dendritic_cells_model_np <- data.frame(tidy(coeftest(Dendritic_cells_model_np)))
p.adjust(Dendritic_cells_model_np[2,5], method="fdr", n=10)

metadata_np_neg_pops$Mast_cells[metadata_np_neg_pops$Mast_cells==0] <- 0.000001
Mast_cells_model_np <- betareg(Mast_cells ~ age, data=metadata_np_neg_pops)
summary(Mast_cells_model_np)
Mast_cells_model_np <- data.frame(tidy(coeftest(Mast_cells_model_np)))
p.adjust(Mast_cells_model_np[2,5], method="fdr", n=10)

metadata_np_neg_pops$Eosinophils[metadata_np_neg_pops$Eosinophils==0] <- 0.000001
Eosinophils_model_np <- betareg(Eosinophils ~ age, data=metadata_np_neg_pops)
summary(Eosinophils_model_np)
Eosinophils_model_np <- data.frame(tidy(coeftest(Eosinophils_model_np)))
p.adjust(Eosinophils_model_np[2,5], method="fdr", n=10)

metadata_np_neg_pops$Neutrophils[metadata_np_neg_pops$Neutrophils==0] <- 0.000001
Neutrophils_model_np <- betareg(Neutrophils ~ age, data=metadata_np_neg_pops)
summary(Neutrophils_model_np)
Neutrophils_model_np <- data.frame(tidy(coeftest(Neutrophils_model_np)))
p.adjust(Neutrophils_model_np[2,5], method="fdr", n=10)

# PERIPHERAL BLOOD SAMPLES (COVID-NEGATIVE BY AGE)
# Adjust comparisons for multiple testing using BH (n=8)
# B cells: 9.03E-14
# Plasma cells: 0.005
# CD8+ T cells: 0.002
# CD4+ T cells: >0.99
# NK cells: 0.69
# Mono/macrophages: 1.95E-5
# Mast cells: 0.15
# Neutrophils: 0.00094

phy.rnaseq.pax.Neg <- subset_samples(phy.rnaseq.pax, corona=="Negative")
metadata_pax_neg <- data.frame(sample_data(phy.rnaseq.pax.Neg))
metadata_pax_neg_pops <- metadata_pax_neg[,c("age","age_cat","B_cells","Plasma_cells", "T_cells_CD8","T_cells_CD4","T_cells_gamma_delta","NK_cells",
                                             "Mono_Macrophages","Dendritic_cells","Mast_cells","Eosinophils","Neutrophils")]

# Exclude cell populations with <25% of samples containing this immune cell population
summary(metadata_pax_neg_pops$T_cells_gamma_delta)
summary(metadata_pax_neg_pops$Dendritic_cells)
summary(metadata_pax_neg_pops$Eosinophils)

metadata_pax_neg_pops$B_cells[metadata_pax_neg_pops$B_cells==0] <- 0.000001
B_cells_model_pax <- betareg(B_cells ~ age, data=metadata_pax_neg_pops)
summary(B_cells_model_pax)
B_cells_model_pax <- data.frame(tidy(coeftest(B_cells_model_pax)))
p.adjust(B_cells_model_pax[2,5], method="fdr", n=8)
# Decreasing proportion of B cells with increasing age (9.03E-14)

metadata_pax_neg_pops$Plasma_cells[metadata_pax_neg_pops$Plasma_cells==0] <- 0.000001
Plasma_cells_model_pax <- betareg(Plasma_cells ~ age, data=metadata_pax_neg_pops)
summary(Plasma_cells_model_pax)
Plasma_cells_model_pax <- data.frame(tidy(coeftest(Plasma_cells_model_pax)))
p.adjust(Plasma_cells_model_pax[2,5], method="fdr", n=8)
# Decreasing proportion of plasma cells with increasing age (p=0.005)

metadata_pax_neg_pops$T_cells_CD8[metadata_pax_neg_pops$T_cells_CD8==0] <- 0.000001
T_cells_CD8_model_pax <- betareg(T_cells_CD8 ~ age, data=metadata_pax_neg_pops)
summary(T_cells_CD8_model_pax)
T_cells_CD8_model_pax <- data.frame(tidy(coeftest(T_cells_CD8_model_pax)))
p.adjust(T_cells_CD8_model_pax[2,5], method="fdr", n=8)
# Decreasing proportion of CD8 T cells with increasing age (p=0.002)

metadata_pax_neg_pops$T_cells_CD4[metadata_pax_neg_pops$T_cells_CD4==0] <- 0.000001
T_cells_CD4_model_pax <- betareg(T_cells_CD4 ~ age, data=metadata_pax_neg_pops)
summary(T_cells_CD4_model_pax)
T_cells_CD4_model_pax <- data.frame(tidy(coeftest(T_cells_CD4_model_pax)))
p.adjust(T_cells_CD4_model_pax[2,5], method="fdr", n=8)

metadata_pax_neg_pops$NK_cells[metadata_pax_neg_pops$NK_cells==0] <- 0.000001
NK_cells_model_pax <- betareg(NK_cells ~ age, data=metadata_pax_neg_pops)
summary(NK_cells_model_pax)
NK_cells_model_pax <- data.frame(tidy(coeftest(NK_cells_model_pax)))
p.adjust(NK_cells_model_pax[2,5], method="fdr", n=8)

metadata_pax_neg_pops$Mono_Macrophages[metadata_pax_neg_pops$Mono_Macrophages==0] <- 0.000001
Mono_Macrophages_model_pax <- betareg(Mono_Macrophages ~ age, data=metadata_pax_neg_pops)
summary(Mono_Macrophages_model_pax)
Mono_Macrophages_model_pax <- data.frame(tidy(coeftest(Mono_Macrophages_model_pax)))
p.adjust(Mono_Macrophages_model_pax[2,5], method="fdr", n=8)
# Increasing proportion of monophages and macrophages with increasing age (p=1.95E-5)

metadata_pax_neg_pops$Mast_cells[metadata_pax_neg_pops$Mast_cells==0] <- 0.000001
Mast_cells_model_pax <- betareg(Mast_cells ~ age, data=metadata_pax_neg_pops)
summary(Mast_cells_model_pax)
Mast_cells_model_pax <- data.frame(tidy(coeftest(Mast_cells_model_pax)))
p.adjust(Mast_cells_model_pax[2,5], method="fdr", n=8)

metadata_pax_neg_pops$Neutrophils[metadata_pax_neg_pops$Neutrophils==0] <- 0.000001
Neutrophils_model_pax <- betareg(Neutrophils ~ age, data=metadata_pax_neg_pops)
summary(Neutrophils_model_pax)
Neutrophils_model_pax <- data.frame(tidy(coeftest(Neutrophils_model_pax)))
p.adjust(Neutrophils_model_pax[2,5], method="fdr", n=8)
# Increasing proportion of neutrophils with increasing age (p=0.00094)

########################################################################################
# Distinct mucosal and systemic immune responses characterize early SARS-CoV-2 infection
########################################################################################

# UPPER RESPIRATORY SAMPLES (COVID-POSITIVE VS. COVID-NEGATIVE)
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

metadata_np_pops <- metadata_np[,c("age","age_cat","corona","B_cells","Plasma_cells", "T_cells_CD8","T_cells_CD4","T_cells_gamma_delta","NK_cells",
                                           "Mono_Macrophages","Dendritic_cells","Mast_cells","Eosinophils","Neutrophils")]

metadata_np_pops$B_cells[metadata_np_pops$B_cells==0] <- 0.000001
B_cells_model_np <- betareg(B_cells ~ age + corona, data=metadata_np_pops)
summary(B_cells_model_np)
B_cells_model_np <- data.frame(tidy(coeftest(B_cells_model_np)))
p.adjust(B_cells_model_np[3,5], method="fdr", n=11)

metadata_np_pops$Plasma_cells[metadata_np_pops$Plasma_cells==0] <- 0.000001
Plasma_cells_model_np <- betareg(Plasma_cells ~ age + corona, data=metadata_np_pops)
summary(Plasma_cells_model_np)
Plasma_cells_model_np <- data.frame(tidy(coeftest(Plasma_cells_model_np)))
p.adjust(Plasma_cells_model_np[3,5], method="fdr", n=11)

metadata_np_pops$T_cells_CD8[metadata_np_pops$T_cells_CD8==0] <- 0.000001
T_cells_CD8_model_np <- betareg(T_cells_CD8 ~ age + corona, data=metadata_np_pops)
summary(T_cells_CD8_model_np)
T_cells_CD8_model_np <- data.frame(tidy(coeftest(T_cells_CD8_model_np)))
p.adjust(T_cells_CD8_model_np[3,5], method="fdr", n=11)

metadata_np_pops$T_cells_CD4[metadata_np_pops$T_cells_CD4==0] <- 0.000001
T_cells_CD4_model_np <- betareg(T_cells_CD4 ~ age + corona, data=metadata_np_pops)
summary(T_cells_CD4_model_np)
T_cells_CD4_model_np <- data.frame(tidy(coeftest(T_cells_CD4_model_np)))
p.adjust(T_cells_CD4_model_np[3,5], method="fdr", n=11)

metadata_np_pops$T_cells_gamma_delta[metadata_np_pops$T_cells_gamma_delta==0] <- 0.000001
T_cells_gamma_delta_model_np <- betareg(T_cells_gamma_delta ~ age + corona, data=metadata_np_pops)
summary(T_cells_gamma_delta_model_np)
T_cells_gamma_delta_model_np <- data.frame(tidy(coeftest(T_cells_gamma_delta_model_np)))
p.adjust(T_cells_gamma_delta_model_np[3,5], method="fdr", n=11)

metadata_np_pops$NK_cells[metadata_np_pops$NK_cells==0] <- 0.000001
NK_cells_model_np <- betareg(NK_cells ~ age + corona, data=metadata_np_pops)
summary(NK_cells_model_np)
NK_cells_model_np <- data.frame(tidy(coeftest(NK_cells_model_np)))
p.adjust(NK_cells_model_np[3,5], method="fdr", n=11)

metadata_np_pops$Mono_Macrophages[metadata_np_pops$Mono_Macrophages==0] <- 0.000001
Mono_Macrophages_model_np <- betareg(Mono_Macrophages ~ age + corona, data=metadata_np_pops)
summary(Mono_Macrophages_model_np)
Mono_Macrophages_model_np <- data.frame(tidy(coeftest(Mono_Macrophages_model_np)))
p.adjust(Mono_Macrophages_model_np[3,5], method="fdr", n=11)
# Increase in monophages and macrophages associated with SARS-CoV-2 (p=0.03)

metadata_np_pops$Dendritic_cells[metadata_np_pops$Dendritic_cells==0] <- 0.000001
Dendritic_cells_model_np <- betareg(Dendritic_cells ~ age + corona, data=metadata_np_pops)
summary(Dendritic_cells_model_np)
Dendritic_cells_model_np <- data.frame(tidy(coeftest(Dendritic_cells_model_np)))
p.adjust(Dendritic_cells_model_np[3,5], method="fdr", n=11)
# Decrease in dendritic cells associated with SARS-CoV-2 (p=0.02)

metadata_np_pops$Mast_cells[metadata_np_pops$Mast_cells==0] <- 0.000001
Mast_cells_model_np <- betareg(Mast_cells ~ age + corona, data=metadata_np_pops)
summary(Mast_cells_model_np)
Mast_cells_model_np <- data.frame(tidy(coeftest(Mast_cells_model_np)))
p.adjust(Mast_cells_model_np[3,5], method="fdr", n=11)

metadata_np_pops$Eosinophils[metadata_np_pops$Eosinophils==0] <- 0.000001
Eosinophils_model_np <- betareg(Eosinophils ~ age + corona, data=metadata_np_pops)
summary(Eosinophils_model_np)
Eosinophils_model_np <- data.frame(tidy(coeftest(Eosinophils_model_np)))
p.adjust(Eosinophils_model_np[3,5], method="fdr", n=11)

metadata_np_pops$Neutrophils[metadata_np_pops$Neutrophils==0] <- 0.000001
Neutrophils_model_np <- betareg(Neutrophils ~ age + corona, data=metadata_np_pops)
summary(Neutrophils_model_np)
Neutrophils_model_np <- data.frame(tidy(coeftest(Neutrophils_model_np)))
p.adjust(Neutrophils_model_np[3,5], method="fdr", n=11)

# PERIPHERAL BLOOD SAMPLES (COVID-POSITIVE VS. COVID-NEGATIVE)
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

metadata_pax_pops <- metadata_pax[,c("age","age_cat","corona","B_cells","Plasma_cells", "T_cells_CD8","T_cells_CD4","T_cells_gamma_delta","NK_cells",
                                   "Mono_Macrophages","Dendritic_cells","Mast_cells","Eosinophils","Neutrophils")]

# Exclude cell populations with <25% of samples containing this immune cell population
summary(metadata_pax_pops$T_cells_gamma_delta)
summary(metadata_pax_pops$Eosinophils)

# Adjust comparisons for multiple testing using BH (n=9)

metadata_pax_pops$B_cells[metadata_pax_pops$B_cells==0] <- 0.000001
B_cells_model_pax <- betareg(B_cells ~ age + corona, data=metadata_pax_pops)
summary(B_cells_model_pax)
B_cells_model_pax <- data.frame(tidy(coeftest(B_cells_model_pax)))
p.adjust(B_cells_model_pax[3,5], method="fdr", n=9)

metadata_pax_pops$Plasma_cells[metadata_pax_pops$Plasma_cells==0] <- 0.000001
Plasma_cells_model_pax <- betareg(Plasma_cells ~ age + corona, data=metadata_pax_pops)
summary(Plasma_cells_model_pax)
Plasma_cells_model_pax <- data.frame(tidy(coeftest(Plasma_cells_model_pax)))
p.adjust(Plasma_cells_model_pax[3,5], method="fdr", n=9)
# Increase in plasma cells associated with SARS-CoV-2 (p=0.03)

metadata_pax_pops$T_cells_CD8[metadata_pax_pops$T_cells_CD8==0] <- 0.000001
T_cells_CD8_model_pax <- betareg(T_cells_CD8 ~ age + corona, data=metadata_pax_pops)
summary(T_cells_CD8_model_pax)
T_cells_CD8_model_pax <- data.frame(tidy(coeftest(T_cells_CD8_model_pax)))
p.adjust(T_cells_CD8_model_pax[3,5], method="fdr", n=9)

metadata_pax_pops$T_cells_CD4[metadata_pax_pops$T_cells_CD4==0] <- 0.000001
T_cells_CD4_model_pax <- betareg(T_cells_CD4 ~ age + corona, data=metadata_pax_pops)
summary(T_cells_CD4_model_pax)
T_cells_CD4_model_pax <- data.frame(tidy(coeftest(T_cells_CD4_model_pax)))
p.adjust(T_cells_CD4_model_pax[3,5], method="fdr", n=9)

metadata_pax_pops$NK_cells[metadata_pax_pops$NK_cells==0] <- 0.000001
NK_cells_model_pax <- betareg(NK_cells ~ age + corona, data=metadata_pax_pops)
summary(NK_cells_model_pax)
NK_cells_model_pax <- data.frame(tidy(coeftest(NK_cells_model_pax)))
p.adjust(NK_cells_model_pax[3,5], method="fdr", n=9)

metadata_pax_pops$Mono_Macrophages[metadata_pax_pops$Mono_Macrophages==0] <- 0.000001
Mono_Macrophages_model_pax <- betareg(Mono_Macrophages ~ age + corona, data=metadata_pax_pops)
summary(Mono_Macrophages_model_pax)
Mono_Macrophages_model_pax <- data.frame(tidy(coeftest(Mono_Macrophages_model_pax)))
p.adjust(Mono_Macrophages_model_pax[3,5], method="fdr", n=9)

metadata_pax_pops$Dendritic_cells[metadata_pax_pops$Dendritic_cells==0] <- 0.000001
Dendritic_cells_model_pax <- betareg(Dendritic_cells ~ age + corona, data=metadata_pax_pops)
summary(Dendritic_cells_model_pax)
Dendritic_cells_model_pax <- data.frame(tidy(coeftest(Dendritic_cells_model_pax)))
p.adjust(Dendritic_cells_model_pax[3,5], method="fdr", n=9)

metadata_pax_pops$Mast_cells[metadata_pax_pops$Mast_cells==0] <- 0.000001
Mast_cells_model_pax <- betareg(Mast_cells ~ age + corona, data=metadata_pax_pops)
summary(Mast_cells_model_pax)
Mast_cells_model_pax <- data.frame(tidy(coeftest(Mast_cells_model_pax)))
p.adjust(Mast_cells_model_pax[3,5], method="fdr", n=9)

metadata_pax_pops$Neutrophils[metadata_pax_pops$Neutrophils==0] <- 0.000001
Neutrophils_model_pax <- betareg(Neutrophils ~ age + corona, data=metadata_pax_pops)
summary(Neutrophils_model_pax)
Neutrophils_model_pax <- data.frame(tidy(coeftest(Neutrophils_model_pax)))
p.adjust(Neutrophils_model_pax[3,5], method="fdr", n=9)

##########################################################################################
# Mucosal immune responses to SARS-CoV-2 infection are similar across pediatric age groups
##########################################################################################

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

phy.rnaseq.np.Pos <- subset_samples(phy.rnaseq.np, corona=="Positive" & age_cat!="Adult")
metadata_np_pos <- data.frame(sample_data(phy.rnaseq.np.Pos))
metadata_np_pos_pops <- metadata_np_pos[,c("age","age_cat","B_cells","Plasma_cells", "T_cells_CD8","T_cells_CD4","T_cells_gamma_delta","NK_cells",
                                           "Mono_Macrophages","Dendritic_cells","Mast_cells","Eosinophils","Neutrophils")]

# Exclude cell populations with <25% of samples containing this immune cell population
summary(metadata_np_pos_pops$Eosinophils)
summary(metadata_np_pos_pops$T_cells_gamma_delta)

# Adjust comparisons for multiple testing using BH (n=9)

metadata_np_pos_pops$B_cells[metadata_np_pos_pops$B_cells==0] <- 0.000001
B_cells_model_np <- betareg(B_cells ~ age, data=metadata_np_pos_pops)
summary(B_cells_model_np)
B_cells_model_np <- data.frame(tidy(coeftest(B_cells_model_np)))
p.adjust(B_cells_model_np[2,5], method="fdr", n=9)

metadata_np_pos_pops$Plasma_cells[metadata_np_pos_pops$Plasma_cells==0] <- 0.000001
Plasma_cells_model_np <- betareg(Plasma_cells ~ age, data=metadata_np_pos_pops)
summary(Plasma_cells_model_np)
Plasma_cells_model_np <- data.frame(tidy(coeftest(Plasma_cells_model_np)))
p.adjust(Plasma_cells_model_np[2,5], method="fdr", n=9)

metadata_np_pos_pops$T_cells_CD8[metadata_np_pos_pops$T_cells_CD8==0] <- 0.000001
T_cells_CD8_model_np <- betareg(T_cells_CD8 ~ age, data=metadata_np_pos_pops)
summary(T_cells_CD8_model_np)
T_cells_CD8_model_np <- data.frame(tidy(coeftest(T_cells_CD8_model_np)))
p.adjust(T_cells_CD8_model_np[2,5], method="fdr", n=9)

metadata_np_pos_pops$T_cells_CD4[metadata_np_pos_pops$T_cells_CD4==0] <- 0.000001
T_cells_CD4_model_np <- betareg(T_cells_CD4 ~ age, data=metadata_np_pos_pops)
summary(T_cells_CD4_model_np)
T_cells_CD4_model_np <- data.frame(tidy(coeftest(T_cells_CD4_model_np)))
p.adjust(T_cells_CD4_model_np[2,5], method="fdr", n=9)

metadata_np_pos_pops$NK_cells[metadata_np_pos_pops$NK_cells==0] <- 0.000001
NK_cells_model_np <- betareg(NK_cells ~ age, data=metadata_np_pos_pops)
summary(NK_cells_model_np)
NK_cells_model_np <- data.frame(tidy(coeftest(NK_cells_model_np)))
p.adjust(NK_cells_model_np[2,5], method="fdr", n=9)

metadata_np_pos_pops$Mono_Macrophages[metadata_np_pos_pops$Mono_Macrophages==0] <- 0.000001
Mono_Macrophages_model_np <- betareg(Mono_Macrophages ~ age, data=metadata_np_pos_pops)
summary(Mono_Macrophages_model_np)
Mono_Macrophages_model_np <- data.frame(tidy(coeftest(Mono_Macrophages_model_np)))
p.adjust(Mono_Macrophages_model_np[2,5], method="fdr", n=9)

metadata_np_pos_pops$Dendritic_cells[metadata_np_pos_pops$Dendritic_cells==0] <- 0.000001
Dendritic_cells_model_np <- betareg(Dendritic_cells ~ age, data=metadata_np_pos_pops)
summary(Dendritic_cells_model_np)
Dendritic_cells_model_np <- data.frame(tidy(coeftest(Dendritic_cells_model_np)))
p.adjust(Dendritic_cells_model_np[2,5], method="fdr", n=9)

metadata_np_pos_pops$Mast_cells[metadata_np_pos_pops$Mast_cells==0] <- 0.000001
Mast_cells_model_np <- betareg(Mast_cells ~ age, data=metadata_np_pos_pops)
summary(Mast_cells_model_np)
Mast_cells_model_np <- data.frame(tidy(coeftest(Mast_cells_model_np)))
p.adjust(Mast_cells_model_np[2,5], method="fdr", n=9)

metadata_np_pos_pops$Eosinophils[metadata_np_pos_pops$Eosinophils==0] <- 0.000001
Eosinophils_model_np <- betareg(Eosinophils ~ age, data=metadata_np_pos_pops)
summary(Eosinophils_model_np)
Eosinophils_model_np <- data.frame(tidy(coeftest(Eosinophils_model_np)))
p.adjust(Eosinophils_model_np[2,5], method="fdr", n=9)

metadata_np_pos_pops$Neutrophils[metadata_np_pos_pops$Neutrophils==0] <- 0.000001
Neutrophils_model_np <- betareg(Neutrophils ~ age, data=metadata_np_pos_pops)
summary(Neutrophils_model_np)
Neutrophils_model_np <- data.frame(tidy(coeftest(Neutrophils_model_np)))
p.adjust(Neutrophils_model_np[2,5], method="fdr", n=9)

# COMPARE RESPONSES TO SARS-CoV-2 ACROSS AGE GROUPS

# Upload module information 
modules <- read_excel("Statistical_Analyses/modules_61.xlsx", sheet="modules_61", na="NA")
library(GSVA)

# ssGSEA data for NP samples
genes_np <- data.frame(tax_table(phy.rnaseq.np))
counts_np <- data.frame(otu_table(phy.rnaseq.np))
data_np <- merge(genes_np, counts_np, by="row.names")
data_np <- data_np[,-c(1,3,4)] # Drop Row.names, Chromosome, Length
nrow(data_np)
# Sum counts for genes that have multiple EnsemblIDs
data_np <- ddply(data_np, "Gene", numcolwise(sum))
nrow(data_np)
row.names(data_np) <- data_np$Gene
data_np$Gene <- NULL
data_np <- as.matrix(data_np)
data_np_log2 <- log2(data_np)
ssgsea_np <- gsva(data_np_log2, modules, method="ssgsea")
ssgsea_np <- as.data.frame(t(ssgsea_np))
colnames(ssgsea_np) <- paste(colnames(ssgsea_np),"NSB",sep=".")
ssgsea_np <- merge(metadata_np, ssgsea_np, by="row.names")
colnames(ssgsea_np) <- gsub(" ", "_", colnames(ssgsea_np))
colnames(ssgsea_np) <- gsub("-", "_", colnames(ssgsea_np))

# Reference = 0-5 years
ssgsea_np$age_cat <- relevel(factor(ssgsea_np$age_cat), ref="0-5 years")
summary(lm(Interferon_Response.NSB ~ batch_num + age_cat + corona + age_cat:corona + sex, data = ssgsea_np)) 
summary(lm(Type_I_Interferon_Signaling.NSB ~ batch_num + age_cat + corona + age_cat:corona + sex, data = ssgsea_np))
summary(lm(Type_II_Interferon_Signaling.NSB ~ batch_num + age_cat + corona + age_cat:corona + sex, data = ssgsea_np))
summary(lm(TNF_Signaling.NSB ~ batch_num + age_cat + corona + age_cat:corona + sex, data = ssgsea_np))
summary(lm(TLR_Signaling.NSB ~ batch_num + age_cat + corona + age_cat:corona + sex, data = ssgsea_np))
summary(lm(Phagocytosis.NSB ~ batch_num + age_cat + corona + age_cat:corona + sex, data = ssgsea_np))
summary(lm(Myeloid_Activation.NSB ~ batch_num + age_cat + corona + age_cat:corona + sex, data = ssgsea_np)) 
summary(lm(Myeloid_Inflammation.NSB ~ batch_num + age_cat + corona + age_cat:corona + sex, data = ssgsea_np))
summary(lm(Inflammasomes.NSB ~ batch_num + age_cat + corona + age_cat:corona + sex, data = ssgsea_np))
summary(lm(NK_Activity.NSB ~ batch_num + age_cat + corona + age_cat:corona + sex, data = ssgsea_np))
summary(lm(Complement_System.NSB ~ batch_num + age_cat + corona + age_cat:corona + sex, data = ssgsea_np))
summary(lm(Coagulation.NSB ~ batch_num + age_cat + corona + age_cat:corona + sex, data = ssgsea_np)) 
summary(lm(Adaptive_Immune_Response.NSB ~ batch_num + age_cat + corona + age_cat:corona + sex, data = ssgsea_np))
summary(lm(BCR_Signaling.NSB ~ batch_num + age_cat + corona + age_cat:corona + sex, data = ssgsea_np)) 
summary(lm(TCR_Signaling.NSB ~ batch_num + age_cat + corona + age_cat:corona + sex, data = ssgsea_np))
summary(lm(Mononuclear_Cell_Migration.NSB ~ batch_num + age_cat + corona + age_cat:corona + sex, data = ssgsea_np)) 
summary(lm(Lymphocyte_Trafficking.NSB ~ batch_num + age_cat + corona + age_cat:corona + sex, data = ssgsea_np)) 
summary(lm(Treg_Differentiation.NSB ~ batch_num + age_cat + corona + age_cat:corona + sex, data = ssgsea_np))
# No significant differences in expression of these modules in NP samples

##################################################################################################
# Systemic immune responses to SARS-CoV-2 infection differ among children, adolescents, and adults
##################################################################################################

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

phy.rnaseq.pax.Pos <- subset_samples(phy.rnaseq.pax, corona=="Positive")
metadata_pax_pos <- data.frame(sample_data(phy.rnaseq.pax.Pos))
metadata_pax_pos_pops <- metadata_pax_pos[,c("age","age_cat","B_cells","Plasma_cells", "T_cells_CD8","T_cells_CD4","T_cells_gamma_delta","NK_cells",
                                             "Mono_Macrophages","Dendritic_cells","Mast_cells","Eosinophils","Neutrophils")]

# Exclude cell populations with <25% of samples containing this immune cell population
summary(metadata_pax_pos_pops$Eosinophils)
summary(metadata_pax_pos_pops$T_cells_gamma_delta)

metadata_pax_pos_pops$B_cells[metadata_pax_pos_pops$B_cells==0] <- 0.000001
B_cells_model_pax <- betareg(B_cells ~ age, data=metadata_pax_pos_pops)
summary(B_cells_model_pax)
B_cells_model_pax <- data.frame(tidy(coeftest(B_cells_model_pax)))
p.adjust(B_cells_model_pax[2,5], method="fdr", n=9)
# Decreasing proportion of B cells with increasing age (p=6.66E-19)

metadata_pax_pos_pops$Plasma_cells[metadata_pax_pos_pops$Plasma_cells==0] <- 0.000001
Plasma_cells_model_pax <- betareg(Plasma_cells ~ age, data=metadata_pax_pos_pops)
summary(Plasma_cells_model_pax)
Plasma_cells_model_pax <- data.frame(tidy(coeftest(Plasma_cells_model_pax)))
p.adjust(Plasma_cells_model_pax[2,5], method="fdr", n=9)
# Decreasing proportion of plasma cells with increasing age (p=0.0011)

metadata_pax_pos_pops$T_cells_CD8[metadata_pax_pos_pops$T_cells_CD8==0] <- 0.000001
T_cells_CD8_model_pax <- betareg(T_cells_CD8 ~ age, data=metadata_pax_pos_pops)
summary(T_cells_CD8_model_pax)
T_cells_CD8_model_pax <- data.frame(tidy(coeftest(T_cells_CD8_model_pax)))
p.adjust(T_cells_CD8_model_pax[2,5], method="fdr", n=9)
# Decreasing proportion of CD8+ T cells with increasing age (p=0.0003)

metadata_pax_pos_pops$T_cells_CD4[metadata_pax_pos_pops$T_cells_CD4==0] <- 0.000001
T_cells_CD4_model_pax <- betareg(T_cells_CD4 ~ age, data=metadata_pax_pos_pops)
summary(T_cells_CD4_model_pax)
T_cells_CD4_model_pax <- data.frame(tidy(coeftest(T_cells_CD4_model_pax)))
p.adjust(T_cells_CD4_model_pax[2,5], method="fdr", n=9)

metadata_pax_pos_pops$NK_cells[metadata_pax_pos_pops$NK_cells==0] <- 0.000001
NK_cells_model_pax <- betareg(NK_cells ~ age, data=metadata_pax_pos_pops)
summary(NK_cells_model_pax)
NK_cells_model_pax <- data.frame(tidy(coeftest(NK_cells_model_pax)))
p.adjust(NK_cells_model_pax[2,5], method="fdr", n=9)
# Increasing proportion of NK cells with increasing age (p=0.02)

metadata_pax_pos_pops$Mono_Macrophages[metadata_pax_pos_pops$Mono_Macrophages==0] <- 0.000001
Mono_Macrophages_model_pax <- betareg(Mono_Macrophages ~ age, data=metadata_pax_pos_pops)
summary(Mono_Macrophages_model_pax)
Mono_Macrophages_model_pax <- data.frame(tidy(coeftest(Mono_Macrophages_model_pax)))
p.adjust(Mono_Macrophages_model_pax[2,5], method="fdr", n=9)
# Increasing proportion of monophages and macrophages with increasing age (p=1.52E-11)

metadata_pax_pos_pops$Dendritic_cells[metadata_pax_pos_pops$Dendritic_cells==0] <- 0.000001
Dendritic_cells_model_pax <- betareg(Dendritic_cells ~ age, data=metadata_pax_pos_pops)
summary(Dendritic_cells_model_pax)
Dendritic_cells_model_pax <- data.frame(tidy(coeftest(Dendritic_cells_model_pax)))
p.adjust(Dendritic_cells_model_pax[2,5], method="fdr", n=9)

metadata_pax_pos_pops$Mast_cells[metadata_pax_pos_pops$Mast_cells==0] <- 0.000001
Mast_cells_model_pax <- betareg(Mast_cells ~ age, data=metadata_pax_pos_pops)
summary(Mast_cells_model_pax)
Mast_cells_model_pax <- data.frame(tidy(coeftest(Mast_cells_model_pax)))
p.adjust(Mast_cells_model_pax[2,5], method="fdr", n=9)

metadata_pax_pos_pops$Neutrophils[metadata_pax_pos_pops$Neutrophils==0] <- 0.000001
Neutrophils_model_pax <- betareg(Neutrophils ~ age, data=metadata_pax_pos_pops)
summary(Neutrophils_model_pax)
Neutrophils_model_pax <- data.frame(tidy(coeftest(Neutrophils_model_pax)))
p.adjust(Neutrophils_model_pax[2,5], method="fdr", n=9)
# Increasing proportion of neutrophils with increasing age (p=2.08E-05)

# COMPARE RESPONSES TO SARS-CoV-2 ACROSS AGE GROUPS

# ssGSEA data for PAX samples
genes_pax <- data.frame(tax_table(phy.rnaseq.pax))
counts_pax <- data.frame(otu_table(phy.rnaseq.pax))
data_pax <- merge(genes_pax, counts_pax, by="row.names")
data_pax <- data_pax[,-c(1,3,4)] # Drop Row.names, Chromosome, Length
nrow(data_pax)
# Sum counts for genes that have multiple EnsemblIDs
data_pax <- ddply(data_pax, "Gene", numcolwise(sum))
nrow(data_pax)
row.names(data_pax) <- data_pax$Gene
data_pax$Gene <- NULL
data_pax <- as.matrix(data_pax)
data_pax_log2 <- log2(data_pax)
ssgsea_pax <- gsva(data_pax_log2, modules, method="ssgsea")
ssgsea_pax <- as.data.frame(t(ssgsea_pax))
colnames(ssgsea_pax) <- paste(colnames(ssgsea_pax),"PAX",sep=".")
ssgsea_pax <- merge(metadata_pax, ssgsea_pax, by="row.names")
colnames(ssgsea_pax) <- gsub(" ", "_", colnames(ssgsea_pax))
colnames(ssgsea_pax) <- gsub("-", "_", colnames(ssgsea_pax))

# Reference = 0-5 years
ssgsea_pax$age_cat <- relevel(factor(ssgsea_pax$age_cat), ref="0-5 years")
summary(lm(Interferon_Response.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax)) 
summary(lm(Type_I_Interferon_Signaling.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax))
summary(lm(Type_II_Interferon_Signaling.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax))
summary(lm(TNF_Signaling.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax))
summary(lm(TLR_Signaling.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax))
summary(lm(Phagocytosis.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax))
# Young children with less upregulation of genes for phagocytosis
summary(lm(Myeloid_Activation.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax)) 
# Young children with less upregulation of genes for myeloid activation
summary(lm(Myeloid_Inflammation.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax))
summary(lm(Inflammasomes.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax))
summary(lm(NK_Activity.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax))
summary(lm(Complement_System.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax))
summary(lm(Coagulation.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax)) 
summary(lm(Adaptive_Immune_Response.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax))
summary(lm(BCR_Signaling.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax)) 
summary(lm(TCR_Signaling.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax))
summary(lm(Mononuclear_Cell_Migration.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax)) 
summary(lm(Lymphocyte_Trafficking.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax)) 
# Young children with less upregulation of genes for lymphocyte trafficking
summary(lm(Treg_Differentiation.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax))

# Reference = Adult
ssgsea_pax$age_cat <- relevel(factor(ssgsea_pax$age_cat), ref="Adult")
summary(lm(Interferon_Response.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax)) 
summary(lm(Type_I_Interferon_Signaling.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax))
summary(lm(Type_II_Interferon_Signaling.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax))
summary(lm(TNF_Signaling.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax))
summary(lm(TLR_Signaling.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax))
summary(lm(Phagocytosis.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax))
summary(lm(Myeloid_Activation.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax)) 
summary(lm(Myeloid_Inflammation.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax))
summary(lm(Inflammasomes.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax))
summary(lm(NK_Activity.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax))
summary(lm(Complement_System.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax))
# Adults have greater upregulation of genes corresponding to the complement system
summary(lm(Coagulation.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax)) 
summary(lm(Adaptive_Immune_Response.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax))
summary(lm(BCR_Signaling.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax)) 
summary(lm(TCR_Signaling.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax))
# Adults with less upregulation of genes for T cell receptor signaling
summary(lm(Mononuclear_Cell_Migration.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax)) 
summary(lm(Lymphocyte_Trafficking.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax)) 
summary(lm(Treg_Differentiation.PAX ~ batch_num + age_cat + corona + age_cat:corona + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = ssgsea_pax))
# Adults with less upregulation of genes for regulatory T cell differentiation

############
# DISCUSSION
############

nrow(metadata_all)
table(metadata_all$corona)

#########
# METHODS
#########

summary(metadata_np$RIN)
summary(metadata_pax$RIN)

phy.rnaseq.pruned <- readRDS("phy.rnaseq.pruned.rds")
phy.rnaseq.np.raw <- subset_samples(phy.rnaseq.pruned, SampleType=="np")
nsamples(phy.rnaseq.np.raw)
summary(sample_sums(phy.rnaseq.np.raw))
phy.rnaseq.pax.raw <- subset_samples(phy.rnaseq.pruned, SampleType=="pax")
nsamples(phy.rnaseq.pax.raw)
summary(sample_sums(phy.rnaseq.pax.raw))

metadata_pax_pca <- metadata_pax[,c("B_cells", "Plasma_cells", "T_cells_CD4", "T_cells_CD8", "T_cells_gamma_delta", "NK_cells",
                                    "Mono_Macrophages", "Dendritic_cells", "Mast_cells", "Eosinophils", "Neutrophils")]
summary(metadata_pax_pca$T_cells_gamma_delta)
metadata_pax_pca$T_cells_gamma_delta <- NULL
summary(metadata_pax_pca$Eosinophils)
metadata_pax_pca$Eosinophils <- NULL
pax_pca <- prcomp(metadata_pax_pca, center = TRUE, scale = TRUE)
print(pax_pca)
summary(pax_pca)

#########
# TABLE 1
#########

nrow(metadata_all)
table(metadata_all$age_cat)
tapply(metadata_all$age, metadata_all$age_cat, summary)
kruskal.test(age ~ age_cat, data=metadata_all)
table(metadata_all$age_cat, metadata_all$sex)
prop.table(table(metadata_all$age_cat, metadata_all$sex),1)
chisq.test(table(metadata_all$age_cat, metadata_all$sex))
table(metadata_all$age_cat, metadata_all$obesity)
prop.table(table(metadata_all$age_cat, metadata_all$obesity),1)
chisq.test(table(metadata_all$age_cat, metadata_all$obesity))
table(metadata_all$age_cat, metadata_all$obesity, useNA="always")
table(metadata_all$age_cat, metadata_all$comorbidity_oth)
prop.table(table(metadata_all$age_cat, metadata_all$comorbidity_oth),1)
fisher.test(table(metadata_all$age_cat, metadata_all$comorbidity_oth))
table(metadata_all$age_cat, metadata_all$comorbidity_oth, useNA="always")
table(metadata_all$age_cat, metadata_all$corona)
prop.table(table(metadata_all$age_cat, metadata_all$corona),1)
chisq.test(table(metadata_all$age_cat, metadata_all$corona))
table(metadata_covid$age_cat, metadata_covid$symptoms)
prop.table(table(metadata_covid$age_cat, metadata_covid$symptoms),1)
fisher.test(table(metadata_covid$age_cat, metadata_covid$symptoms))
table(metadata_covid$age_cat, metadata_covid$fever)
prop.table(table(metadata_covid$age_cat, metadata_covid$fever),1)
chisq.test(table(metadata_covid$age_cat, metadata_covid$fever))
table(metadata_covid$age_cat, metadata_covid$cough)
prop.table(table(metadata_covid$age_cat, metadata_covid$cough),1)
chisq.test(table(metadata_covid$age_cat, metadata_covid$cough))
table(metadata_covid$age_cat, metadata_covid$headache)
prop.table(table(metadata_covid$age_cat, metadata_covid$headache),1)
fisher.test(table(metadata_covid$age_cat, metadata_covid$headache))
table(metadata_covid$age_cat, metadata_covid$abd_pain)
prop.table(table(metadata_covid$age_cat, metadata_covid$abd_pain),1)
fisher.test(table(metadata_covid$age_cat, metadata_covid$abd_pain))
table(metadata_covid$age_cat, metadata_covid$myalgias)
prop.table(table(metadata_covid$age_cat, metadata_covid$myalgias),1)
fisher.test(table(metadata_covid$age_cat, metadata_covid$myalgias))
table(metadata_covid$age_cat, metadata_covid$congestion)
prop.table(table(metadata_covid$age_cat, metadata_covid$congestion),1)
fisher.test(table(metadata_covid$age_cat, metadata_covid$congestion))
table(metadata_covid$age_cat, metadata_covid$rhinorrhea)
prop.table(table(metadata_covid$age_cat, metadata_covid$rhinorrhea),1)
fisher.test(table(metadata_covid$age_cat, metadata_covid$rhinorrhea))
table(metadata_covid$age_cat, metadata_covid$anosmia)
prop.table(table(metadata_covid$age_cat, metadata_covid$anosmia),1)
fisher.test(table(metadata_covid$age_cat, metadata_covid$anosmia))
table(metadata_covid$age_cat, metadata_covid$dysgeusia)
prop.table(table(metadata_covid$age_cat, metadata_covid$dysgeusia),1)
fisher.test(table(metadata_covid$age_cat, metadata_covid$dysgeusia))
tapply(metadata_covid_vl$vl_copies, metadata_covid_vl$age_cat, summary)
kruskal.test(vl_copies ~ age_cat, data=metadata_covid_vl)

##########
# TABLE S1
##########

table(metadata_all$corona)
table(metadata_all$age_cat, metadata_all$corona)

metadata_0to5 <- subset(metadata_all, age_cat=="0-5 years")
tapply(metadata_0to5$age, metadata_0to5$corona, summary)
wilcox.test(age ~ corona, data=metadata_0to5)
table(metadata_0to5$sex, metadata_0to5$corona, useNA="always")
prop.table(table(metadata_0to5$sex, metadata_0to5$corona),2)
fisher.test(metadata_0to5$sex, metadata_0to5$corona)
table(metadata_0to5$obesity, metadata_0to5$corona, useNA="always")
prop.table(table(metadata_0to5$obesity, metadata_0to5$corona),2)
fisher.test(metadata_0to5$obesity, metadata_0to5$corona)
table(metadata_0to5$comorbidity_oth, metadata_0to5$corona, useNA="always")
prop.table(table(metadata_0to5$comorbidity_oth, metadata_0to5$corona),2)
fisher.test(metadata_0to5$comorbidity_oth, metadata_0to5$corona)

metadata_6to13 <- subset(metadata_all, age_cat=="6-13 years")
tapply(metadata_6to13$age, metadata_6to13$corona, summary)
wilcox.test(age ~ corona, data=metadata_6to13)
table(metadata_6to13$sex, metadata_6to13$corona, useNA="always")
prop.table(table(metadata_6to13$sex, metadata_6to13$corona),2)
fisher.test(metadata_6to13$sex, metadata_6to13$corona)
table(metadata_6to13$obesity, metadata_6to13$corona, useNA="always")
prop.table(table(metadata_6to13$obesity, metadata_6to13$corona),2)
fisher.test(metadata_6to13$obesity, metadata_6to13$corona)
table(metadata_6to13$comorbidity_oth, metadata_6to13$corona, useNA="always")
prop.table(table(metadata_6to13$comorbidity_oth, metadata_6to13$corona),2)
fisher.test(metadata_6to13$comorbidity_oth, metadata_6to13$corona)

metadata_14to20 <- subset(metadata_all, age_cat=="14-20 years")
tapply(metadata_14to20$age, metadata_14to20$corona, summary)
wilcox.test(age ~ corona, data=metadata_14to20)
table(metadata_14to20$sex, metadata_14to20$corona, useNA="always")
prop.table(table(metadata_14to20$sex, metadata_14to20$corona),2)
fisher.test(metadata_14to20$sex, metadata_14to20$corona)
table(metadata_14to20$obesity, metadata_14to20$corona, useNA="always")
prop.table(table(metadata_14to20$obesity, metadata_14to20$corona),2)
fisher.test(metadata_14to20$obesity, metadata_14to20$corona)
table(metadata_14to20$comorbidity_oth, metadata_14to20$corona, useNA="always")
prop.table(table(metadata_14to20$comorbidity_oth, metadata_14to20$corona),2)
fisher.test(metadata_14to20$comorbidity_oth, metadata_14to20$corona)

metadata_adult <- subset(metadata_all, age_cat=="Adult")
tapply(metadata_adult$age, metadata_adult$corona, summary)
wilcox.test(age ~ corona, data=metadata_adult)
table(metadata_adult$sex, metadata_adult$corona, useNA="always")
prop.table(table(metadata_adult$sex, metadata_adult$corona),2)
fisher.test(metadata_adult$sex, metadata_adult$corona)
table(metadata_adult$obesity, metadata_adult$corona, useNA="always")
prop.table(table(metadata_adult$obesity, metadata_adult$corona),2)
fisher.test(metadata_adult$obesity, metadata_adult$corona)
table(metadata_adult$comorbidity_oth, metadata_adult$corona, useNA="always")
prop.table(table(metadata_adult$comorbidity_oth, metadata_adult$corona),2)
fisher.test(metadata_adult$comorbidity_oth, metadata_adult$corona)

##########
# TABLE S2
##########

table(metadata_covid$sars2_variant)