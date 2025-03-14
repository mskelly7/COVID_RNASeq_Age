# BRAVE-MESSI RNA Sequencing Analysis

Author: Aditya Mohan, Matthew Kelly <a href="https://orcid.org/0000-0001-8819-2315" target="orcid.widget" rel="noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon"></a>  
Last update: February 22, 2025

This repository contains the files and code necessary to replicate the analyses presented in the manuscript '_Age-associated differences in mucosal and systemic host responses to SARS-CoV-2 infection_', which has been published in Nature Communications (https://rdcu.be/ecWQQ). The overall objective of this manuscript was to evaluate associations between age and both local and systemic host responses to SARS-CoV-2. 

## Overview

This repository contains the data files ([`Data_Files/`](Data_Files/)), scripts ([`Scripts/`](Scripts/)), and outputs ([`Outputs/`](Outputs/)) for each of the RNA sequencing analyses presented in this manuscript. The raw RNA sequencing files are publicly available in the Gene Expression Omnibus ([GSE231409](https://www.ncbi.nlm.nih.gov/bioproject/?term=GSE231409)). The sequencing dataset used for SARS-CoV-2 lineage assignment is available in the Sequence Read Archive ([PRJNA1024980](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA1024980)). 

## File Description

- [`Data_Files`](Data_Files/)/

  - [`BRAVE_RNASeq_counts.csv.gz`](Data_Files/BRAVE_RNASeq_counts.csv.gz): compressed file containing gene count matrix (all samples)
  - [`BRAVE_RNASeq_Metadata.csv`](Data_Files/BRAVE_RNASeq_Metadata.csv): metadata file (all samples)
  - [`BRAVE_RNASeq_Data_Dictionary.xlsx`](Data_Files/BRAVE_RNASeq_Data_Dictionary.xlsx): data dictionary defining variables contained in metadata file
  - [`phy.rnaseq.np.rds`](Data_Files/phy.rnaseq.np.rds): phyloseq object containing metadata, gene count matrix, and gene annotations for upper respiratory samples included in manuscript analyses
  - [`phy.rnaseq.pax.rds`](Data_Files/phy.rnaseq.pax.rds): phyloseq object containing metadata, gene count matrix, and gene annotations for peripheral blood samples included in manuscript analyses
  - [`LM22.txt`](Data_Files/CIBERSORTx/LM22.txt): LM22 signature matrix (required for CIBERSORTx submission)
  - [`LM22_source_GEPs.txt`](Data_Files/CIBERSORTx/LM22_source_GEPs.txt): source GEPs for LM22 (required for CIBERSORTx submission)
  - [`CIBERSORTx_Results.txt`](Data_Files/CIBERSORTx/CIBERSORTx_Results.txt): CIBERSORTx output containing imputed immune cell populations for upper respiratory and peripheral blood samples included in manuscript analyses
  - [`nCounter_Host_Response_Gene_List.xlsx`](Data_Files/nCounter_Host_Response/nCounter_Host_Response_Gene_List.xlsx): nanoString Host Response Panel genes and annotations
  - [`modules_61.gmt`](Data_Files/nCounter_Host_Response/modules_61.gmt): list of genes for nanoString Host Response Panel modules (.gmt file)
  - [`modules_61.xlsx`](Data_Files/nCounter_Host_Response/modules_61.xlsx): list of genes for nanoString Host Response Panel modules (.xlsx file)

- [`Scripts`](Scripts/)/

  - [`RNASeq_Preprocessing.R`](Scripts/RNASeq_Preprocessing.R): processing of raw gene count matrix and metadata; outputs: [`phy.rnaseq.np.rds`](Data_Files/phy.rnaseq.np.rds), [`phy.rnaseq.pax.rds`](Data_Files/phy.rnaseq.pax.rds)
  - [`RNASeq_Manuscript_Text.R`](Scripts/RNASeq_Manuscript_Text.R): results presented in manuscript text and Table 1; inputs: [`phy.rnaseq.np.rds`](Data_Files/phy.rnaseq.np.rds), [`phy.rnaseq.pax.rds`](Data_Files/phy.rnaseq.pax.rds)
  - [`RNASeq_DESeq2.R`](Scripts/RNASeq_DESeq2.R): DESeq2 models for differential expression analyses; inputs: [`phy.rnaseq.np.rds`](Data_Files/phy.rnaseq.np.rds), [`phy.rnaseq.pax.rds`](Data_Files/phy.rnaseq.pax.rds)
  - [`RNASeq_Modules.R`](Scripts/RNASeq_Modules.R): gene set enrichment analyses using FGSEA; inputs: [`modules_61.gmt`](Data_Files/nCounter_Host_Response/modules_61.gmt), DESeq2 ("genes") output files
  - [`RNASeq_Figure1.R`](Scripts/RNASeq_Figure1.R): Figure 1; inputs: [`phy.rnaseq.np.rds`](Data_Files/phy.rnaseq.np.rds), [`phy.rnaseq.pax.rds`](Data_Files/phy.rnaseq.pax.rds), FGSEA ("modules") output files in [`1_COVID_Neg_by_Age`](Output_Files/1_COVID_Neg_by_Age/)
  - [`RNASeq_Figure2.R`](Scripts/RNASeq_Figure2.R): Figure 2; inputs: DESeq2 ("genes") and FGSEA ("modules") output files in [`2_COVID_Pos_vs_Neg`](Output_Files/2_COVID_Pos_vs_Neg/)
  - [`RNASeq_Figure3.R`](Scripts/RNASeq_Figure3.R): Figure 3; inputs: [`phy.rnaseq.np.rds`](Data_Files/phy.rnaseq.np.rds), [`phy.rnaseq.pax.rds`](Data_Files/phy.rnaseq.pax.rds), FGSEA ("modules") output files in [`3_COVID_Pos_by_Age`](Output_Files/3_COVID_Pos_by_Age/)
  - [`RNASeq_Figure4.R`](Scripts/RNASeq_Figure4.R): Figure 4; inputs: FGSEA ("modules") output files in [`4_Symptoms`](Output_Files/4_Symptoms/)
  - [`RNASeq_Figure5.R`](Scripts/RNASeq_Figure5.R): Figure 5; inputs: [`phy.rnaseq.np.rds`](Data_Files/phy.rnaseq.np.rds), [`phy.rnaseq.pax.rds`](Data_Files/phy.rnaseq.pax.rds), [`modules_61.xlsx`](Data_Files/nCounter_Host_Response/modules_61.xlsx)
  - [`RNASeq_FigureS1.R`](Scripts/RNASeq_FigureS1.R): Supplemental Figure 1; inputs: DESeq2 ("genes") output files in [`1_COVID_Neg_by_Age`](Output_Files/1_COVID_Neg_by_Age/)
  - [`RNASeq_FigureS2.R`](Scripts/RNASeq_FigureS2.R): Supplemental Figure 2; inputs: [`phy.rnaseq.np.rds`](Data_Files/phy.rnaseq.np.rds), [`phy.rnaseq.pax.rds`](Data_Files/phy.rnaseq.pax.rds)
  - [`RNASeq_FigureS3.R`](Scripts/RNASeq_FigureS3.R): Supplemental Figure 3; inputs: DESeq2 ("genes") output files in [`3_COVID_Pos_by_Age`](Output_Files/3_COVID_Pos_by_Age/)
  - [`RNASeq_FigureS4.R`](Scripts/RNASeq_FigureS4.R): Supplemental Figure 4; inputs: DESeq2 ("genes") output files in [`3_COVID_Pos_by_Age`](Output_Files/3_COVID_Pos_by_Age/)

- [`Output_Files`](Output_Files/)/ 

  - [`1_COVID_Neg_by_Age`](Output_Files/1_COVID_Neg_by_Age/): output files from analyses of SARS-CoV-2-uninfected subjects by age category
  - [`2_COVID_Pos_vs_Neg`](Output_Files/2_COVID_Pos_vs_Neg/): output files from analyses comparing SARS-CoV-2-infected and uninfected subjects
  - [`3_COVID_Pos_by_Age`](Output_Files/3_COVID_Pos_by_Age/): output files from analyses of SARS-CoV-2-infected subjects by age
  - [`4_Symptoms`](Output_Files/4_Symptoms/): output files from analyses of comorbidities and symptoms among SARS-CoV-2-infected subjects
  - [`5_NP_vs_PAX_Correlations`](Output_Files/5_NP_vs_PAX_Correlations/): output files from analyses evaluating correlations between upper respiratory and peripheral blood gene set expression
  - [`Figures`](Output_Files/Figures/): figures and supplemental figures
