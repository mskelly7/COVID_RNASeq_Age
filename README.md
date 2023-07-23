# BRAVE-MESSI RNA Sequencing Analysis

Author: Aditya Mohan, Matthew Kelly <a href="https://orcid.org/0000-0001-8819-2315" target="orcid.widget" rel="noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon"></a>  
Date: July 23, 2023

This repository contains the files and code necessary to replicate the analyses presented in the manuscript '_Differential host responses within the upper respiratory tract and peripheral blood of children, adolescents, and adults with SARS-CoV-2 infection_', which has been submitted to medRxiv as a preprint (_________________). The overall objective of this manuscript was to evaluate associations between age and both local and systemic host responses to SARS-CoV-2. 

## Overview

This repository contains the data files ([`Data_Files/`](Data_Files/)), scripts ([`Scripts/`](Scripts/)), and outputs ([`Outputs/`](Outputs/)) for each of the RNA sequencing analyses presented in this manuscript. The raw RNA sequencing files are publicly available in the Gene Expression Omnibus ([GSE231409](https://www.ncbi.nlm.nih.gov/bioproject/?term=GSE231409)). 

## File Description

- [`Data_Files`](Data_Files/)/

  - [`BRAVE_RNASeq_counts.csv.gz`](Data_Files/BRAVE_RNASeq_counts.csv.gz): compressed file containing final gene count matrix
  - [`BRAVE_RNASeq_Metadata.csv`](Data_Files/BRAVE_RNASeq_Metadata.csv): metadata file
  - [`BRAVE_RNASeq_Data_Dictionary.xlsx`](Data_Files/BRAVE_RNASeq_Data_Dictionary.xlsx): data dictionary 
  - [`phy.rnaseq.np.rds`](Data_Files/phy.rnaseq.np.rds): phyloseq object containing metadata, gene count matrix, and gene annotations for NP samples included in analyses
  - [`phy.rnaseq.pax.rds`](Data_Files/phy.rnaseq.pax.rds): phyloseq object containing metadata, gene count matrix, and gene annotations for peripheral blood samples included in analyses

  - [`LM22.txt`](Data_Files/CIBERSORTx/LM22.txt): LM22 signature matrix 
  - [`LM22_source_GEPs.txt`](Data_Files/CIBERSORTx/LM22_source_GEPs.txt): source GEPs for LM22
  - [`CIBERSORTx_Results.txt`](Data_Files/CIBERSORTx/CIBERSORTx_Results.txt): output from CIBERSORTx containing imputed immune cell populations

  - [`nCounter_Host_Response_Gene_List.xlsx`](Data_Files/nCounter_Host_Response/nCounter_Host_Response_Gene_List.xlsx): nanoString Host Response Panel genes and annotations
  - [`modules_61.gmt`](Data_Files/nCounter_Host_Response/modules_61.gmt): list of genes for nanoString Host Response Panel modules (.gmt file)
  - [`modules_61.xlsx`](Data_Files/nCounter_Host_Response/modules_61.xlsx): list of genes for nanoString Host Response Panel modules (.xlsx file)

- [`Scripts`](Scripts/)/

  - [`RNASeq_Preprocessing.R`](Scripts/RNASeq_Preprocessing.R): preprocessing of gene count matrix and metadata; outputs: [`phy.rnaseq.np.rds`](Data_Files/phy.rnaseq.np.rds) and [`phy.rnaseq.pax.rds`](Data_Files/phy.rnaseq.pax.rds)
  - [`RNASeq_Manuscript_Text.R`](Scripts/RNASeq_Manuscript_Text.R): generation of results presented in manuscript text and Table 1
  - [`RNASeq_DESeq2.R`](Scripts/RNASeq_DESeq2.R): DESeq2 models for differential expression analyses of upper respiratory and peripheral blood samples
  - [`RNASeq_Modules.R`](Scripts/RNASeq_Modules.R): gene set enrichment analyses for upper respiratory and peripheral blood samples using FGSEA
  - [`RNASeq_Figure1.R`](Scripts/RNASeq_Figure1.R): Figure 1
  - [`RNASeq_Figure2.R`](Scripts/RNASeq_Figure2.R): Figure 2
  - [`RNASeq_Figure3.R`](Scripts/RNASeq_Figure3.R): Figure 3
  - [`RNASeq_Figure4.R`](Scripts/RNASeq_Figure4.R): Figure 4
  - [`RNASeq_Figures5_&_S5.R`](Scripts/RNASeq_Figures5_&_S5.R): Figure 5 and Supplemental Figure 5
  - [`RNASeq_FigureS1.R`](Scripts/RNASeq_FigureS1.R): Supplemental Figure 1
  - [`RNASeq_FigureS2.R`](Scripts/RNASeq_FigureS2.R): Supplemental Figure 2
  - [`RNASeq_FigureS3.R`](Scripts/RNASeq_FigureS3.R): Supplemental Figure 3
  - [`RNASeq_FigureS4.R`](Scripts/RNASeq_FigureS4.R): Supplemental Figure 4

- [`Output_Files`](Output_Files/)/

  - [`1_COVID_Neg_by_Age`](Output_Files/1_COVID_Neg_by_Age/): results from DESeq2 and FGSEA for analyses of upper respiratory and peripheral blood samples among SARS-CoV-2-uninfected subjects
  
