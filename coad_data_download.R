# load library ------------------------------------------------------------

library(MultiAssayExperiment)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(tidyverse)


working_dir = "./"

setwd(working_dir)
genome <- "hg38"
data_dir = "../GDCdata"
# download transcription profiling data-----------------------------------------

query.rna <- GDCquery(project = "TCGA-COAD",
                          data.category = "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification",
                          workflow.type = "STAR - Counts",
                          sample.type = c("Primary Tumor","Solid Tissue Normal"),
                           legacy=FALSE)

GDCdownload(query.rna,directory=data_dir)


rna <-GDCprepare(query.rna,
                 directory=data_dir,
                save = TRUE,
                save.filename="./data/COAD_rna_hg38.rda")
