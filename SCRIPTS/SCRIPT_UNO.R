#########################################################
##SCRIPT PARA HACER TODA LA CARGA DE DATOS Y ESAS COSAS##
#########################################################

#CARGAR TODAS LAS LIBRERÍAS Y PAQUETES
install.packages("BiocManager")
library(BiocManager)
install.packages("data.table")
library(data.table)
BiocManager::install("biomaRt")
library(biomaRt)
BiocManager::install("edgeR")
library(edgeR)
library(limma)
BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(survival)
library(enrichplot)
#las otras fueron cargadas desde PACKAGES

#CARGA DE LOS DATOS TCGA
query1<-GDCquery(project = "TCGA-BRCA",
                 data.category= "Transcriptome Profiling",
                 data.type= "Gene Expression Quantification",
                 workflow.type = "STAR - Counts",
                 sample.type = c("Primary Tumor", "Solid Tissue Normal"),
                 experimental.strategy = "RNA-Seq")
