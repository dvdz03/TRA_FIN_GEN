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
install.packages("limma")
library(limma)
BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(survival)
library(enrichplot)
library(clusterProfiler)
library(ggplot2)
library(DESeq2)
library(multiMiR)
library(dplyr)
library(factoextra)
library(FactoMineR)
library(genefilter)
library(gProfileR)
library(igraph)
library(lattice)
library(modelr)
library(networkD3)
library(org.Hs.eg.db)
library(patchwork)
library(pheatmap)
library(tidyverse)

#las otras fueron cargadas desde PACKAGES

#CARGA DE LOS DATOS TCGA
query1<-GDCquery(project = "TCGA-BRCA",
                 data.category= "Transcriptome Profiling",
                 data.type= "Gene Expression Quantification",
                 workflow.type = "STAR - Counts",
                 sample.type = c("Primary Tumor", "Solid Tissue Normal"),
                 experimental.strategy = "RNA-Seq")
GDCdownload(query1, files.per.chunk = 50)
datos1<-GDCprepare(query1)

