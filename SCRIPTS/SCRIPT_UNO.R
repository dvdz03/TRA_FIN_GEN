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

####
counts<-assay(datos1)
counts_filt<- counts[rowSums(counts)>10, ]
dge<- DGEList(counts_filt)
dge<- calcNormFactors(dge)

###
metadatos<- colData(datos1)
subtipos<- metadatos$paper_BRCA_Subtype_PAM50
genes_top<- counts_filt[order(apply(counts_filt, 1, var), decreasing= TRUE)[1:20], ]#tiene NA

#
head(subtipos)#hay una NA ahí
table(subtipos, useNA = "always")
muestras_buenas <- !is.na(subtipos)
genes_top_limpios <- genes_top [ ,muestras_buenas]
subtipos_limpios <- subtipos [muestras_buenas]
subtipos_factor <- as.factor (subtipos_limpios)

## VOLCANO PLOT
grupos <- ifelse (metadatos$sample_type == "Primary Tumor", "Tumor", "Normal")
dis <- model.matrix(~0 + grupos)
colnames(dis) <- c("Normal", "Tumor")
fit <- lmFit (counts_filt, dis)
contraste <- makeConstrasts(Tumor_vs_Normal = Tumor-Normal, levels = dis)#es de limma, falta instalar
