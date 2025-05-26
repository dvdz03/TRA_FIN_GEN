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
counts<-assay(datos1) #extrae la matriz de cuentas (genes x muestras)
counts_filt<- counts[rowSums(counts)>10, ] # Filtrar los genes con baja expresión
dge<- DGEList(counts_filt)
dge<- calcNormFactors(dge)

###
metadatos<- colData(datos1) #METADATOS
subtipos<- metadatos$paper_BRCA_Subtype_PAM50 #Extrae los subtipos moleculares de BRCA
genes_top<- counts_filt[order(apply(counts_filt, 1, var), decreasing= TRUE)[1:20], ]#tiene NA

#
head(subtipos)#hay una NA ahí
table(subtipos, useNA = "always")
muestras_buenas <- !is.na(subtipos)
genes_top_limpios <- genes_top [ ,muestras_buenas]
subtipos_limpios <- subtipos [muestras_buenas]
subtipos_factor <- as.factor (subtipos_limpios)

#Cambiamos el codigo por los nombres de esos genes:
library(biomaRt)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Lista original con versión
ensg_ids <- c("ENSG00000153002.12", "ENSG00000210082.2", "ENSG00000276168.1", "ENSG00000211896.7",
              "ENSG00000108821.14", "ENSG00000198804.2", "ENSG00000198886.2", "ENSG00000115414.21",
              "ENSG00000164692.18", "ENSG00000168542.16", "ENSG00000211592.8", "ENSG00000110484.7",
              "ENSG00000198938.2", "ENSG00000198712.1", "ENSG00000198727.2", "ENSG00000115461.5",
              "ENSG00000135222.6", "ENSG00000012223.13", "ENSG00000198763.3", "ENSG00000100604.13")

# Quitar versión
ensg_ids_no_version <- sub("\\..*", "", ensg_ids)

# Obtener nombres de genes
gene_names <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                    filters = 'ensembl_gene_id',
                    values = ensg_ids_no_version,
                    mart = ensembl)

print(gene_names)
# Unir original con el resultado (ID sin versión con versión)
library(dplyr)
library(gridExtra)
df_final <- data.frame(original = ensg_ids,
                       ensembl_gene_id = ensg_ids_no_version,
                       stringsAsFactors = FALSE) %>%
  left_join(gene_names, by = "ensembl_gene_id") %>%
  select(original, ensembl_gene_id, hgnc_symbol)

# Guardar csv
write.csv(df_final, "resultados/los_20_genes_finales.csv", row.names = FALSE)

# Guardar tabla en png para incluir en R Markdown
png("resultados/los_20_genes_finales.png", width = 800, height = 400)
grid.table(df_final)
dev.off()


## VOLCANO PLOT
library(edgeR)
library(limma)
library(ggplot2)
library(ggrepel) 

# Matriz de cuentas ya filtrada: counts_filt
grupos <- ifelse(metadatos$sample_type == "Primary Tumor", "Tumor", "Normal")
grupos <- factor(grupos)
design <- model.matrix(~0 + grupos)
colnames(design) <- levels(grupos)

# Ajuste con limma
fit <- lmFit(counts_filt, design)
contraste <- makeContrasts(Tumor_vs_Normal = Tumor - Normal, levels = design)
fit2 <- contrasts.fit(fit, contraste)
fit2 <- eBayes(fit2)
resultados <- topTable(fit2, coef = "Tumor_vs_Normal", number = Inf)

# Etiquetar los 20 genes más significativos
resultados$gene <- rownames(resultados)
top20 <- resultados[order(resultados$adj.P.Val), ][1:20, ]

# Guardar el volcano plot con etiquetas
png("resultados/volcano_plot_BRCA.png", width = 1000, height = 800, res = 150)

ggplot(resultados, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = ifelse(abs(logFC) > 1 & adj.P.Val < 0.05, "Significativo", "No Significativo")), size = 2) +
  scale_color_manual(values = c("Significativo" = "red", "No Significativo" = "blue"), name = "Significancia") +
  labs(title = "Volcano Plot BRCA", x = "log2 FC", y = "-log10 p") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_text_repel(data = top20, aes(label = gene), size = 3, max.overlaps = 30) +
  theme_minimal()

dev.off()




grupos <- ifelse (metadatos$sample_type == "Primary Tumor", "Tumor", "Normal")
dis <- model.matrix(~0 + grupos)
colnames(dis) <- c("Normal", "Tumor")
fit <- lmFit (counts_filt, dis)
contraste <- makeContrasts(Tumor_vs_Normal = Tumor-Normal, levels = dis)
fit2 <- contrasts.fit(fit, contraste)
fit2 <- eBayes(fit2)
resultados <- topTable(fit2, number = Inf, coef = "Tumor_vs_Normal")
str(resultados)#esto es solo para ver si si tiene lo que se necesita para hacer el volcano, logFC...

png("resultados/volcano_plot_BRCA.png", width = 1000, height = 800, res = 150)

ggplot(resultados, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = ifelse(abs(logFC) > 1 & adj.P.Val < 0.05, "Significativo", "No Significativo")), size = 2) +
  scale_color_manual(values = c("Significativo" = "red", "No Significativo" = "blue"), name = "Significancia") +
  labs(title = "Volcano Plot BRCA", x = "log2 FC", y = "-log10 p") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black")

dev.off()







