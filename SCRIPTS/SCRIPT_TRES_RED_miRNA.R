library(multiMiR)
library(dplyr)
library(igraph)

genes <- c("CPB1", "MT-RNR2", "RN7SL1", "IGHG1", "COL1A1", "MT-CO1", "MT-ND4",
           "FN1", "COL1A2", "COL3A1", "IGKC", "SCGB2A2", "MT-CO3", "MT-CO2", 
           "MT-CYB", "IGFBP5", "CSN2", "LTF", "MT-ND2", "CHGA")

miRNAs_data <- get_multimir(org = "hsa", target = genes, table = "predicted")

# Forzar conversión a data.frame base
miRNA_df <- as.data.frame(miRNAs_data@data)

# Filtrar con subset base para evitar problemas con dplyr (score > 0.8)
miRNA_df_filtered <- subset(miRNA_df, score > 0.8)

# Seleccionar columnas usando base R para evitar error de select
edges <- unique(miRNA_df_filtered[, c("mature_mirna_acc", "target_symbol")])

# Renombrar columnas para igraph
colnames(edges) <- c("miRNA", "gene")

# Crear nodos con base R
nodes <- unique(c(edges$miRNA, edges$gene))
nodes_df <- data.frame(name = nodes)

# Crear grafo dirigido
g <- graph_from_data_frame(d = edges, vertices = nodes_df, directed = TRUE)

# Plotear

png("resultados/red_miRNA_genes.png", width = 1200, height = 900)
plot(g,
     layout = layout_with_fr,
     vertex.size = 10,
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     vertex.color = ifelse(V(g)$name %in% genes, "salmon", "lightblue"),
     edge.arrow.size = 0.3,
     main = "Red de interacciones miRNA - genes")

dev.off()
