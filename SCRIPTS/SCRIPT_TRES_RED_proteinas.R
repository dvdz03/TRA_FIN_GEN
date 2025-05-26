##RED DE PROTEINAS, NO TIENE MUCHO SENTIDO BIOLOGICO

# Instalar y cargar paquetes necesarios
install.packages("STRINGdb")
library(STRINGdb)
library(igraph)

# Tu lista de símbolos de genes 
genes <- c("CPB1", "MT-RNR2", "RN7SL1", "IGHG1", "COL1A1", "MT-CO1", "MT-ND4", "FN1",
           "COL1A2", "COL3A1", "IGKC", "SCGB2A2", "MT-CO3", "MT-CO2", "MT-CYB",
           "IGFBP5", "CSN2", "LTF", "MT-ND2", "CHGA")

# Crear objeto STRINGdb para humano (species=9606)
string_db <- STRINGdb$new(version="11", species=9606, score_threshold=400, input_directory="")

# Mapear tus genes a IDs STRING
mapped <- string_db$map(data.frame(gene=genes), "gene", removeUnmappedRows=TRUE)

# Obtener interacciones conocidas entre tus genes
interactions <- string_db$get_interactions(mapped$STRING_id)

# Filtrar para solo interacciones entre tus genes
interactions <- interactions[interactions$from %in% mapped$STRING_id & interactions$to %in% mapped$STRING_id, ]

# Crear grafo igraph
red_proteinas_20_genes <- graph_from_data_frame(interactions[, c("from", "to")], directed=FALSE)

# Cambiar nombres de nodos a símbolos de genes
V(red_proteinas_20_genes)$name <- mapped$gene[match(V(red_proteinas_20_genes)$name, mapped$STRING_id)]

# Plot básico
plot(red_proteinas_20_genes,
     vertex.size=15,
     vertex.label.cex=0.8,
     vertex.color="lightblue",
     edge.width=1,
     main="Red de Interacciones Proteicas (STRING)")

# Guardar la figura si quieres
png("resultados/red_proteinas.png", width=800, height=600)
plot(red_proteinas_20_genes,
     vertex.size=15,
     vertex.label.cex=0.8,
     vertex.color="lightblue",
     edge.width=1,
     main="Red de Interacciones Proteicas (STRING)")
dev.off()


