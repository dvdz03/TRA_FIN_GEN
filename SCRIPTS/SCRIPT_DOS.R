#################################
#### ESTE ES PARA EL HEATMAP ####
#################################

install.packages("pheatmap")
library(pheatmap)


paraheat <- data.frame(Subtype = subtipos_factor)
rownames(paraheat) <- colnames(genes_top_limpios)
pheatmap(mat = genes_top_limpios,
         annotation_col = paraheat,
         show_colnames = FALSE,
         color = colorRampPalette(c("orange","black","darkorchid"))(100))

# Anotación por subtipo
paraheat <- data.frame(Subtype = subtipos_factor)
rownames(paraheat) <- colnames(genes_top_limpios)

# Guardar el heatmap como PNG
png("resultados/heatmap_subtipos.png", width = 1200, height = 1000)
pheatmap(mat = genes_top_limpios,
         annotation_col = paraheat,
         show_colnames = FALSE,
         color = colorRampPalette(c("orange","black","darkorchid"))(100))
dev.off()


# Guardar imagen, verificando errores, ya que el primero salia como "archivo dañado"
png("resultados/heatmap_subtipos.png", width = 1200, height = 1000)

tryCatch({
  pheatmap(mat = genes_top_limpios,
           annotation_col = paraheat,
           show_colnames = FALSE,
           color = colorRampPalette(c("orange","black","darkorchid"))(100))
}, error = function(e) {
  dev.off()  # Cierra el dispositivo en caso de error
  stop("Error al generar el heatmap: ", e$message)
})

dev.off()


#Para intepretar estos resulatdos, tambien sacamos tablas de:
# Asegurarse de que la matriz esté en formato data.frame
df_genes <- as.data.frame(genes_top_limpios)
df_genes$Gene <- rownames(genes_top_limpios)
df_genes <- df_genes[, c(ncol(df_genes), 1:(ncol(df_genes)-1))]  # Colocar la columna Gene al inicio


#La expresión génica de los 20 genes con mayor varianza (genes_top_limpios)
write.csv(df_genes, "resultados/genes_top_20.csv", row.names = FALSE)

#La anotación de los subtipos de cada muestra (paraheat), que sirve para agrupar o visualizar en heatmaps.
write.csv(paraheat, "resultados/subtipos_muestras.csv", row.names = TRUE)

#visualizarlo de forma numerica los subtipos
# Limpiar subtipos (quitar NA) y crear factor
muestras_buenas <- !is.na(subtipos)
subtipos_limpios <- subtipos[muestras_buenas]
subtipos_factor <- as.factor(subtipos_limpios)

# Contar cuántas muestras hay por subtipo 
tabla_subtipos <- table(subtipos_factor)
tabla_subtipos
# Convertir a data.frame,para explicar
df_tabla_subtipos <- as.data.frame(tabla_subtipos)
colnames(df_tabla_subtipos) <- c("Subtipo", "Numero_muestras")

write.csv(df_tabla_subtipos, file = "resultados/conteo_subtipos.csv", row.names = FALSE)
print(df_tabla_subtipos)

