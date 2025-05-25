#################################
#### ESTE ES PARA EL HEATMAP ####
#################################

paraheat <- data.frame(Subtype = subtipos_factor)
rownames(paraheat) <- colnames(genes_top_limpios)
pheatmap(mat = genes_top_limpios,
         annotation_col = paraheat,
         show_colnames = FALSE,
         color = colorRampPalette(c("orange","black","darkorchid"))(100))
