################################
#### ESTE ES PARA ONTOLOGÍA ####
################################

genes_significativos <- resultados %>% 
  filter(abs(logFC)>1 & adj.P.Val > 0.05)%>%
  rownames ()
#para saber que tipo de genes tengo
head(genes_significativos)
keytypes(org.Hs.eg.db)# esto es para ver las opciones que puedo hacer

genes_buenos <- sub("\\..*","", rownames(resultados))
genes_ya <- bitr(genes_buenos,
                 fromType = "ENSEMBL",
                 toType = "ENTREZID",
                 OrgDb = "org.Hs.eg.db")
head(genes_ya)
cosago <- enrichGO( gene = genes_ya$ENTREZID,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    readable = TRUE)
head(cosago)
funbio <- barplot(cosago,
        showCategory = 20, 
        title= "ontología (Proceso biológico)",
        font.size = 7)
funbio

molecular <- enrichGO( gene = genes_ya$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       readable = TRUE)
mole <- barplot(molecular,
                showCategory = 20,
                title = "ontología (Función Molecular)",
                font.size = 7)
mole

elotro <- enrichGO(gene = genes_ya$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   readable = TRUE)
otro <- barplot(elotro, 
                showCategory = 20,
                title = "ontología (COMPONENTE CELULAR)",
                font.size = 7)

otro

