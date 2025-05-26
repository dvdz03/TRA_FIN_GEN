##############################
#### ESTE ES PARA LA RED #####
##############################

miRNAs <- get_multimir(org = "hsa", target = "TP53", table = "predicted")
miRNABONITO <- miRNAs@data %>% filter(score > 0.8)
plot(graph_from_data_frame(miRNABONITO), layout = layout_with_fr)
head(miRNABONITO)
#AHÍ ESTÁ LA RED, SOLO FALTA EDITARLA PARA QUE SE VEA BIEN Y TODO ESO 

miRNAs1 <- get_multimir (org = "hsa", target = genes_ya$ENTREZID[1:10], table = "predicted" )
genes_ya$ENTREZID[1:20]
miRNABONITO1 <- miRNAs1@data %>% filter(score > 0.8)
plot(graph_from_data_frame(miRNABONITO1), layout = layout_nicely, font.size = 0.1, )
plot(graph_from_data_frame(miRNABONITO1),
     layout = layout_nicely,
     vertex.size = 10,
     vertex.color = "lightblue",
     vertex.label.cex = 0.3,
     vertex.label.dist = 1.5,
     edge.width = 0.5,
     edge.arrow.size = 0.2,
     main = "Red")
