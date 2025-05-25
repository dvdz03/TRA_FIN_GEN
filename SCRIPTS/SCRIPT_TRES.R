##############################
#### ESTE ES PARA LA RED #####
##############################

miRNAs <- get_multimir(org = "hsa", target = "TP53", table = "predicted")
miRNABONITO <- miRNAs@data %>% filter(score > 0.8)
plot(graph_from_data_frame(miRNABONITO), layout = layout_with_fr)
head(miRNABONITO)
#AHÍ ESTÁ LA RED, SOLO FALTA EDITARLA PARA QUE SE VEA BIEN Y TODO ESO 