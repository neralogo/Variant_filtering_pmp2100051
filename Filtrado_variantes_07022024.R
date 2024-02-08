#Test filtrado variantes con un solo trío. Cargar los tres arhivos de padre
#madre e hijo.

library(dplyr)

#Después de cargar la base de Genes de TEA y la de biomart, eliminamos la única
#instancia de gen que no tiene nombre. Tampoco consta en la base de biomart por
#su ensemble_id, por eso la elimino.
GENES_DIAG_TEA <- na.omit(GENES_DIAG_TEA)

#Filtramos los genes del listado de mart_export que coincidan con el nombre
#del conjunto de genes de TEA.
filtered_mart_export <- mart_export %>%
  filter(Gene.name %in% GENES_DIAG_TEA$GENE_NAME)

#Creamos una BBDD con los genes que faltan entre los filtrados y el conjunto de 
#genes TEA inicial. Son 15 los que no aparecen en mart_export por su nombre.
genes_faltantes <- anti_join(GENES_DIAG_TEA, filtered_mart_export, by = c("GENE_NAME" = "Gene.name"))
#Nos aseguramos que estén todos en mayúscula para que los encuentre en la columna
#de sinónimos
genes_faltantes$GENE_NAME[2] <- toupper(genes_faltantes$GENE_NAME[2])

#Generamos un dataset temporal para guardar las filas de estos 15 genes que encontramos
#en sinónimos.
temp <- mart_export %>%
  filter(Gene.Synonym %in% genes_faltantes$GENE_NAME)

#Juntamos ambos datasets
genes_intersect <- bind_rows(filtered_mart_export, temp)

#Eliminamos el dataset temporal
rm(temp)

#Quitamos la columna de nombres de genes sinónimos
genes_intersect$Gene.Synonym <- NULL

#Quitamos los duplicados de los genes incluidos por los sinónimos.
genes_no_dup <- genes_intersect %>% distinct(Gene.name, .keep_all = TRUE)

#
system("bedtools intersect -a Exoma_HGUGM.bed -b GENE_BBDD.bed > intersect.bed")


