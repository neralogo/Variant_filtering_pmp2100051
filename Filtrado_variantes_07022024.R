#Test filtrado variantes con un solo trío. Cargar los tres arhivos de padre
#madre e hijo.

library(dplyr)
install.packages("devtools")
devtools::install_github("PhanstielLab/bedtoolsr")

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

new_column_names <- c("ENSEMBLE_ID", "Start", "End", "Chr", "GENE_NAME")

colnames(genes_no_dup)[colnames(genes_no_dup) %in% c("Gene.stable.ID", 
       "Gene.start..bp.", "Gene.end..bp.", "Chromosome.scaffold.name", 
       "Gene.name")] <- new_column_names


GENES_FILTRADO <- genes_no_dup %>%
  select(Chr, Start, End, GENE_NAME, ENSEMBLE_ID)

GENES_FILTRADO$Start <- GENES_FILTRADO$Start - 1
GENES_FILTRADO$Chr <- paste0("chr", GENES_FILTRADO$Chr)


#Guardar archivo con formato .bed del conjunto de genes a filtrar con bedtools
write.table(GENES_FILTRADO, "GENES_FILTRADO.bed", sep = "\t", quote = F, 
            col.names = F, row.names = F)

G01.GEA.10.HI.split.tab$Start <- G01.GEA.10.HI.split.tab$Start -1 
write.table(G01.GEA.10.HI.split.tab, "Prueba_pac10.bed", sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)

#
system("bedtools intersect -a Prueba_pac10.bed -b GENES_FILTRADO.bed > intersect.bed")


