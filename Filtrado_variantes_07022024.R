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

new_column_names <- c("ENSEMBLE_ID", "Start", "End", "Chr", "GENE_NAME")

colnames(genes_no_dup)[colnames(genes_no_dup) %in% c("Gene.stable.ID", 
       "Gene.start..bp.", "Gene.end..bp.", "Chromosome.scaffold.name", 
       "Gene.name")] <- new_column_names


GENES_FILTRADO <- genes_no_dup %>%
  select(Chr, Start, End, GENE_NAME, ENSEMBLE_ID)

GENES_FILTRADO$Chr <- paste0("chr", GENES_FILTRADO$Chr)


#Guardar archivo con formato .bed del conjunto de genes a filtrar con bedtools
write.table(GENES_FILTRADO, "GENES_FILTRADO.bed", sep = "\t", quote = F, 
            col.names = F, row.names = F)

G01.GEA.10.HI.split.tab$Start <- G01.GEA.10.HI.split.tab$Start -1 
G01.GEA.10.HI.split.tab$Chr <- paste0("chr", G01.GEA.10.HI.split.tab$Chr)
write.table(G01.GEA.10.HI.split.tab, "GEA10HI.bed", sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)

G01.GEA.10.MA.split.tab$Start <- G01.GEA.10.MA.split.tab$Start -1 
G01.GEA.10.MA.split.tab$Chr <- paste0("chr", G01.GEA.10.MA.split.tab$Chr)
write.table(G01.GEA.10.MA.split.tab, "GEA10MA.bed", sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)

G01.GEA.10.PA.split.tab$Start <- G01.GEA.10.PA.split.tab$Start -1 
G01.GEA.10.PA.split.tab$Chr <- paste0("chr", G01.GEA.10.PA.split.tab$Chr)
write.table(G01.GEA.10.PA.split.tab, "GEA10PA.bed", sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)

rm(G01.GEA.10.HI.split.tab)
rm(G01.GEA.10.MA.split.tab)
rm(G01.GEA.10.PA.split.tab)

system("bedtools intersect -a GEA10_HI.bed -b GENES_FILTRADO.bed > GEA10_HI_intersect.bed")
system("bedtools intersect -a GEA10_MA.bed -b GENES_FILTRADO.bed > GEA10_MA_intersect.bed")
system("bedtools intersect -a GEA10_PA.bed -b GENES_FILTRADO.bed > GEA10_PA_intersect.bed")


temporal <- colnames(G01.GEA.10.HI.split.tab)
colnames(GEA10_HI_intersect) <- temporal
colnames(GEA10_MA_intersect) <- temporal
colnames(GEA10_PA_intersect) <- temporal
rm(temporal)
rm(genes_faltantes)
rm(genes_intersect)
rm(filtered_mart_export)
rm(GENES_DIAG_TEA)
rm(mart_export)
rm(mart_export_csv)
rm(genes_no_dup)

# Function to categorize variants based on conditions
categorize_variant <- function(row, mother_data, father_data) {
  variant_key <- paste0(row[1], "_", row[2])  # Assuming Start is in the first column and End is in the second
  
  is_de_novo <- is.na(match(variant_key, paste0(mother_data$Start, "_", mother_data$End))) &&
    is.na(match(variant_key, paste0(father_data$Start, "_", father_data$End)))
  
  is_heredada_MA <- !is_de_novo && !is.na(match(variant_key, paste0(mother_data$Start, "_", mother_data$End)))
  is_heredada_PA <- !is_de_novo && !is.na(match(variant_key, paste0(father_data$Start, "_", father_data$End)))
  
  if (is_de_novo) {
    return("de novo")
  } else if (is_heredada_MA && is_heredada_PA) {
    return("heredada_MA_PA")
  } else if (is_heredada_MA) {
    return("heredada_MA")
  } else if (is_heredada_PA) {
    return("heredada_PA")
  }
}

# Assuming "Start" and "End" are in the first and second columns of patient_data
# Apply the categorize_variant function to create a new column in patient_data
GEA10_HI_intersect$Variant_Category <- apply(GEA10_HI_intersect[, c("Start", "End")], 1, 
                                             categorize_variant, GEA10_MA_intersect, GEA10_PA_intersect)

# Print the updated patient_data
print(GEA10_HI_intersect$Variant_Category)

# Move the "Variant_Category" column to the 7th position
GEA10_HI_intersect <- GEA10_HI_intersect %>%
  select(1:6, Variant_Category, everything())

