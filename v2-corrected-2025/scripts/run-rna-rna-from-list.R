# Cargar la función rna_rna desde otro archivo si está guardada por separado
# source("rna_rna_function.R")

# Librerías necesarias
library(readr)
source("encori2.R")

# Leer lista de ncRNAs desde archivo CSV o TXT
# Asegúrate de que tu archivo tenga una columna llamada 'RNA'
ncrna_file <- "Results_Ivan/IPF_table_Differential_expression_analysis_NR.csv"  # Cambia esto si usas .txt o si tiene otro nombre
ncrna_file <- "Results_Ivan/Control_table_Differential_expression_analysis_NR.csv"  # Cambia esto si usas .txt o si tiene otro nombre

ncRNAs <- read_csv(ncrna_file, show_col_types = FALSE)
ncRNAs2 <- ncRNAs[ncRNAs$adj.P.Val < 0.05 & abs(ncRNAs$logFC) > 0.7, ]

# Ejecutar consulta
resultadosIPF <- rna_rna(geneType = "lncRNA", RNA = ncRNAs2$RNA)
resultadosControl <- rna_rna(geneType = "lncRNA", RNA = ncRNAs2$RNA)

# Guardar resultados
if (!is.null(resultados)) {
  write_csv(resultados, "rna_interactions_resultadosIPF2.csv")
  cat("✅ Resultados guardados en 'rna_interactions_resultadosIPF2.csv'\n")
} else {
  cat("⚠️ No se obtuvieron resultados para los ncRNAs proporcionados.\n")
}

# Guardar resultados
if (!is.null(resultados)) {
  write_csv(resultados, "rna_interactions_resultadosControl.csv")
  cat("✅ Resultados guardados en 'rna_interactions_resultadoscControl.csv'\n")
} else {
  cat("⚠️ No se obtuvieron resultados para los ncRNAs proporcionados.\n")
}

# Para el caso de los microRNAs hay que hacer una conversión de formato primero

# Crear nueva columna con formato ENCORI IPF
ncRNAs2$miRNA_encori <- tolower(ncRNAs2$RNA)
ncRNAs2$miRNA_encori <- ifelse(grepl("^mir", ncRNAs2$miRNA_encori),
                               gsub("^mir", "hsa-miR-", ncRNAs2$miRNA_encori),
                               ifelse(grepl("^let", ncRNAs2$miRNA_encori),
                                      gsub("^let", "hsa-let-", ncRNAs2$miRNA_encori),
                                      ncRNAs2$miRNA_encori))
ncRNAs2_miRNAs_subset <- ncRNAs2[grepl("^hsa", ncRNAs2$miRNA_encori), ]
resultadosIPFmiRNAs <- rna_rna(geneType = "miRNA", RNA = ncRNAs2_miRNAs_subset$miRNA_encori) # No hay
resultadosIPFmiRNAs <- rna_rna(geneType = "miRNA", RNA = "hsa-miR-1244")

if (!is.null(resultadosIPFmiRNAs)) {
  write_csv(resultadosIPFmiRNAs, "rna_interactions_resultadosIPFmiRNAs.csv")
  cat("✅ Resultados guardados en 'rna_interactions_resultadosIPF2miRNAs.csv'\n")
} else {
  cat("⚠️ No se obtuvieron resultados para los ncRNAs proporcionados.\n")
}

# Crear nueva columna con formato ENCORI Control
ncRNAs2$miRNA_encori <- tolower(ncRNAs2$RNA)
ncRNAs2$miRNA_encori <- ifelse(grepl("^mir", ncRNAs2$miRNA_encori),
                               gsub("^mir", "hsa-miR-", ncRNAs2$miRNA_encori),
                               ifelse(grepl("^let", ncRNAs2$miRNA_encori),
                                      gsub("^let", "hsa-let-", ncRNAs2$miRNA_encori),
                                      ncRNAs2$miRNA_encori))
ncRNAs2_sub <- ncRNAs2[grepl("^hsa", ncRNAs2$miRNA_encori), ]
ncRNAs2_sub <- ncRNAs2_sub |> 
  dplyr::mutate(
    RNA_3p = paste0(miRNA_encori, "-3p"),
    RNA_5p = paste0(miRNA_encori, "-5p")
  )
ncRNAs_expanded <- c(ncRNAs2_sub$miRNA_encori, ncRNAs2_sub$RNA_3p, ncRNAs2_sub$RNA_5p)


resultadosControlmiRNAs <- rna_rna(geneType = "miRNA", RNA = ncRNAs_expanded) 
resultadosControlmiRNAs <- rna_rna(geneType = "miRNA", RNA = c("hsa-miR-21-5p")) 

