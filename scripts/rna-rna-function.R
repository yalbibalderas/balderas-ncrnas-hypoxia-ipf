rna_rna <- function(assembly="hg38",
                    geneType,
                    RNA,
                    interNum=1,
                    expNum=1,
                    cellType="all"){
  rna_interactions <- NULL
  for (gen in geneType){
    for (rna in RNA){
      for (cell in cellType){
        
        # Dentro del bucle, justo antes de usar el link:
        encoded_rna <- URLencode(rna, reserved = TRUE)
        
        link <- paste("https://rnasysu.com/encori/api/RNARNA/?",
                      "assembly=",assembly,
                      "&geneType=",gen,
                      "&RNA=",rna,
                      "&interNum=",interNum,
                      "&expNum=",expNum,
                      "&cellType=",cell, sep = "")

        rna_int <- utils::read.csv(url(link), comment.char = "#", sep = "\t", row.names = NULL)

        if(rna_int[1,1] != "The RNA parameter haven't been set correctly! Or the input of RNA parameter is not available!"
           & !is.na(rna_int[1,1])){
          rna_interactions <- rbind(rna_interactions, rna_int)
        }

      }
    }
  }
  if (is.null(rna_interactions)){
    print("Data not available or incorrect parameters")
  } else{
    BiocGenerics::colnames(rna_interactions)[1] <- "geneID"
    return(rna_interactions)
  }
}
