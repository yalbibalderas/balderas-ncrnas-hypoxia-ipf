############################################################################################################
# 1. Set working directory
# setwd("path-to-files")
setwd("/Volumes/GoogleDrive/My Drive/CEL_ARNOLDO_INER")
############################################################################################################

############################################################################################################
# 2. Load necessary packages and functions
############################################################################################################
# To Read and process data files
library(affy)
library(dplyr)
library(tidyr)
library(reshape2)
library(limma)
library(oligo)
library(pd.clariom.d.human)
library(affycoretools)
library(M3C)
library(rstatix)

# To Produce plots
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(ggrepel) # for nice annotations
library(gridExtra)
library(ComplexHeatmap)
library(RColorBrewer)
library(corrplot)
library(gridGraphics)
library(grid)
library(gridExtra)

# Function used to compute z-score values
scale_rows <-function (x)
{
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m)/s)
}

############################################################################################################
# 3. Load arrays, perform quality control and produce normalized expression matrix 

# 3.1 Read and ordering CEL files
celfiles <- list.celfiles(path=getwd(),full.names=TRUE)
readraw=read.celfiles(celfiles)
raw = readraw

# 3.2 Preprocessing data: Load; Bckg_correction; Normalization; Summarization
eset.rma= rma(raw,normalize=TRUE, background=TRUE) 
eset_full= annotateEset(eset.rma, pd.clariom.d.human)

# 3.3 Select most variable genes, those above the 50 percentile.
filtered_results <-featurefilter(exprs(eset_full), percentile = 50, method='MAD', topN=0)
eset_full_filt <- eset_full[rownames(filtered_results$filtered_data), ]

# 3.4 Create Experimental design table
array <- gsub(".CEL", "", sampleNames(readraw))
# Reorder samples for exp_design, according to last sample ID from Yahir
new_order <- order(factor(sub("_.*", "", array), levels = c("F1", "F2", "F3", "C1", "C2", "C3")))
array <- array[new_order]
# Define group tag vector
group <- rep(c("IPF", "Ctr"), each = 18)
# Define replicates tag vector
replicate <- rep(c(1,2,3), times = 12)
# Define treatment tag vector
treatment <- rep(c("HOx", "nOx"), each = 3, times =6)
# Define sample_name tag vector
sample_name <- paste(group, treatment, replicate, sep = "_")
# Create Data frame of experimental design 
exp_design <- data.frame(array, group, replicate, treatment, sample_name)
# Write out exp_design table as csv
write.csv(exp_design, file="Experimental_design.csv")

# 3.5 Reorder Samples in normalized expression data according to their order in exp_design)
eset_full_filt <- eset_full_filt[, new_order]
sampleNames(eset_full_filt) <- exp_design$array


################################################################################
###Differential expression analysis	
# input: Pre-processed expression Matrix and experimental design data frame
################################################################################
# 4.1 Write Function to perform Differential expression Analysis
difExp.limma<- function(mat,conditions,Levels,indxcontr=c(1:length(contrasts)),averageCont,P,FC, separate.contrasts=FALSE) {
  #Design Matrix 
  Levels= Levels
  
  f <- factor(conditions, levels=Levels)
  design <- model.matrix(~0+f)
  colnames(design)= make.names(Levels)
  fit <- lmFit(mat, design)
  
  P<- P #select threshold for P.value
  FC<- FC
  M <- log(FC,2) #select threshold for logFC
  
  #Make contrasts vector
  Leveldop<- Levels
  Leveldop<- make.names(Leveldop)
  n=length(Leveldop)
  contrasts <-c()
  for (i in 1:(n-1))
  {
    cont_1 <- Leveldop[i]
    contrast <- c()
    for (j in (i+1):n)
    {
      contrast <-c(contrast,(paste(Leveldop[j],cont_1,sep="-")))
    }
    contrasts <-c(contrasts,contrast)
    
  }
  contrasts <- contrasts[indxcontr]
  # Make contrast matrix
  cont.matrix= makeContrasts( contrasts=contrasts ,levels=design)
  if(averageCont==TRUE){cont.matrix =rowMeans(cont.matrix)}else{cont.matrix=cont.matrix}
  
  #compute lineal model for contrasts
  fit2 = contrasts.fit(fit, cont.matrix)
  eb = eBayes(fit2)
  
  # Select t.test output by contrast or F-value (ANNOVA style)
  if(separate.contrasts==TRUE)
  {
    DiffExp_Table = NULL
    for (i in 1:length(contrasts))
    {
      TTFch= topTable(eb, coef =i, number=nrow(mat))
      TTFch$Contrasts = contrasts[i]
      TTFch$Gene = rownames(TTFch)
      DiffExp_Table = rbind(DiffExp_Table,TTFch)
      
    }
    return(DiffExp_Table)
    
  }else{
    TTFch= topTable(eb,number=nrow(mat))
    dtestsig=decideTests(eb, method="separate",p.value=P, lfc=M, adjust.method="fdr")
    dtestsig_all_times=data.frame(dtestsig[apply(dtestsig,1,function(row) any(row)!=0),])
    outlist=list(Contrasts=contrasts,ebResults=eb, dtestsig=dtestsig_all_times, topTab=TTFch )
    return(outlist)	
  }	
}

# 4.2 Perform Differential expression Analysis
topTabTableGene=list()
eset4DExpList = list()
for (gr in 1:length(unique(group)))
{
  
  ### Create expression subset by treatment
  lables_genes_sub=exp_design[grep(unique(group)[gr],exp_design$group),]
  eset = eset_full_filt[,grep(unique(group)[gr],exp_design$group)]

  ### Select ncRNA features from Annotation data and remove NAs 
  annota= na.omit(eset@featureData@data)
  annota = annota[grep("NR", annota$ID),]
  
  ### Select ncRNA features in expression subset, that matches the annotation
  esetAnn=eset[featureNames(eset) %in% annota$PROBEID,]
  ### Collapse expression values by taking the mean across probes assignet to the same gene 
  eset4DExp=avereps(esetAnn, ID=annota[,"SYMBOL"])
  
  ################################################################################
  ## This code invokes the difExp.limma function, 
  ## which executes the DE computations for the selected dataset and,
  ## returns the resulting tables
  Dex.tab=difExp.limma(
    mat = eset4DExp, # Normalized expression matrix for all included samples
    conditions = lables_genes_sub$treatment, # Experimental groups to be compared
    Levels = rev(unique(lables_genes_sub$treatment)), # Unique Experimental groups to be compared
    indxcontr=c(1),
    averageCont=FALSE, 
    P=1, # p.value cutoff
    FC=0, # log2FCh cutoff
    separate.contrasts=TRUE) # Boolean to select output as table for each paired contrast OR annova style)
  
  ### Write out matrix of normalized expression values of the selected group of samples
  eset4DExpList[[gr]] <- eset4DExp
  filenameGenes=paste(unique(group)[gr],"_NormExp.txt",sep="")
  write.table(eset4DExp, filenameGenes, quote=FALSE,sep="\t")
  
  ### List of Genes F.p values 
  topTabTableGene[[gr]] = Dex.tab
  
  ### Write out Table of DEA results for each contrast
  tempListName= as.character(unique(group)[gr]) # define prefix for the output file name
  write.csv(topTabTableGene[[gr]], file = paste(tempListName,".table_DE_genes_NR.csv",sep=""))
}

########################################################################################################
###Figure A: heatmap: Normalized expression of differentially expressed ncRNAs

#Reading DE tables
IPF_fibroblasts = read.csv("IPF.table_DE_genes_NR.csv", row.names = 1, check.names = FALSE)
Ctr_fibroblasts = read.csv("Ctr.table_DE_genes_NR.csv", row.names = 1, check.names = FALSE)

## Reading Normalized expression matrix
IPF_expmat = read.delim("IPF_NormExp.txt")
Ctr_expmat = read.delim("Ctr_NormExp.txt")

# Prepare matrix for heatmap
ht_list = list()
DEA_list = list(IPF_fibroblasts, Ctr_fibroblasts)
mat_list = list(IPF_expmat, Ctr_expmat)
widths = c(6, 10)
heights = c(10, 25)
fontsizes = c(13, 7)
for(k in 1:length(DEA_list)){
  
  DEA_table = DEA_list[[k]]
  DEA_table_Up = DEA_table[DEA_table$adj.P.Val <0.05 & DEA_table$logFC > 0.7,  ]
  DEA_table_Down = DEA_table[DEA_table$adj.P.Val <0.05 & DEA_table$logFC < -0.7,  ]
  
  mat = mat_list[[k]][ c(DEA_table_Up$Gene, DEA_table_Down$Gene), ]
  mat_scale = scale_rows(mat)
  
  # Create vector of samples labels as input to be splitted by group in heatmap
  colsplit = exp_design[exp_design$group == unique(exp_design$group)[k], ]$treatment
  colsplit = factor(colsplit, levels= c("HOx","nOx"))
  
  # Create vector of samples labels as input to be splitted by group in heatmap
  rowsplit = c(rep("Up-regulated", times = length(DEA_table_Up$Gene)), 
               rep("Down-regulated", times = length(DEA_table_Down$Gene)))
  rowsplit = factor(rowsplit, levels= c("Up-regulated","Down-regulated"))
  
  # Draw heatmap
  ht_list[[k]] = Heatmap(mat_scale, name = "Z-score",
               row_split = rowsplit,
               row_names_side = "left",
               row_names_gp = gpar(fontsize = fontsizes[k]),
               row_title_gp = gpar(fill = c("red4", "steelblue4"), col = "white" , font=2, fontsize = 13),
               column_split = colsplit,
               column_labels = colnames(mat_scale),
               column_names_gp = gpar(col = c("black"), fontsize = 10),
               column_title_gp = gpar(fill = c("orange3", "gray70"), font=2, fontsize = 16),
               show_row_dend = FALSE,
               show_column_dend = FALSE,
               cluster_columns = FALSE,
               cluster_row_slices = FALSE,
               width = unit(widths[k], "cm"), 
               height = unit(heights[k], "cm"))
  ### Plot single Heatmap
  # png(file="Hmap_DE_IPF_fibroblasts_NR_NEW.png", units="in", width=7, height=7, res=300)	
  # print(ht)
  # dev.off()
}


grab_grob <- function(){
  grid.grab()
}

drawGridHeatmap  <- function(hm) {
  draw(hm)
  grab_grob()
}

# Use rev function to invert the order of heatmaps on the canvas.
gl <- lapply(rev(ht_list), drawGridHeatmap) 

png(file= "Hmap_DE_NR_NEW.png", units="in", width=20, height=12, res=400)
grid.newpage()
grid.arrange(grobs=gl, ncol=2, clip=TRUE)
dev.off()

########################################################################################################
###Figure B: Volcano Plot
# vplot_list = list()
DEA_list = list(IPF_fibroblasts, Ctr_fibroblasts)
plot_names <- c("IPF-fibroblasts", "Control-fibroblasts")
widths = c(5, 10)
heights = c(10, 25)
fontsizes = c(13, 7)

for(k in 1:2)
{
  
  volcano_df <- DEA_list[[k]]
  
  # Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
  volcano_df <- volcano_df %>%
    mutate(expression_status = case_when(
      logFC > 0.7 & adj.P.Val < 0.05 ~ "UP",  # Change threshold as needed
      logFC < -0.7 & adj.P.Val < 0.05 ~ "DOWN",
      TRUE ~ "NotSig"  # For logFC between -1 and 1
    ))
  
  # Create a new column "delabel" to de, that will contain the name of the top 30 differentially expressed genes (NA in case they are not)
  volcano_df$delabel <- ifelse(volcano_df$Gene %in% head(volcano_df[order(volcano_df$adj.P.Val) & volcano_df$expression_status != "NotSig", "Gene"], 20) | volcano_df$Gene =="MEG3", volcano_df$Gene, NA)
  
  xlim <- max (max(abs(volcano_df$logFC)), 5)
  ymax <-  max(-log10(volcano_df$adj.P.Val)) + 2
  
  # Plot the volcano 
  vplot <- ggplot(data = volcano_df, aes(x = logFC, y = -log10(adj.P.Val), col = expression_status, label = delabel)) + 
    geom_vline(xintercept = c(-0.7, 0.7), col = "grey55", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "grey55", linetype = 'dashed') +
    geom_point(size = 2) +
    theme_minimal() +
  
    scale_color_manual(values = c("#153f65", "grey55", "#bb0c00"), # to set the colours of our variable
                       labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
    coord_cartesian(ylim = c(0, ymax), xlim = c(-xlim, xlim)) + # since some genes can have minuslog10padj of inf, we set these limits
    labs(color = 'Severe', #legend_title
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
    scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
    ggtitle(paste0(plot_names[k],' under low oxygen vs normal oxygen levels')) + # Plot title
    
    theme(plot.title = element_text(lineheight=.5, face="bold", size=17, hjust=0),
          axis.text.x = element_text(face="bold", color="black",size=17),
          axis.title.x = element_text(color="black", size=20, face="bold"),
          axis.text.y = element_text(face="bold", color="black", size=17),
          axis.title.y = element_text(color="black", size=20, face="bold"),
          legend.position="none",
          panel.border = element_rect(colour = "black", fill=NA, size=2))+
    
    geom_text_repel(max.overlaps = 21, nudge_x = -1,  nudge_y = -0.5) # To show all labels 
  
  png(file= paste0(plot_names[k], "_VolcanoPlot_DE_NR_NEW.png"), units="in", width=8, height=8, res=400)
  print(vplot)
  dev.off()
}