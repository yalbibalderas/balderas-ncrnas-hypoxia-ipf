
# Import file from rise
rise <- read.table("Network/rise_human_all.txt", sep = "\t", header=F)
columns <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "rise_id", "score", "strand1", 
              "strand2", "gene_id1", "gene_name1", "gene_id2", "gene_name2", "gene_type1", "gene_type2",
              "genomic_context1", "genomic_context2", "method", "species", "cell_line", "pubmed_id")
colnames(rise) <- columns
head(rise)
dups <- rise$gene_name2[duplicated(rise$gene_name2) %in% rise$gene_name1]


# Import ncRNAs 
fibrotic <- read.table("Network/DEG_Fibrotic_Hipoxia_vs_Fibrotic_Normoxia_logFC_0.7_NR_NCBI_non-coding_RNA.txt", sep = "\t", header=T) 
normal <- read.table("Network/DEG_Normal_Hipoxia_vs_Normal_Normoxia_logFC_0.7_NR_NCBI_non-coding_RNA.txt", sep = "\t", header=T)
colnames(fibrotic)

# Fibrotic
merged_df1 <- merge(rise, fibrotic, by.x = "gene_name1", by.y = "SYMBOL")
merged_df2 <- merge(rise, fibrotic, by.x = "gene_name2", by.y = "SYMBOL")
merged_df <- rbind(merged_df1, merged_df2)
dup_rows <- duplicated(merged_df)
df_unique <- subset(merged_df, !dup_rows)
write.csv(df_unique, file="Network/fibrotic_rise.csv")

# Normal
merged_df1 <- merge(rise, normal, by.x = "gene_name1", by.y = "SYMBOL")
merged_df2 <- merge(rise, normal, by.x = "gene_name2", by.y = "SYMBOL")
merged_df <- rbind(merged_df1, merged_df2)
dup_rows <- duplicated(merged_df)
df_unique <- subset(merged_df, !dup_rows)
write.csv(df_unique, file="Network/normal_rise.csv")







