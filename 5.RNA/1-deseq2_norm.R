setwd("Path_to/FINAL_RNA/")
library(DESeq2)
library(ggplot2)

#Samples:
# "T18__P10610_116"
# "T19__P10610_117"
# "N18_Proximal__P10610_129"
# "N19_Proximal__P10610_130"
# "N19_Distal__P10610_148" 

rm(list=ls())
# Input 2 files:
# [1] Raw count matrix   (ensemblid, sample1, sample2, .. sampleN,)
# [2] groups matrix      (samples, sex, treatment1, treatment2 ... treatmentn) (only one column after samples is required)
# note: samples names from [1] should be identical to samples column in [2]
count_matrix <- read.table("../merged_gene_counts.csv",header = T, sep = "\t")
groups  <- read.table("../C.Williams_18_02_sample_info.csv", header = T, sep = "\t")

#keep wanted samples, in this case it will affect the normalization
wanted <- groups[grepl(pattern = "18", groups$User_ID, fixed = T) |
                   grepl(pattern = "19", groups$User_ID, fixed = T),]

#remove the other normal
wanted <- wanted[wanted$NGI_ID != "P10610_148",]

counts <- count_matrix[c("ENSEMBL_ID",wanted$NGI_ID)]

groups <- wanted
#assign row names
groups <- data.frame(groups, row.names = 'NGI_ID')
counts <- data.frame(counts, row.names = 'ENSEMBL_ID')

#remove low counts genes
counts <- counts[rowSums(counts) >= 10,]

#assign groups
groups$Status <- substr(groups$User_ID,1,1)

#put WT to the left, to be the reference
groups$Status <- factor(groups$Status, levels = c('N', 'T'))


matrix <- counts
head(matrix)
colData <- groups
groups

#row.names(colData) <- colnames(matrix)
dds <- DESeqDataSetFromMatrix(countData = matrix,
                              colData = colData,
                              design= ~ Status )  #which column of to be used in the comparison genotype..group..sex

dds <- estimateSizeFactors(dds)
dds <- DESeq(dds,)

#getting the results
resultsNames(dds)

normCounts <- counts(dds, normalized=TRUE)
res <- results(dds)

res <- res[rownames(normCounts), ]

Final_matrix <- cbind(normCounts,res)
resOrdered <- data.frame(Final_matrix)

rm("dds","colData","matrix","res", "counts",
   "normCounts", "count_matrix", "groups", "wanted","Final_matrix")


# --------------------------- Genes' names annotation  --------------------------

library(rtracklayer)

build_gencode_lookup <- function(gtf_path) {
  gtf <- rtracklayer::import(gtf_path, format = "gtf")
  df <- as.data.frame(gtf)
  rm(gtf)
  # Keep only gene-level entries to avoid duplicates
  df <- df[df$type == "gene", 
           c("gene_id", "gene_name", "gene_type", 
             "seqnames", "strand", "start", "end")]
  # Strip version from gene_id to match typical DESeq2 rownames
  df$gene_id_stripped <- gsub("\\.[0-9]+$", "", df$gene_id)
  return(df)
}

# Run ones and save the table for future use, as it takes time to run
# ---------------- 
ref_annotaion <- build_gencode_lookup("Path_to/annotation/gencode.v49.annotation.gtf")
write.table(ref_annotaion, file = "gencode_v49_lookup.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
# ---------------- 

ref_annotaion <- read.table("gencode_v49_lookup.tsv", header = TRUE, sep = "\t")
  
add_gene_info_gtf <- function(resOrdered, lookup) {
  labling <- gsub("\\.[0-9]+$", "", rownames(resOrdered))
  
  idx <- match(labling, lookup$gene_id_stripped)
  
  n_missing <- sum(is.na(idx))
  if (n_missing > 0) warning(n_missing, " gene IDs had no match in GTF.")
  
  annotation <- lookup[idx, c("gene_name", "gene_type", 
                              "seqnames", "strand", "start", "end")]
  
  resOrdered <- cbind(ensembl_id = rownames(resOrdered), 
                      annotation, resOrdered)
 
  return(resOrdered)
}


Final_matrix2 <- add_gene_info_gtf(resOrdered, ref_annotaion)

head(Final_matrix2,3)

dim(Final_matrix2)



#update ID with samples names
for (el in 1:nrow(wanted))
{
 if (wanted[el,1] %in% colnames(Final_matrix2))
 {
   colnames(Final_matrix2)[colnames(Final_matrix2)==wanted[el,1]] <- wanted[el,2]
     # paste0(wanted[el,2],"__",wanted[el,1])
 }
}


#remove NA as the previous reference was HG37
dim(Final_matrix2[is.na(Final_matrix2$ensembl_gene_id),])

Final_matrix2 <- Final_matrix2[!is.na(Final_matrix2$ensembl_gene_id),]
              
#SAVE
write.csv(file = 'Deseq2_log2fold-NvsT_geneloci_hg38_filtered_normalized.csv', Final_matrix2) #norm2 after removing distal n P10610_148 and removing deprecated genes

#CLEAN
#rm(list = ls())


