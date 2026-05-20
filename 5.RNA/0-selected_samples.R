setwd("/Path_to/FINAL_RNA/")

count_matrix <- read.table("../merged_gene_counts.csv",header = T, quote = ",")
sample_info  <- read.table("../C.Williams_18_02_sample_info.csv", header = T, quote = "\t")

wanted <- sample_info[grepl(pattern = "18", sample_info$User_ID, fixed = T) |
            grepl(pattern = "19", sample_info$User_ID, fixed = T),]
wanted_matrix <- count_matrix[wanted$NGI_ID]
rownames(wanted_matrix) <- count_matrix$ENSEMBL_ID

rownames(wanted) <- wanted$NGI_ID

wanted <- t(wanted)

wanted_matrix <- rbind(wanted, wanted_matrix)

#update IDs
colnames(wanted_matrix) <- paste0(wanted_matrix["User_ID",],"__",colnames(wanted_matrix))

write.csv(wanted_matrix, "wanted.samples.RNA.csv")

rm(list = ls())
