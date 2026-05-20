# The analysis has 3 modalities: RNA expression (R), CNA (D), and mutations (M).
# The goal is to integrate them together to see if there are any correlations between them.
# For example, if a gene is amplified in the CNA data, is it also overexpressed in the RNA data?
# Or if a gene is mutated, does it show any changes in expression or copy number?


## plot with CopyNumberPLots
#############################
setwd("/Path_to/FINAL_RNA/")

library(data.table)
library(GenomicRanges)
library(karyoploteR)
library(CopyNumberPlots)
library(tidyverse)
library(tidyr)

rm(list = ls())
#######################################################
##        1. Integrate gene expression
##       With small mutations (Cosmic_dup_withHP)
#######################################################


#Load the normalized RNA reads
rna_exp_matrix <- read.csv("Deseq2_log2fold-NvsT_geneloci_hg38_filtered_normalized.csv", header = T, sep = ",")

#_______________
# Filttrations
#_______________
#keep only annotated genes
rna_exp_matrix <- rna_exp_matrix[ !is.na(rna_exp_matrix$ensembl_id), ]
nrow(rna_exp_matrix)

#remove low counts per sample
 rna_exp_matrix <- rna_exp_matrix[rna_exp_matrix$N18_Proximal + rna_exp_matrix$T18 >= 5 |
                                  rna_exp_matrix$N19_Proximal + rna_exp_matrix$T19 >= 5, ]
 nrow(rna_exp_matrix)

#calculate log ratio for gene expression per sample
rna_exp_matrix["lr18"] = log2((rna_exp_matrix$T18 + 0.01) / (rna_exp_matrix$N18_Proximal + 0.01))
rna_exp_matrix["lr19"] = log2((rna_exp_matrix$T19 + 0.01) / (rna_exp_matrix$N19_Proximal + 0.01))

#plot
#plot_ly(data = data.frame(rna_exp_matrix), x = ~lr18, y = ~lr19 )


#keep only canonical chromosomes
canonical_chr <- paste0("chr", seq(1,22))
canonical_chr <- c(canonical_chr, "chrX") #female samples
rna_exp_matrix <- rna_exp_matrix[rna_exp_matrix$seqnames %in% canonical_chr,]
nrow(rna_exp_matrix)

rna_exp_matrix <- makeGRangesFromDataFrame(df = rna_exp_matrix, seqnames.field = "seqnames", start.field = "start",
                                         end.field = "end", strand.field = "strand", keep.extra.columns = T)

lr_samples <- data.frame(rna_exp_matrix)[c("gene_name", "ensembl_id" , "lr18", "lr19")]


# #plot overall gene expression
 kp <- plotKaryotype(genome = "hg38", chromosomes = "chr20")
 plotLRR(kp, rna_exp_matrix, out.of.range = "density", lrr.column = "lr18", points.col = "blue", r0 = 0, r1 = 0.45)
 plotLRR(kp, rna_exp_matrix, out.of.range = "density", lrr.column = "lr19", points.col = "red", r0 = 0.55, r1 = 1)


####
# Load small mutations dup, with cosmic
####

## Linking samples with dup mutations and diff gene expression
mut18_dups_withHP_cosmic <- read.table("/Path_to/mut18_dups_withHP.csv", header = T, sep = " ")
mut19_dups_withHP_cosmic <- read.table("/Path_to/mut19_dups_withHP.csv", header = T, sep = " ")


# check these genes experssion
mut18_dups_withHP_cosmic$lr18 <- lr_samples$lr18[match(mut18_dups_withHP_cosmic$SYMBOL, lr_samples$gene_name)]
mut19_dups_withHP_cosmic$lr19 <- lr_samples$lr19[match(mut19_dups_withHP_cosmic$SYMBOL, lr_samples$gene_name)]


#write.table(mut18_dups_withHP_cosmic, "mut18_dups_withHP_cosmic_and_exp.csv", sep = "\t")
#write.table(mut19_dups_withHP_cosmic, "mut19_dups_withHP_cosmic_and_exp.csv", sep = "\t")
 
 

##################################################
##        2. Integrate CNA analysis
## intersect with bins from TITAN genomics CNA
##################################################


# plot Copy numbers per chromosome, and color the gene if it has single or multiple mutations

read_cna <- function(cnapath) {
  cna <- read.table(cnapath, header = T, quote = "\t")
  #cna <- na.omit(cna)
  cna <- makeGRangesFromDataFrame(df = cna,
                                   seqnames.field = "chr",
                                  start.field = "Position",
                                  end.field = "Position",
                                  keep.extra.columns = T)
  return(cna)
}

cna_p18 <- read_cna ("/Path_to/CNA/P18/titanCNA_ploidy3/tumor_sample_1_cluster2.titan.txt")
cna_p18 <- cna_p18[,c("Corrected_Ratio", "PhaseSet", "LogRatio", "Corrected_logR", "Corrected_Copy_Number")]

cna_p19 <- read_cna ("/Path_to/CNA/P19/titanCNA_ploidy3/tumor_sample_1_cluster2.titan.txt")
cna_p19 <- cna_p19[,c("Corrected_Ratio", "PhaseSet", "LogRatio", "Corrected_logR", "Corrected_Copy_Number")]


#csv files are generated from a previous step.
muts_18 <- read.table("/Path_to/muts_18_labeled.csv", sep = "\t", header = T)
muts_19 <- read.table("/Path_to/muts_19_labeled.csv", sep = "\t", header = T)


filter_color_mutation_list <- function(muts) {
  muts <- muts[!(colnames(muts) %in% c("cov", "per_hp","hgnc_symbol"))]
          
  muts <- muts %>% filter(all_gene_mut>=2) %>% dplyr::distinct(Hugo_Symbol, .keep_all = T)
          
  muts <- muts[!is.na(muts$end),]
          
  muts <- makeGRangesFromDataFrame(muts, keep.extra.columns = T, start.field = "start",
                                              end.field = "end", seqnames.field = "seqnames")
  muts$color <- NA
  muts[muts$type=="single"]$color <- "red"
  muts[muts$type=="double"]$color <- "darkgreen"
  return(muts)
  }

muts_18 <- filter_color_mutation_list(muts_18)
muts_19 <- filter_color_mutation_list(muts_19)


##############################################################
###           Get the mutations RNA expression             ###
##############################################################


muts_18$lr18 <- lr_samples$lr18[match(muts_18$ENSGene, lr_samples$ensembl_id)]
muts_18$lr18[is.na(muts_18$lr18)] <- 0

muts_19$lr19 <- lr_samples$lr19[match(muts_19$ENSGene, lr_samples$ensembl_id)]
muts_19$lr19[is.na(muts_19$lr19)] <- 0

#can also be by gene name
#muts_18$lr18 <- lr_samples$lr18[match(muts_18$Hugo_Symbol, lr_samples$gene_name)]


plot_selected_karyotype <- function(loc = "chr20",
                                    rna_matrix= selected_gene, rna_sample= "lr19", rna_col = "blue", 
                                    cna_matrix= cna,          cna_sample= "cna_p19", cna_col = "red",
                                    genes_color =  "#333333",
                                    y_max = 3, y_min=-3, plot_exp = T) {

  if (grepl("-", loc[1], fixed = T)) {
    loc <- toGRanges(loc)
    kp  <- plotKaryotype(genome = "hg38", zoom =loc, plot.type = 3)
  }
  else {
    kp  <- plotKaryotype(genome = "hg38", chromosomes=loc, plot.type = 3)
  }

  kpAddBaseNumbers(kp)

  #copy numbers
  plotLRR(kp, cna_matrix , out.of.range = "points", lrr.column = cna_sample, points.cex = 0.4,
          ymax = y_max, ymin= y_min,  points.col = cna_col, labels = NA, line.at.0.col = "gray")
  #RNA
  if (plot_exp) {
  plotLRR(kp, rna_matrix , out.of.range = "points", lrr.column = rna_sample, points.cex = 0.8,
          ymax = y_max, ymin= y_min, points.col = rna_col, out.of.range.col = rna_col, labels = NA, line.at.0.col = "gray")
 }
  kpPlotMarkers(kp, data=rna_matrix, label.color = rna_matrix$color,
                labels = rna_matrix$Hugo_Symbol, label.margin=0.5, data.panel = 2) #,marker.parts=c(0.8,3.1,0.01))
 
}


png("RNA_cp_18.png", width = 12, height = 6,  units = "in",  bg = "white", res = 400)
plot_selected_karyotype(loc= "chr8", 
                        rna_matrix = muts_18, rna_sample = "lr18", rna_col = "black",
                        cna_matrix = cna_p18, cna_sample = "LogRatio", cna_col = "#daeaf5")
dev.off()

png("RNA_cp_19.png", width = 12, height = 6,  units = "in",  bg = "white", res = 400)
plot_selected_karyotype(loc= "chr8", 
                        rna_matrix = muts_19, rna_sample = "lr19", rna_col = "black",
                        cna_matrix = cna_p19, cna_sample = "LogRatio", cna_col = "#daeaf5")
dev.off()

