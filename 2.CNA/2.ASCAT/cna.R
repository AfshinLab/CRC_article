library(maftools)

setwd("/Path/to/1-ASCAT/P18-or-P19")

local_gtMarkers <- function (t_bam = NULL, n_bam = NULL, build = "hg19", prefix = NULL, 
          add = TRUE, mapq = 10, sam_flag = 1024, loci = NULL, fa = NULL, 
          op = NULL, zerobased = FALSE, nthreads = 4, verbose = TRUE) 
{
  if (is.null(t_bam)) 
    stop("Missing tumor BAM file!")
  bam = c(t_bam)
  if (!is.null(n_bam)) {
    bam = c(bam, n_bam)
  }
  if (is.null(loci)) {
    if (build == "hg19") {
      download.file(url = "https://github.com/CompEpigen/ezASCAT/blob/main/inst/extdata/GRCh37_SNP6.tsv.gz?raw=true", 
                    destfile = "GRCh37_SNP6.tsv.gz")
      loci = "GRCh37_SNP6.tsv.gz"
    }
    else {
      #download.file(url = "https://github.com/CompEpigen/ezASCAT/blob/main/inst/extdata/GRCh38_SNP6.tsv.gz?raw=true", 
      #              destfile = "GRCh38_SNP6.tsv.gz")
      loci = "/Path/to/1-ASCAT/GRCh38_SNP6.tsv.gz"
    }
  }
  loci = data.table::fread(input = loci)
  colnames(loci)[1:2] = c("Chr", "start")
  if (!is.null(prefix)) {
    if (add) {
      loci$Chr = paste(prefix, loci$Chr, sep = "")
    }
    else {
      loci$Chr = gsub(pattern = prefix, replacement = "", 
                      x = loci$Chr, fixed = TRUE)
    }
  }
  data.table::setDF(x = loci)
  if (zerobased) {
    loci$start = as.numeric(loci$start) + 1
  }
  op_files = lapply(bam, function(x) {
    bam_ext = substr(x = basename(x), start = nchar(basename(path = x)) - 
                       3, nchar(basename(x)))
    if (bam_ext != ".bam" && bam_ext != "cram") {
      stop("Input file is not a BAM or CRAM file: ", x)
    }
    if (!file.exists(x)) {
      stop("BAM file does not exist: ", x)
    }
    gsub(pattern = "\\.bam$", replacement = "", x = basename(x), 
         ignore.case = TRUE)
  })
  if (is.null(op)) {
    op = as.character(unlist(op_files))
    op_files = lapply(op, function(x) {
      paste0(x, "_nucleotide_counts")
    })
    op_files = as.character(unlist(op_files))
  }
  else {
    if (length(op) != length(bam)) {
      stop("No. of output file names must be equal to no. of BAM files.")
    }
    op_files = paste0(op, "_nucleotide_counts")
  }
  if (all(file.exists(op_files))) {
    warning("Counts are already generated!")
    res = lapply(seq_along(op_files), function(x) {
      data.table::fread(file = paste0(op_files[x], ".tsv"), 
                        sep = "\t", header = TRUE)
    })
    names(res) = op
    return(res)
  }
  if (is.null(fa)) {
    fa = "NULL"
  }
  loci = split(loci, loci$Chr)
  loci_files = lapply(1:length(loci), function(idx) {
    chrname = names(loci)[idx]
    lfile = tempfile(pattern = paste0(chrname, "_"), fileext = paste0("_loci.tsv"))
    data.table::fwrite(x = loci[[idx]][, c(1:2)], file = lfile, 
                       col.names = FALSE, sep = "\t", row.names = FALSE)
    lfile
  })
  if (verbose) {
    cat("Fetching readcounts from BAM files..\n")
  }
  res = list()
  bam_idxstats = list()
  for (b in bam) {
    if (verbose) {
      cat("Processing", basename(b), ":\n")
    }
    bam_counts = parallel::mclapply(loci_files, function(lfile) {
      chr = unlist(data.table::tstrsplit(basename(path = lfile), 
                                         split = "_", keep = 1))
      if (verbose) {
        system(paste("echo ' current chromosome:", chr, 
                     "'"))
      }
      opcount = tempfile(pattern = paste0(chr, "_", basename(b)), 
                         fileext = ".tsv")
      withCallingHandlers(suppressWarnings(invisible(.Call("snpc", 
                                                           b, lfile, mapq, sam_flag, fa, opcount, PACKAGE = "maftools"))))
      paste0(opcount, ".tsv")
    }, mc.cores = nthreads)
    idxstat = apply(data.table::fread(file = bam_counts[[1]], 
                                      nrow = 1, sep = "\t"), 1, paste, collapse = " ")
    bam_idxstats[[length(bam_idxstats) + 1]] = idxstat
    res[[length(res) + 1]] = data.table::rbindlist(lapply(bam_counts, 
                                                          data.table::fread), use.names = TRUE, fill = TRUE)
    lapply(bam_counts, unlink)
  }
  names(res) = op
  lapply(seq_along(res), function(idx) {
    cat(paste0(bam_idxstats[[idx]], "\n"), file = paste0(op_files[[idx]], 
                                                         ".tsv"))
    data.table::fwrite(x = res[[idx]], file = paste0(op_files[[idx]], 
                                                     ".tsv"), append = TRUE, sep = "\t", na = "NA", quote = FALSE, 
                       col.names = TRUE)
  })
}

#Highlighted as takes long time and was done, the output file is saved
#careful: if the bam files have the same name then the output is ONE file only
counts = local_gtMarkers(    t_bam = "T.final.phased.cram",
                              n_bam = "N.final.phased.cram",
                              build = "hg38", prefix = 'chr', nthreads = 16)


saveRDS(object = counts, file = "counts.rds")

#ASCAT analysis
#Prepare input 

library(ASCAT)
ascat.bc = maftools::prepAscat(t_counts = "T.final.phased.cram_nucleotide_counts.tsv",
                               n_counts = "N.final.phased.cram_nucleotide_counts.tsv",
                               sample_name = "tumor")

ascat.bc = ASCAT::ascat.loadData(
  Tumor_LogR_file = "tumor.tumour.logR.txt",
  Tumor_BAF_file = "tumor.tumour.BAF.txt",
  Germline_LogR_file = "tumor.normal.logR.txt",
  Germline_BAF_file = "tumor.normal.BAF.txt",
  chrs = c(1:22, "X"),
  sexchromosomes = c("X")
)

ASCAT::ascat.plotRawData(ASCATobj = ascat.bc, img.prefix = "tumor")
ascat.bc = ASCAT::ascat.aspcf(ascat.bc)
ASCAT::ascat.plotSegmentedData(ascat.bc)
ascat.output = ASCAT::ascat.runAscat(ascat.bc) 


sessionInfo()
