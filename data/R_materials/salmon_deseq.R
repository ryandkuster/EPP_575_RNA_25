library(tidyverse)
library(tximport)
library(GenomicFeatures)
library(pheatmap)
library(DESeq2)
library(R.utils)

BiocManager::install("DEGreport")
library(DEGreport)

if (requireNamespace("rstudioapi", quietly = TRUE)) {
  # Get active document context
  doc_context <- rstudioapi::getActiveDocumentContext()
  # Extract directory path
  script_dir <- dirname(doc_context$path)
} else {
  script_dir <- normalizePath(path = commandArgs()[1])
}
script_dir = dirname(script_dir)
setwd(script_dir)

# give the full path to the "salmon_results" (it should end in "salmon_results"
salmon_dir <- paste0(script_dir, "/salmon_output")

# load the gff3 file, then create a transcript database/dataframe for use with deseq
download.file("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/735/GCA_000001735.2_TAIR10.1/GCA_000001735.2_TAIR10.1_genomic.gff.gz", destfile = "GCA_000001735.2_TAIR10.1_genomic.gff.gz", method = "wget")
gunzip("GCA_000001735.2_TAIR10.1_genomic.gff.gz", remove = TRUE)

txdb <- makeTxDbFromGFF("GCA_000001735.2_TAIR10.1_genomic.gff")
keytypes(txdb)
k <- keys(txdb, keytype = "CDSNAME")
str(k)

txdf = AnnotationDbi::select(txdb, k, "GENEID", "CDSNAME")

samples <- read_csv(paste0(script_dir, "/salmon_output/salmon_data.csv"))
Qfiles <- file.path(salmon_dir, samples$quant_file)

# this step imports the count data from salmon
txi <- tximport(files = Qfiles, type = "salmon", txOut = TRUE)
head(txi$counts)
colnames(txi$counts) <- samples$sample_id

# convert fields to factors
samples$treatment = factor(samples$treatment)
samples$time = factor(samples$time)
samples$replicate = factor(samples$replicate)

# now we convert the txi object into a deseq-formatted object
dds <- DESeqDataSetFromTximport(txi = txi, colData = samples, design = ~ treatment + time + treatment:time)
dds <- DESeq(dds, test="LRT", reduced = ~ treatment + time)
dds <- dds[which(mcols(dds)$fullBetaConv),]
res_LRT <- results(dds)

# plot dispersion
plotDispEsts(dds)
vsd <- vst(dds)
plotPCA(vsd, intgroup = c("time"))
plotPCA(vsd, intgroup = c("replicate"))
plotPCA(vsd, intgroup = c("treatment"))

################################################################################
# summarize results
res <- results(dds)
head(res)
summary(res)

# create a contrast with lfcThreshold and alpha cutoff (first list item is condition from samples object)
# here are contrasts we can do:

res_treatment <- results(dds, alpha = 0.05, contrast = c("treatment", "Col", "minusEGTA"))
res_time <- results(dds, alpha = 0.05, contrast = c("time", "3h", "0h"))

plotMA(res_treatment, ylim=c(-12,12))
plotMA(res_time, ylim=c(-12,12))

res_sig <- as.data.frame(res_treatment[ which(res_treatment$padj < 0.05),])
res_sig <- res_sig[order(res_sig$padj, -abs(res_sig$log2FoldChange)),]

write.csv(res_sig, file="salmon_featurecounts_sig_results.csv")




