library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(pathview)
library(DOSE)


BiocManager::install("org.At.tair.db", character.only = TRUE)
BiocManager::install("gage")
library(gage)
library(pathview)
library("org.At.tair.db", character.only = TRUE)

# Arabidopsis Database for GO
keytypes(org.At.tair.db)


############################################################
# import DGE data
############################################################
setwd("~/Downloads/EPP_575_RNA_25/DESEQ/")
df = read.csv("star_featurecounts_results.csv", header=TRUE)
names(df)[names(df) == "X"] <- "gene"
df$Name <- sub("^gene-", "", df$gene)
l2f_gene_list <- df$log2FoldChange
names(l2f_gene_list) <- df$Name
gene_list<-na.omit(l2f_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)


############################################################
# GO approach
############################################################
gsea <- gseGO(geneList=gene_list, 
              ont ="ALL", 
              keyType = "TAIR",  
              minGSSize = 3, 
              maxGSSize = 1000, 
              pvalueCutoff = 0.01, 
              verbose = TRUE, 
              OrgDb = org.At.tair.db,
              pAdjustMethod = "BH")

summary(as.data.frame(gsea))
dotplot(gsea, showCategory=10, split=".sign") + facet_grid(.~.sign)

# optional, but useful, plots
ridgeplot(gsea) + labs(x = "enrichment distribution")

for (i in 1:10) {
  p <- gseaplot(gsea, by = "all", title = gsea$Description[i], geneSetID = i)
  print(p)
}

terms <- gsea$Description[1:4]
pmcplot(terms, 2010:2024, proportion=FALSE)

############################################################
# gage approach
############################################################
kg.ath <- kegg.gsets(species="ath")
kegg_df = kg.ath[["kg.sets"]]
unique_entries <- unique(unlist(kegg_df))
binary_matrix <- sapply(kegg_df, function(x) as.integer(unique_entries %in% x))
result_df <- as.data.frame(binary_matrix, row.names = unique_entries)
colnames(result_df) <- names(kegg_df)
result_df <- tibble::rownames_to_column(result_df, var = "Name")
merged_df <- merge(df, result_df, by="Name", all.x=TRUE)


############################################################
# KEGG approach
############################################################
kegg_species <- clusterProfiler:::kegg_species_data()
ath <- search_kegg_organism("ath", by='kegg_code')

# over-representation test (only uses gene names)
res_sig = subset(df, abs(df$log2FoldChange) > 2)
res_sig = subset(res_sig, res_sig$padj < 0.05)
gene <- res_sig$Name
length(gene)
"AT4G25480" %in% gene

kk <- enrichKEGG(gene,
                 organism="ath",
                 pvalueCutoff=0.05,
                 pAdjustMethod="BH")
head(kk)

browseKEGG(kk, 'ath04075')
browseKEGG(kk, 'ath04626')


# gseKEGG gene set enrichment analysis (gsea) (takes in log2FC with gene names)
kk2 <- gseKEGG(geneList      = gene_list,
               organism      = "ath",
               pvalueCutoff  = 0.05,
               verbose       = FALSE,
               pAdjustMethod = "fdr")

summary(as.data.frame(kk2))
dotplot(kk2, showCategory=10, split=".sign") + facet_grid(.~.sign)

############################################################
# work-around for kegg ath database not importing by default
############################################################
korg <- cbind("ktax.id" = "T00041", "tax.id" = "3702", "kegg.code" = "ath",
              "scientific.name" = "Arabidopsis thaliana", "common.name" = "thale cress",
              "entrez.gnodes" = "1", "kegg.geneid" = "816394", "ncbi.geneid" = "816394",
              "ncbi.proteinid" = "NP_001325249.1 ", "uniprot" = "A0A178VXX1")

pathview(gene.data  = gene_list,
         pathway.id = "ath04075",
         species    = "ath")
