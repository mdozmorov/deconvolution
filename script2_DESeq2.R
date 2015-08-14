# Load the data ant get DEGs
readDDS <- function(DIR) {
  sampleTable <- read.csv(paste(DIR, "sampleTable.csv", sep=""))
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                         directory = DIR,
                                         design = ~ condition)
  ddsHTSeq$condition <- factor(ddsHTSeq$condition)
  return(ddsHTSeq)
}

DIR <- "RNA-seq/04_htseq/"
ddsHTSEQ <- readDDS(DIR)
# write.table(assay(ddsHTSEQ), "RNA-seq/globalCounts38_80.txt", sep="\t", col.names=NA, quote=FALSE)

dds <- DESeq(ddsHTSEQ)
res <- results(dds, contrast = c("condition", "Case", "Control"))
res <- res[complete.cases(res), ]

DEGs <- res[res$padj < 0.25 & abs(res$log2FoldChange) > 0, ]
genes <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description'), filters='ensembl_gene_id', values=rownames(DEGs), mart = mart)
DEGs <- left_join(data.frame(GENE = rownames(DEGs), DEGs), genes, by = c("GENE" = "ensembl_gene_id"))
write.xlsx(DEGs, "RNA-seq/Tables/Table_csSAM.xlsx", sheetName ="DESeq2_counts", append=TRUE)
