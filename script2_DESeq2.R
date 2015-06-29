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
res.sig <- res[res$padj < 0.1 & abs(res$log2FoldChange) > 0, ]

dim(res)
summary(res)
genes <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description'), filters='ensembl_gene_id', values=rownames(res.sig), mart=mart, uniqueRows=T)
ressig.gene <- merge(as.data.frame(res.sig), genes, by.x="row.names", by.y=1, all.x=T, sort=FALSE)
ressig.gene <- data.frame(ressig.gene, AVEXP=rowMeans(assay(ddsHTSEQ)[ ressig.gene$Row.names, , drop=F]), SD=rowSds(assay(ddsHTSEQ)[ ressig.gene$Row.names, , drop=F]))
# Write results to file
write.xlsx(ressig.gene, "RNA-seq/Tables/Table_csSAM.xlsx", sheetName ="DESeq2_counts", append=TRUE)
