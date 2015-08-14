library(samr)
library(csSAM)
library(DESeq2)
library(CellMix)
library(reshape2)

library("xlsx")
library(biomaRt)
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

## Cell type-specific differences
cc <- as.matrix(t(read.table("RNA-seq/data/Cell_matrix.txt", sep="\t", row.names=1, header=TRUE, stringsAsFactors=FALSE, check.names=TRUE)))
apply(cc[ grepl("Control", rownames(cc)), ], 2, mean) # Control mean
apply(cc[ grepl("Case", rownames(cc)), ], 2, mean) # Case mean
apply(cc[ grepl("Case", rownames(cc)), ], 2, mean) - apply(cc[ grepl("Control", rownames(cc)), ], 2, mean)
apply(cc[ grepl("Control", rownames(cc)), ], 2, sd) # Control SD
apply(cc[ grepl("Case", rownames(cc)), ], 2, sd) # Case SD

cc.pval <- list()
for (i in 1:ncol(cc)) {
  res <- wilcox.test(cc[ grepl("Control", rownames(cc)), i], cc[ grepl("Case", rownames(cc)), i], alternative = "two.sided")$p.value
  cc.pval <- c(cc.pval, list(res))
  names(cc.pval)[length(cc.pval)] <- colnames(cc)[i]
}
cc.pval

## Heterogeneous SAM analysis, FPKM
G <- as.matrix(read.table("RNA-seq/data/FPKM_matrix.txt", sep="\t", row.names=1, header=TRUE, stringsAsFactors=FALSE))
y <- ifelse(grepl("Control", colnames(G)), 1, 2)
data.samr <- list(x=G, y=y, geneid=as.character(rownames(G)), logged2=FALSE)
set.seed(1)
samr.obj <- samr(data.samr, resp.type="Two class unpaired")
delta.table <- samr.compute.delta.table(samr.obj)
delta <- 0.4 # FDR 0.2, after FDR goes up to 0.7
samr.plot(samr.obj,delta)
siggenes.table <- samr.compute.siggenes.table(samr.obj, delta, data.samr, delta.table)
siggenes.up <- siggenes.table$genes.up
siggenes.dn <- siggenes.table$genes.lo
siggenes <- rbind(siggenes.up, siggenes.dn)
genes <- getBM(attributes = c('ensembl_gene_id', 'entrezgene', 'hgnc_symbol', 'description'), filters='ensembl_gene_id', values=siggenes[, "Gene Name"], mart = mart)
siggenes <- left_join(data.frame(siggenes, AVEXP=rowMeans(G[ siggenes[, "Gene Name"], , drop=F]), SD=rowSds(G[ siggenes[, "Gene Name"], , drop=F])), genes, by=c("Gene.Name" = "ensembl_gene_id"))
write.xlsx(siggenes, "RNA-seq/Tables/Table_csSAM.xlsx", sheetName ="SAM_FPKM", append=TRUE)

## Heterogeneous SAM analysis, raw counts
G.raw <- as.matrix(read.table("RNA-seq/data/globalCounts38_80.txt", sep="\t", row.names=1, header=TRUE, stringsAsFactors = FALSE))
# sum(apply(G.raw, 1, function(x) sum(x <= 0) >= 1 )) # Check how many zeros
# G.raw <- G.raw[ !apply(G.raw, 1, function(x) sum(x <= 0) >= 1 ), ]
data.samr.raw <- list(x=G.raw, y=y, geneid=as.character(rownames(G.raw)), logged2=FALSE)
set.seed(1)
samr.obj.raw <- samr(data.samr.raw, resp.type="Two class unpaired", assay.type="seq")
delta.table.raw <- samr.compute.delta.table(samr.obj.raw)
# Using all raw data, nothing is significant. Only 7 genes at FDR 0.4365820
# Using filtered data, 16 genes at FDR 0.4764184
delta <- 1  # FDR 0.4, after in goes to 0.68
samr.plot(samr.obj.raw, delta)
siggenes.table <- samr.compute.siggenes.table(samr.obj.raw, delta, data.samr.raw, delta.table.raw)
siggenes.up <- siggenes.table$genes.up # Nothing
siggenes.dn <- siggenes.table$genes.lo
genes <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description'), filters='ensembl_gene_id', values=siggenes.dn[, "Gene Name"], mart = mart)
siggenes <- left_join(data.frame(siggenes.dn, AVEXP=rowMeans(G.raw[ siggenes.dn[, "Gene Name"], ,drop=F]), SD=rowSds(G.raw[ siggenes.dn[, "Gene Name"], , drop=F])), genes, by=c("Gene.Name" = "ensembl_gene_id"))
write.xlsx(siggenes, "RNA-seq/Tables/Table_csSAM.xlsx", sheetName ="SAM_counts", append=TRUE)

## Cell type-specific csSAM analysis, FPKM
set.seed(1)
res.cssam <- csSamWrapper(t(G), cc, y, fileName = "RNA-seq/Figures/Figure_csSAM_FPKM_default.pdf")
# set.seed(1)
# csSamWrapper(G, cc, y, nonNeg = TRUE, fileName = "RNA-seq/Figures/Figure_csSAM_FPKM_nonneg.pdf")
# set.seed(1)
# csSamWrapper(G, cc, y, nonNeg = TRUE, standardize = FALSE, medianCenter = FALSE, fileName = "RNA-seq/Figures/Figure_csSAM_FPKM_nonneg_nostandardize_nomediancenter.pdf")
# set.seed(1)
# csSamWrapper(G, cc, y, nonNeg = FALSE, standardize = FALSE, medianCenter = FALSE, fileName = "RNA-seq/Figures/Figure_csSAM_FPKM_nononneg_nostandardize_nomediancenter.pdf")

## Cell type-specific csSAM results
names(res.cssam) # "deconv","fdr.csSAM","fdr.SAM","sigGene.csSAM","fileName"
res.cssam.sig <- t(res.cssam$sigGene.csSAM)
colnames(res.cssam.sig) <- colnames(cc) # B-cells are index 4
rownames(res.cssam.sig) <- rownames(G)
# B-cells
siggenes <- res.cssam.sig[ res.cssam.sig[, 4] < 0.5, 4, drop=FALSE]
genes <- getBM(attributes=c('ensembl_gene_id', 'entrezgene', 'hgnc_symbol', 'description'), filters='ensembl_gene_id', values=rownames(siggenes), mart=mart)#, uniqueRows=T)
siggenes <- merge(data.frame(siggenes, AVEXP=rowMeans(G[ rownames(siggenes), , drop=F]), SD=rowSds(G[ rownames(siggenes), , drop=F])), genes, by.x="row.names", by.y="ensembl_gene_id", all.x=TRUE, sort=FALSE)
write.xlsx(siggenes, "RNA-seq/Tables/Table_csSAM.xlsx", sheetName = "B-cells", append=TRUE)
hist(log2(rowMeans(G)))
wilcox.test(log2(rowMeans(G)), log2(rowMeans(G[ rownames(siggenes), ])))

# Monocytes
siggenes <- res.cssam.sig[ res.cssam.sig[, 2] < 0.67, 2, drop=FALSE]
genes <- getBM(attributes=c('ensembl_gene_id', 'entrezgene', 'hgnc_symbol', 'description'), filters='ensembl_gene_id', values=rownames(siggenes), mart=mart)#, uniqueRows=T)
siggenes <- merge(data.frame(siggenes, AVEXP=rowMeans(G[ rownames(siggenes), , drop=F]), SD=rowSds(G[ rownames(siggenes), , drop=F])), genes, by.x="row.names", by.y="ensembl_gene_id", all.x=TRUE, sort=FALSE)
write.xlsx(siggenes, "RNA-seq/Tables/Table_csSAM.xlsx", sheetName = "Monocytes", append = TRUE)


## Cell type-specific csSAM analysis, random data
source("simulation_high.r")
source("mtx.rand.R")
set.seed(1)
G.rnd <- matrix(gener.cond(counts = rowMeans(G), noise = 0.1, nrepl = 20), ncol = 20)
# res.sam.rnd <- csSamWrapper(mtx.rand(t(G), randomize = "mix"), mtx.rand(cc, randomize = "mix"), sample(y),fileName = "RNA-seq/Figures/Figure_csSAM_FPKM_rand_mix.pdf")
# set.seed(1)
# res.sam.rnd <- csSamWrapper(mtx.rand(t(G), randomize = "rnd"), mtx.rand(cc, randomize = "mix"), sample(y), fileName = "RNA-seq/Figures/Figure_csSAM_FPKM_rand_rnd.pdf")
set.seed(1)
res.sam.rnd <- csSamWrapper(t(G.rnd), mtx.rand(cc, randomize = "mix"), sample(y), fileName = "RNA-seq/Figures/Figure_csSAM_FPKM_randNB.pdf")


# ## Correlation of cell types with FPKM counts
# # No meaningful results
# cor1 <- apply(G, 1, function(x) rcorr(x, cc[1, ])[[1]][1, 2])
# cor2 <- apply(G, 1, function(x) rcorr(x, cc[2, ])[[1]][1, 2])
# cor3 <- apply(G, 1, function(x) rcorr(x, cc[3, ])[[1]][1, 2])
# cor4 <- apply(G, 1, function(x) rcorr(x, cc[4, ])[[1]][1, 2])
# 
# cor <- cbind(cor1, cor2, cor3, cor4)
# colnames(cor) <- colnames(cc)
# cor <- cor[ apply(cor, 1, function(x) sum(x > 0.5) >= 1 ), ]
# genes <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description'), filters='ensembl_gene_id', values=rownames(cor), mart = mart)
# cor <- merge(cor, genes, by.x="row.names", by.y="ensembl_gene_id", all.x=TRUE, sort=FALSE)
# write.xlsx(cor, "RNA-seq/Tables/Table_csSAM.xlsx", sheetName = "Cor", append = TRUE)
