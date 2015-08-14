library(limma)
library(edgeR)

sampleTable <- read.csv(paste(DIR, "sampleTable.csv", sep=""))
G.raw <- as.matrix(read.table("RNA-seq/data/globalCounts38_80.txt", sep="\t", row.names=1, header=TRUE, stringsAsFactors = FALSE))
G.raw <- G.raw[, match(sampleTable$sampleName, colnames(G.raw))]

## Heterogeneous Limma-voom analysis
G.raw <- G.raw[ !apply(G.raw, 1, function(x) sum(x <= 0) >= 1 ), ]
# Prepare design matrix
Group <- as.factor(sampleTable$condition)
Group <- relevel(Group, ref="Control")
design <- model.matrix(~Group)
colnames(design) <- c("Control", "CasevsControl")
# Differential expression
dge <- DGEList(counts = G.raw)
dge <- calcNormFactors(dge)
v <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
topTable(fit, coef="CasevsControl", adjust="none")
DEGs <- topTable(fit, coef = "CasevsControl", number = nrow(dge), adjust.method = "none", p.value = 0.01)
genes <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description'), filters='ensembl_gene_id', values=rownames(DEGs), mart = mart)
DEGs <- left_join(data.frame(GENE = rownames(DEGs), DEGs), genes, by = c("GENE" = "ensembl_gene_id"))
write.xlsx(DEGs, "RNA-seq/Tables/Table_csSAM.xlsx", sheetName ="Limma_voom", append=TRUE)
# Analysis with combined voom and sample quality weights
plotMDS(v, labels=y, col=y, main="MDS plot")
legend("bottomright", legend=c("Control", "Case"), col=1:2, pch=15)
vwts <- voomWithQualityWeights(dge, design=design, normalization="none", plot=TRUE)
vfit2 <- lmFit(vwts)
vfit2 <- eBayes(vfit2)
topTable(vfit2,coef="CasevsControl",sort.by="P")
DEGs <- topTable(vfit2, coef = "CasevsControl", number = nrow(dge), adjust.method = "none", p.value = 0.01)
genes <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description'), filters='ensembl_gene_id', values=rownames(DEGs), mart = mart)
DEGs <- left_join(data.frame(GENE = rownames(DEGs), DEGs), genes, by = c("GENE" = "ensembl_gene_id"))
write.xlsx(DEGs, "RNA-seq/Tables/Table_csSAM.xlsx", sheetName ="Limma_voomwts", append=TRUE)

# edgeR classic heterogeneous analysis
y <- DGEList(counts=G.raw, group=Group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
topTags(et)
DEGs <- as.data.frame(topTags(et, n = nrow(y), adjust.method = "none"))
DEGs <- DEGs[ DEGs$PValue < 0.01, ]
genes <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description'), filters='ensembl_gene_id', values=rownames(DEGs), mart = mart)
DEGs <- left_join(data.frame(GENE = rownames(DEGs), DEGs), genes, by = c("GENE" = "ensembl_gene_id"))
write.xlsx(DEGs, "RNA-seq/Tables/Table_csSAM.xlsx", sheetName ="edgeR_classic", append=TRUE)
# edgeR glm heterogeneous analysis
y <- DGEList(counts=G.raw)
design <- model.matrix(~Group)
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)
DEGs <- as.data.frame(topTags(lrt, n = nrow(y), adjust.method = "none"))
DEGs <- DEGs[ DEGs$PValue < 0.01, ]
genes <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description'), filters='ensembl_gene_id', values=rownames(DEGs), mart = mart)
DEGs <- left_join(data.frame(GENE = rownames(DEGs), DEGs), genes, by = c("GENE" = "ensembl_gene_id"))
write.xlsx(DEGs, "RNA-seq/Tables/Table_csSAM.xlsx", sheetName ="edgeR_glm", append=TRUE)

# NOISeq heterogeneous analysis
library("NOISeq")
mtx <- NOISeq::readData(data = count.matrix, factors=sampleTable)
NOISeq.test <- noiseqbio(mtx, k = 0.5, norm = "tmm", factor = "condition")
DEGs = degenes(NOISeq.test, q = 0.8, M = NULL)
genes <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description'), filters='ensembl_gene_id', values=rownames(DEGs), mart = mart)
DEGs <- left_join(data.frame(GENE = rownames(DEGs), DEGs), genes, by = c("GENE" = "ensembl_gene_id"))
write.xlsx(DEGs, "RNA-seq/Tables/Table_csSAM.xlsx", sheetName ="NOISeq", append=TRUE)
