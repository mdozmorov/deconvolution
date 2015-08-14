## Load data
G <- as.matrix(read.table("RNA-seq/data/FPKM_matrix.txt", sep="\t", row.names=1, header=TRUE, stringsAsFactors=FALSE))
cc <- as.matrix(t(read.table("RNA-seq/data/Cell_matrix.txt", sep="\t", row.names=1, header=TRUE, stringsAsFactors=FALSE, check.names=TRUE)))
y <- ifelse(grepl("Control", colnames(G)), 1, 2)

## CellMix DSection
# FPKM data
Y <- G
p0 <- t(cc)
# Using 2500 simulations
res.dsection <- DSection(Y, p0, groups=y, W0=10, W_proposal=100, nSamples=500, nBurnIn=2000, summarize=TRUE, verbose=TRUE)
# saveRDS(res.dsection, "RNA-seq/data/res.dsection.Rds")
res.dsection <- readRDS("RNA-seq/data/res.dsection.Rds")
# # Using 1250 simulations
# set.seed(1)
# res.dsection <- CellMix::DSection(Y, p0, groups=y, W0=10, W_proposal=100, nSamples=250, nBurnIn=1000, summarize=TRUE, verbose=TRUE)
# saveRDS(res.dsection, "RNA-seq/data/res.dsection.small.Rds")
# res.dsection <- readRDS("RNA-seq/data/res.dsection.small.Rds")

# # Count data
# Y <- as.matrix(read.table("RNA-seq/globalCounts38_80.txt", sep="\t", row.names=1, header = TRUE, stringsAsFactors = FALSE))
# Y <- Y[ !apply(Y, 1, function(x) sum(x <= 0) >= 1), ]
# # Using 1250 simulations
# set.seed(1)
# res.dsection <- CellMix::DSection(Y, p0, groups=y, W0=10, W_proposal=100, nSamples=250, nBurnIn=1000, summarize=TRUE, verbose=TRUE)
# saveRDS(res.dsection, "RNA-seq/data/res.dsection.count.Rds")
# res.dsection <- readRDS("RNA-seq/data/res.dsection.count.Rds")


# # Random data
# set.seed(1)
# Y.rnd <- mtx.rand(G, randomize = "mix")
# p0.rnd <- t(mtx.rand(cc, randomize = "mix"))
# y.rnd <- sample(y)
# res.dsection.rnd <- CellMix::DSection(Y.rnd, p0.rnd, groups=y.rnd, W0=10, W_proposal=100, nSamples=250, nBurnIn=1000, summarize=TRUE, verbose=TRUE)
# saveRDS(res.dsection.rnd, "RNA-seq/data/res.dsection.rnd.Rds")
# res.dsection <- readRDS("RNA-seq/data/res.dsection.rnd.Rds")
# Random data, row simulations
set.seed(1)
Y.rnd <- G.rnd
p0.rnd <- t(mtx.rand(cc, randomize = "row"))
y.rnd <- sample(y)
res.dsection.rnd <- CellMix::DSection(Y.rnd, p0.rnd, groups=y.rnd, W0=10, W_proposal=100, nSamples=250, nBurnIn=1000, summarize=TRUE, verbose=TRUE)
saveRDS(res.dsection.rnd, "RNA-seq/data/res.dsection.rnd.NB.Rds")
res.dsection <- readRDS("RNA-seq/data/res.dsection.rnd.NB.Rds")
# # Random filtered data, 1250 simulations
# set.seed(1)
# Y.rnd <- t(mtx.rand(G[, !apply(G, 2, function(x) sum(x <= 0) >= 1) ], randomize = "mix"))
# p0.rnd <- t(mtx.rand(cc, randomize = "mix"))
# y.rnd <- sample(y)
# res.dsection.rnd <- CellMix::DSection(Y.rnd, p0.rnd, groups=y.rnd, W0=10, W_proposal=100, nSamples=250, nBurnIn=1000, summarize=TRUE, verbose=TRUE)
# saveRDS(res.dsection.rnd, "RNA-seq/data/res.dsection.rnd.filt.Rds")
# res.dsection <- readRDS("RNA-seq/data/res.dsection.rnd.filt.Rds")

# # Exploring what's inside
# names(res.dsection)
# # "MCData"     "x_LS"       "lambda_LS"  "groups"     "call"       "parameters" "p0" 
# dim(res.dsection$MCData$x)
# # 100   4   2
# res.dsection$MCData$x[1,1,1]
# # 0.263055
# dim(res.dsection$x_LS)
# # 100   4   2
# res.dsection$x_LS[1,1,1]
# # -0.002183831
# dim(res.dsection$lambda_LS)
# 

# Cell type-specific differential expression p-values
n1 <- 10 # Group 1 sample size
n2 <- 10 # Group 2 sample size
res.ttest <- (res.dsection$MCData$x[, , 1] - res.dsection$MCData$x[, , 2]) / as.numeric(sqrt(((1/n1) + (1/n2)) / sqrt(res.dsection$lambda_LS))) # T-statistics

res.pval <- apply(res.ttest, 2, function(x) pt(x, n1 + n2 - 2)) # P-value
res.qval <- apply(res.pval, 2, function(x) p.adjust(x, method="BH")) # Multiple testing correction

apply(res.qval, 2, function(x) length(x[ x < 0.1 ])) # How many significant genes per cell type
res.sig <- apply(res.qval, 2, function(x) (x[ x < 0.1 ])) # Collect them

res.sig.genes <- list() # Annotate differentially expressed genes
for (i in 1:length(res.sig)) {
  genes <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description'), filters='ensembl_gene_id', values=names(res.sig[[i]]), mart = mart)
  res.sig.genes <- c(res.sig.genes, list(left_join(data.frame(V1=names(res.sig[[i]]), PVAL=res.sig[[i]], AVEXP=rowMeans(Y[ names(res.sig[[i]]), , drop=F]), SD=rowSds(Y[ names(res.sig[[i]]), ,drop=F])), genes, by=c("V1" = "ensembl_gene_id"))))
  names(res.sig.genes)[length(res.sig.genes)] <- names(res.sig)[i]
  write.xlsx(res.sig.genes[[i]], "RNA-seq/Tables/Table_csSAM.xlsx", sheetName = paste("DS", names(res.sig)[i], sep="_"), append = TRUE)
}

res.sig.genes[[1]]
res.sig.genes[[2]]
res.sig.genes[[3]]
res.sig.genes[[4]]

# intersect(res.sig.genes[[1]]$hgnc_symbol, c(siggenes.up$hgnc_symbol, siggenes.dn$hgnc_symbol))
# intersect(res.sig.genes[[2]]$hgnc_symbol, c(siggenes.up$hgnc_symbol, siggenes.dn$hgnc_symbol))
# intersect(res.sig.genes[[3]]$hgnc_symbol, c(siggenes.up$hgnc_symbol, siggenes.dn$hgnc_symbol))
# intersect(res.sig.genes[[4]]$hgnc_symbol, c(siggenes.up$hgnc_symbol, siggenes.dn$hgnc_symbol))
