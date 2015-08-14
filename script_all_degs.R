# http://www.biomedcentral.com/content/supplementary/1471-2105-14-91-s1.pdf

sampleTable <- read.csv("RNA-seq/04_htseq/sampleTable.csv")
class <- sampleTable$condition
count.matrix <- as.matrix(read.table("RNA-seq/data/globalCounts38_80.txt", sep="\t", row.names=1, header=TRUE, stringsAsFactors = FALSE))
count.matrix <- count.matrix[, match(sampleTable$sampleName, colnames(count.matrix))]
count.matrix <- count.matrix[ !apply(count.matrix, 1, function(x) sum(x <= 0) >= 1 ), ]


library(edgeR)
edgeR.dgelist = DGEList(counts = count.matrix, group = factor(class))
edgeR.dgelist = calcNormFactors((edgeR.dgelist), method = "TMM")
edgeR.dgelist = estimateCommonDisp(edgeR.dgelist)
edgeR.dgelist = estimateTagwiseDisp(edgeR.dgelist, trend = "movingave")
edgeR.test = exactTest(edgeR.dgelist)
edgeR.pvalues = edgeR.test$table$PValue
edgeR.adjpvalues = p.adjust(edgeR.pvalues, method = "BH")

library(DESeq) 
DESeq.cds = newCountDataSet(countData = count.matrix, conditions = factor(class)) 
DESeq.cds = estimateSizeFactors(DESeq.cds) 
DESeq.cds = estimateDispersions(DESeq.cds, sharingMode = "maximum", method = "pooled", fitType = "local") 
DESeq.test = nbinomTest(DESeq.cds, "Control", "Case") 
DESeq.pvalues = DESeq.test$pval 
DESeq.adjpvalues = p.adjust(DESeq.pvalues, method = "BH")

library(edgeR) 
library(NBPSeq) 
NBPSeq.dgelist = DGEList(counts = count.matrix, group = factor(class))
NBPSeq.dgelist = calcNormFactors(NBPSeq.dgelist, method = "TMM") 
NBPSeq.norm.factors = as.vector(NBPSeq.dgelist$samples$norm.factors) 
NBPSeq.test = nbp.test(counts = count.matrix, grp.ids = class, grp1 = "Control", grp2 = "Case", norm.factors = NBPSeq.norm.factors, model.disp = "NBP") 
NBPSeq.pvalues = NBPSeq.test$p.values 
NBPSeq.adjpvalues = NBPSeq.test$q.values

library(baySeq) 
baySeq.cd = new("countData", data = count.matrix, replicates = class, groups = list(NDE = rep(1, length(class)), DE = class)) 
baySeq.cd@libsizes = getLibsizes(baySeq.cd, estimationType = "edgeR") 
baySeq.cd = getPriors.NB(baySeq.cd, samplesize = 5000, equalDispersions = TRUE, estimation = "QL", cl = NULL) 
baySeq.cd = getLikelihoods.NB(baySeq.cd, prs = c(0.5, 0.5), pET = "BIC", cl = NULL) 
baySeq.posteriors.DE = exp(baySeq.cd@posteriors)[, 2] 
baySeq.table = topCounts(baySeq.cd, group = "DE", FDR = 1) 
baySeq.FDR = baySeq.table$FDR[match(rownames(count.matrix), rownames(baySeq.table))]

library(EBSeq)
sizes = MedianNorm(count.matrix) 
EBSeq.test = EBTest(Data = count.matrix, Conditions = factor(class), sizeFactors = sizes, maxround = 10) 
EBSeq.ppmat = GetPPMat(EBSeq.test) 
EBSeq.probabilities.DE = EBSeq.ppmat[, "PPDE"] 
EBSeq.lFDR = 1 - EBSeq.ppmat[, "PPDE"] 
EBSeq.FDR = rep(NA, length(EBSeq.lFDR)) 
for (i in 1:length(EBSeq.lFDR)) { EBSeq.FDR[i] = mean(EBSeq.lFDR[which(EBSeq.lFDR <= EBSeq.lFDR[i])]) }

library(edgeR) 
# http://www.stat.purdue.edu/~doerge/software/TSPM.R
source("TSPM.R") 
TSPM.dgelist = DGEList(counts = count.matrix, group = factor(class)) 
TSPM.dgelist = calcNormFactors(TSPM.dgelist, method = "TMM") 
norm.lib.sizes = as.vector(TSPM.dgelist$samples$norm.factors) * as.vector(TSPM.dgelist$samples$lib.size) 
TSPM.test = TSPM(counts = count.matrix, x1 = factor(class), x0 = rep(1, length(class)), lib.size = norm.lib.sizes) 
TSPM.pvalues = TSPM.test$pvalues 
TSPM.adjpvalues = TSPM.test$padj

library(samr) 
SAMseq.test = SAMseq(count.matrix, as.numeric(factor(class)), resp.type = "Two class unpaired", geneid = rownames(count.matrix), genenames = rownames(count.matrix), nperms = 100, nresamp = 20, fdr.output = 1) 
SAMseq.result.table = rbind(SAMseq.test$siggenes.table$genes.up, SAMseq.test$siggenes.table$genes.lo) 
SAMseq.score = rep(0, nrow(count.matrix)) 
SAMseq.score[match(SAMseq.result.table[, 1], rownames(count.matrix))] = as.numeric(SAMseq.result.table[, 3]) 
SAMseq.FDR = rep(1, nrow(count.matrix)) 
SAMseq.FDR[match(SAMseq.result.table[, 1], rownames(count.matrix))] = as.numeric(SAMseq.result.table[, 5])/100

library(edgeR) 
# http://bioinfo.cipf.es/noiseq/doku.php?id=downloads
source("noiseq.r") 
library("NOISeq")
mtx <- NOISeq::readData(data = count.matrix)
nf = calcNormFactors(count.matrix) 
libsizes = apply(count.matrix, 2, sum) 
common.libsize = prod(libsizes^(1/length(libsizes)))
normfactors = nf * libsizes/common.libsize 
norm.matrix = sweep(count.matrix, 2, normfactors, "/") 
NOISeq.test = noiseq(norm.matrix[, class == 1], norm.matrix[, class == 2], repl = "bio", k = 0.5, norm = "n", long = 1000) 
NOISeq.probabilities = NOISeq.test$probab

library(limma) 
nf = calcNormFactors(count.matrix, method = "TMM") 
voom.data = voom(count.matrix, design = model.matrix(~factor(class)), lib.size = colSums(count.matrix) * nf) 
voom.data$genes = rownames(count.matrix) 
voom.fitlimma = lmFit(voom.data, design = model.matrix(~factor(class))) 
voom.fitbayes = eBayes(voom.fitlimma) 
voom.pvalues = voom.fitbayes$p.value[, 2] 
voom.adjpvalues = p.adjust(voom.pvalues, method = "BH")

library(DESeq) 
library(limma) 
DESeq.cds = newCountDataSet(countData = count.matrix, conditions = factor(class)) 
DESeq.cds = estimateSizeFactors(DESeq.cds) 
DESeq.cds = estimateDispersions(DESeq.cds, method = "blind", fitType = "local") 
DESeq.vst = getVarianceStabilizedData(DESeq.cds) 
DESeq.vst.fitlimma = lmFit(DESeq.vst, design = model.matrix(~factor(class))) 
DESeq.vst.fitbayes = eBayes(DESeq.vst.fitlimma) 
DESeq.vst.pvalues = DESeq.vst.fitbayes$p.value[, 2] 
DESeq.vst.adjpvalues = p.adjust(DESeq.vst.pvalues, method = "BH")

library(ShrinkBayes) 
library(edgeR)
nf = calcNormFactors(count.matrix, method = "TMM") * colSums(count.matrix)/exp(mean(log(colSums(count.matrix)))) 
count.matrix = round(sweep(count.matrix, 2, nf, "/")) 
group = factor(class) 
form = y ~ 1 + group 
ShrinkSeq.shrinkres = ShrinkSeq(form = form, dat = count.matrix, shrinkfixed = "group", mixtdisp = FALSE, shrinkdisp = TRUE, fams = "zinb", ncpus = 1) 
ShrinkSeq.fitzinb = FitAllShrink(forms = form, dat = count.matrix, fams = "zinb", shrinksimul = ShrinkSeq.shrinkres, ncpus = 1) 
ShrinkSeq.npprior = NonParaUpdatePrior(fitall = ShrinkSeq.fitzinb, modus = "fixed", shrinkpara = "group", maxiter = 15, ncpus = 1) 
ShrinkSeq.nppostshr = NonParaUpdatePosterior(ShrinkSeq.fitzinb, ShrinkSeq.npprior, ncpus = 1) 
ShrinkSeq.lfdrless = SummaryWrap(ShrinkSeq.nppostshr, thr = 0, direction = "lesser")
ShrinkSeq.lfdrgreat = SummaryWrap(ShrinkSeq.nppostshr, thr = 0, direction = "greater") 
ShrinkSeq.FDR = BFDR(ShrinkSeq.lfdrless, ShrinkSeq.lfdrgreat)
