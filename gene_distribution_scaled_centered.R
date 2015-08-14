# mtx <- read.table("RNA-seq/globalCounts.txt", sep="\t", row.names=1, header=T, stringsAsFactors = FALSE)
# saveRDS(mtx[1:1000, ], "RNA-seq/globalCounts_1000subset.Rds")
mtx <- readRDS("RNA-seq/globalCounts_1000subset.Rds")
summary(mtx[, 1])
summary(as.numeric(mtx[1, ]))
library(vioplot)
vioplot(as.numeric(mtx[1, ]))

mtx <- mtx[ apply(mtx, 1, sd) != 0, ] # Remove genes with SD=0 (likely all zeros)
mtx <- t(apply(mtx, 1, scale)) # Scale/center _genes_ (horizontally)

library(reshape2)
mtx <- melt(t(mtx)) # Stack genes one after the other

pdf("RNA-seq/gene_distribution_scaled_centered_pdf")
vioplot(mtx$value)
title("Scaled/centered gene expression distribution")
dev.off()
# ========================================================================
# http://vinaykmittal.blogspot.com/2013/10/fpkmrpkm-normalization-caveat-and-upper.html
# Normalize raw counts
mtx <- read.table("RNA-seq/globalCounts.txt", sep="\t", row.names=1, header=T, stringsAsFactors = FALSE)
mtx <- mtx[ apply(mtx, 1, function(x) { sum(x == 0) != length(x) }), ] # Remove genes with all zero counts
mtx.3q <- summary(mtx)[5, ] # Extract 75% quartile
mtx.3q <- as.numeric(sub("3rd Qu.:", "", mtx.3q)) 

# ========================================================================
# FPKM counts
mtx <- read.table("RNA-seq/genes.fpkm_table", sep="\t", header=T, row.names=1, stringsAsFactors=F) # 65,199 genes total
sum(apply(mtx, 1, function(x) sum(x <= 0) >= 1 )) # 30,866 are all 0; 41,077 are 0 50% of more times; 46,502 have at least one 0; 52,551 have at least one expression equal to 1
mtx <- mtx[ !apply(mtx, 1, function(x) sum(x <= 0) >= 1 ), ] # 34,333 are non-zero; 24,122 left after removing 50% or more zero genes; 18,697 left after removing genes with at least one 0; 12,648 have expression above 1
write.table(mtx, "RNA-seq/Tables/Table_S1_FPKM_matrix.txt", sep="\t", quote=FALSE, col.names=NA)

boxplot(mtx, ylim=c(0, 100))

pdf("RNA-seq/Figures/Figure_S1_FPKM_distribution.pdf", height=5)
hist(rowMeans(mtx), n=100000, xlim=c(0,25), xlab="Average FPKM")
dev.off()
summary(rowMeans(mtx))

for (i in 1:10) {
  barplot(as.numeric(mtx[i,]), n=20)
  cat ("Press [enter] to continue")
  line <- readline()
}

