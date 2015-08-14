
# There are several sources for the RPKM / FPKM / TPM formulas:

# TPM function: see "https://www.biostat.wisc.edu/bmi776/lectures/rnaseq.pdf

# https://www.biostars.org/p/11378/
# The formula for FPKM is 10^9 * C / (N * L), 
# with C is the number of mappable reads that fell onto the gene's exons
# N the total number of mappable reads in the experiment and 
# L the number of base pairs in the exon 
# http://genometoolbox.blogspot.com/2014/07/reads-per-kilobase-per-million-mapped.html
# 
# https://bedtools.readthedocs.org/en/latest/content/tools/coverage.html

# https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/

# I think that I trust Lior Pachter more for this:
# https://liorpachter.wordpress.com/2014/04/30/estimating-number-of-transcripts-from-rna-seq-measurements-and-why-i-believe-in-paywall/


###############################################################################
# Here is code for "exonic.gene.sizes" object:
  # Will need this for the FPKM transformation, to get the "effLength" measure:
  # Found here: https://www.biostars.org/p/83901/ on 2 Feb 2015
  # Warning: this takes TIME!!! ...like an hour or so...
  library(GenomicFeatures)
  gtf.file <- "genes38.gtf"
  txdb <- makeTxDbFromGFF( gtf.file ,format="gtf", organism="Homo Sapiens")
  # then collect the exons per gene id
  exons.list.per.gene <- exonsBy(txdb,by="gene")
  # then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
  exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})
###############################################################################

###############################################################################
# This is the FPKM code I put together originally; but now I think this is actually RPKM:
# fPKM function:
FPKM <- function(FPKM.this, exonic.gene.sizes=exonic.gene.sizes){
  col.sums <- apply( FPKM.this, 2 , sum)    # That's right... get the SUM for the ENTIRE sample over ALL genes
  effLength <- unlist(exonic.gene.sizes)     # NOTE: "exonic.gene.sizes" was painfully/time-consumingly gotten from "create_DESeq" scrpt
  # Don't forget to "filter" out genes from "effLength" that were filtered earlier:
  effLength.0 <- effLength[row.names(FPKM.this)]
  n.genes <- dim(FPKM.this)[1]
  n.subs <- dim(FPKM.this)[2]
  rna.fpkm.raw.mc <- matrix(0,ncol=n.subs,nrow=n.genes)
  for (i in 1:n.subs){
    for( j in 1:n.genes){
      # This LOOKS like "RPKM" to me:
      rna.fpkm.raw.mc[j,i] <- (10E9) * ( as.double(FPKM.this[j,i]) / (as.double(col.sums[i]) * as.double(effLength.0[j])))
    } # 1 to number of genes
  } # 1 to number of subjects
  colnames(rna.fpkm.raw.mc) <- colnames(FPKM.this)
  rownames(rna.fpkm.raw.mc) <- names(effLength.0)
  return(rna.fpkm.raw.mc)
  } # End FPKM function
# Now CALL the FPKM function:
# mc.this <- FPKM(FPKM.this) # FPKM.this is created in the "Filter_by..." and "Remove_outlier..." scripts
###############################################################################


