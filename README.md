Scripts for performing cell type-specific gene expression analysis
---

Tested under R version 3.2.0 (2015-04-16)

The following packages should be available:

attached base packages:

 [1] parallel  stats4    compiler  graphics  grDevices datasets  utils     stats     methods   base     

other attached packages:

 [1] biomaRt_2.24.0            xlsx_0.5.7                xlsxjars_0.6.1            rJava_0.9-6              
 [5] CellMix_1.6.2             GSEABase_1.30.2           graph_1.46.0              annotate_1.46.0          
 [9] XML_3.98-1.2              AnnotationDbi_1.30.1      stringr_1.0.0             NMF_0.22                 
[13] Biobase_2.28.0            cluster_2.0.2             rngtools_1.2.4            pkgmaker_0.25.9          
[17] registry_0.2              DESeq2_1.8.1              RcppArmadillo_0.5.200.1.0 Rcpp_0.11.6              
[21] GenomicRanges_1.20.5      GenomeInfoDb_1.4.1        IRanges_2.2.4             S4Vectors_0.6.0          
[25] BiocGenerics_0.14.0       csSAM_1.2.4               samr_2.0                  matrixStats_0.14.1       
[29] impute_1.42.0             ggplot2_1.0.1             dplyr_0.4.2               BiocInstaller_1.18.3     

loaded via a namespace (and not attached):

 [1] splines_3.2.0         foreach_1.4.2         gtools_3.5.0          Formula_1.2-1         assertthat_0.1       
 [6] latticeExtra_0.6-26   RSQLite_1.0.0         lattice_0.20-31       quadprog_1.5-5        digest_0.6.8         
[11] RColorBrewer_1.1-2    XVector_0.8.0         colorspace_1.2-6      preprocessCore_1.30.0 plyr_1.8.3           
[16] lpSolve_5.6.11        bibtex_0.4.0          genefilter_1.50.0     xtable_1.7-4          corpcor_1.6.7        
[21] scales_0.2.5          whisker_0.3-2         BiocParallel_1.2.5    nnet_7.3-9            proto_0.3-10         
[26] survival_2.38-2       magrittr_1.5          doParallel_1.0.8      MASS_7.3-41           foreign_0.8-63       
[31] beeswarm_0.2.0        tools_3.2.0           gridBase_0.4-7        munsell_0.4.2         locfit_1.5-9.1       
[36] limSolve_1.5.5.1      lambda.r_1.1.7        RCurl_1.95-4.6        futile.logger_1.4.1   grid_3.2.0           
[41] iterators_1.0.7       bitops_1.0-6          gtable_0.1.2          codetools_0.2-11      DBI_0.3.1            
[46] reshape2_1.4.1        R6_2.0.1              gridExtra_0.9.1       Hmisc_3.16-0          futile.options_1.0.0 
[51] dendextend_0.18.3     stringi_0.5-2         geneplotter_1.46.0    rpart_4.1-9           acepack_1.3-3.3  