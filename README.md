# Uncovering functional lncRNAs by scRNA-seq with ELATUS
Repository containing all the scripts and code necessary for generating the results in the manuscript [Uncovering functional lncRNAs by scRNA-seq with ELATUS](https://www.nature.com/articles/s41467-024-54005-7).

# System Requirements
## Hardware requirements
Some of the scRNA-seq datasets included in the manuscript are quite computationally heavy to process. Therefore a HPC with enough memory, cores and RAM is required to end-to-end reproduce all the results. 

## Software requirements
### OS Requirements
The different analysis have been performed on the following systems:
+ Linux: CentOS 7.15.1804

### ScRNA-seq preprocessing pipelines 
+ Cell Ranger 3.0.1
+ STAR 2.7.9
+ Kallisto 0.46.1
+ Bustools 0.40.0
+ Salmon 1.4.0
+ Kallisto 0.46.2 (pulmonary fibrotic dataset only)
+ Cell Ranger 5.0.1 (single-cell multiome only)

### R dependencies
```{r}
install.packages("Seurat")
install.packages("patchwork")
install.packages("ggplot2")
install.packages("RColorBrewer")
install.packages("cluster")
install.packages("UpSetR")
install.packages("ggpubr")
install.packages("reshape")
BiocManager::install("tximport")
BiocManager::install("scran")
BiocManager::install("scater")
BiocManager::install("SingleCellExperiment")
BiocManager::install("BUSpaRse")
BiocManager::install("DropletUtils")
BiocManager::install("scDblFinder")
BiocManager::install("scRNAseq")
install.packages("Signac")
```

## Session information
```{r}
> sessionInfo()
R version 4.1.2 (2021-11-01)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /beegfs/easybuild/CentOS/7.5.1804/Skylake/software/FlexiBLAS/3.0.4-GCC-11.2.0/lib64/libflexiblas.so.3.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] Signac_1.9.0                reshape_0.8.9              
 [3] scDblFinder_1.8.0           ggpubr_0.4.0               
 [5] scRNAseq_2.8.0              UpSetR_1.4.0               
 [7] cluster_2.1.4               RColorBrewer_1.1-3         
 [9] DropletUtils_1.14.2         BUSpaRse_1.8.0             
[11] scater_1.22.0               scran_1.22.1               
[13] scuttle_1.4.0               SingleCellExperiment_1.16.0
[15] SummarizedExperiment_1.24.0 Biobase_2.54.0             
[17] GenomicRanges_1.46.1        GenomeInfoDb_1.30.1        
[19] IRanges_2.28.0              S4Vectors_0.32.4           
[21] BiocGenerics_0.40.0         MatrixGenerics_1.6.0       
[23] matrixStats_0.62.0          ggplot2_3.4.2              
[25] tximport_1.22.0             patchwork_1.1.2            
[27] SeuratObject_4.1.3          Seurat_4.0.1               

loaded via a namespace (and not attached):
  [1] rappdirs_0.3.3                rtracklayer_1.54.0           
  [3] scattermore_0.8               R.methodsS3_1.8.2            
  [5] tidyr_1.2.1                   bit64_4.0.5                  
  [7] irlba_2.3.5.1                 DelayedArray_0.20.0          
  [9] R.utils_2.12.0                data.table_1.14.2            
 [11] rpart_4.1.16                  KEGGREST_1.34.0              
 [13] RCurl_1.98-1.9                AnnotationFilter_1.18.0      
 [15] generics_0.1.3                GenomicFeatures_1.46.5       
 [17] ScaledMatrix_1.2.0            cowplot_1.1.1                
 [19] RSQLite_2.2.18                RANN_2.6.1                   
 [21] future_1.28.0                 bit_4.0.4                    
 [23] spatstat.data_3.0-0           xml2_1.3.3                   
 [25] httpuv_1.6.6                  assertthat_0.2.1             
 [27] viridis_0.6.2                 hms_1.1.2                    
 [29] promises_1.2.0.1              fansi_1.0.3                  
 [31] restfulr_0.0.15               progress_1.2.2               
 [33] dbplyr_2.2.1                  igraph_1.3.5                 
 [35] DBI_1.1.3                     htmlwidgets_1.5.4            
 [37] spatstat.geom_3.0-6           purrr_1.0.2                  
 [39] ellipsis_0.3.2                backports_1.4.1              
 [41] dplyr_1.0.10                  biomaRt_2.50.3               
 [43] deldir_1.0-6                  sparseMatrixStats_1.6.0      
 [45] vctrs_0.6.5                   ensembldb_2.18.4             
 [47] ROCR_1.0-11                   abind_1.4-5                  
 [49] cachem_1.0.6                  withr_2.5.0                  
 [51] BSgenome_1.62.0               progressr_0.11.0             
 [53] sctransform_0.3.5             GenomicAlignments_1.30.0     
 [55] prettyunits_1.1.1             goftest_1.2-3                
 [57] ExperimentHub_2.2.1           lazyeval_0.2.2               
 [59] crayon_1.5.2                  edgeR_3.36.0                 
 [61] pkgconfig_2.0.3               nlme_3.1-160                 
 [63] vipor_0.4.5                   ProtGenerics_1.26.0          
 [65] rlang_1.1.1                   globals_0.16.1               
 [67] lifecycle_1.0.3               miniUI_0.1.1.1               
 [69] filelock_1.0.2                BiocFileCache_2.2.1          
 [71] rsvd_1.0.5                    AnnotationHub_3.2.2          
 [73] polyclip_1.10-0               lmtest_0.9-40                
 [75] Matrix_1.5-1                  carData_3.0-5                
 [77] Rhdf5lib_1.16.0               zoo_1.8-11                   
 [79] beeswarm_0.4.0                ggridges_0.5.4               
 [81] png_0.1-7                     viridisLite_0.4.1            
 [83] rjson_0.2.21                  bitops_1.0-7                 
 [85] R.oo_1.25.0                   KernSmooth_2.23-20           
 [87] rhdf5filters_1.6.0            Biostrings_2.62.0            
 [89] blob_1.2.3                    DelayedMatrixStats_1.16.0    
 [91] stringr_1.4.1                 parallelly_1.32.1            
 [93] spatstat.random_3.1-3         rstatix_0.7.0                
 [95] ggsignif_0.6.3                beachmat_2.10.0              
 [97] scales_1.2.1                  memoise_2.0.1                
 [99] magrittr_2.0.3                plyr_1.8.7                   
[101] ica_1.0-3                     zlibbioc_1.40.0              
[103] compiler_4.1.2                dqrng_0.3.0                  
[105] BiocIO_1.4.0                  fitdistrplus_1.1-8           
[107] Rsamtools_2.10.0              cli_3.6.2                    
[109] XVector_0.34.0                listenv_0.8.0                
[111] pbapply_1.5-0                 MASS_7.3-58.1                
[113] mgcv_1.8-40                   tidyselect_1.2.0             
[115] stringi_1.7.8                 yaml_2.3.5                   
[117] BiocSingular_1.10.0           locfit_1.5-9.6               
[119] ggrepel_0.9.1                 grid_4.1.2                   
[121] fastmatch_1.1-3               tools_4.1.2                  
[123] future.apply_1.9.1            parallel_4.1.2               
[125] bluster_1.4.0                 metapod_1.2.0                
[127] gridExtra_2.3                 plyranges_1.14.0             
[129] Rtsne_0.16                    digest_0.6.29                
[131] BiocManager_1.30.19           shiny_1.7.2                  
[133] Rcpp_1.0.9                    car_3.1-0                    
[135] broom_1.0.1                   BiocVersion_3.14.0           
[137] later_1.3.0                   RcppAnnoy_0.0.20             
[139] httr_1.4.4                    AnnotationDbi_1.56.2         
[141] colorspace_2.0-3              XML_3.99-0.11                
[143] tensor_1.5                    reticulate_1.26              
[145] splines_4.1.2                 RcppRoll_0.3.0               
[147] uwot_0.1.14                   statmod_1.4.37               
[149] spatstat.utils_3.0-1          sp_1.5-0                     
[151] xgboost_1.6.0.1               plotly_4.10.0                
[153] xtable_1.8-4                  jsonlite_1.8.2               
[155] zeallot_0.1.0                 R6_2.5.1                     
[157] pillar_1.8.1                  htmltools_0.5.3              
[159] mime_0.12                     glue_1.6.2                   
[161] fastmap_1.1.0                 BiocParallel_1.28.3          
[163] BiocNeighbors_1.12.0          interactiveDisplayBase_1.32.0
[165] codetools_0.2-18              utf8_1.2.2                   
[167] lattice_0.20-45               spatstat.sparse_3.0-0        
[169] tibble_3.1.8                  curl_4.3.3                   
[171] ggbeeswarm_0.6.0              leiden_0.4.3                 
[173] survival_3.4-0                limma_3.50.3                 
[175] munsell_0.5.0                 rhdf5_2.38.1                 
[177] GenomeInfoDbData_1.2.7        HDF5Array_1.22.1             
[179] reshape2_1.4.4                gtable_0.3.1                 
[181] spatstat.core_2.4-4
```







