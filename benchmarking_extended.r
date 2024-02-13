########################################################### Generic #################################################################################################
#####################################################################################################################################################################
libraries <- c("Seurat","patchwork", "tximport", "ggplot2", "scran" , "scater", "SingleCellExperiment", "BUSpaRse", "DropletUtils", "RColorBrewer", "cluster", "UpSetR","scRNAseq", "ggpubr","scDblFinder","reshape")

lapply(libraries, require, character.only = TRUE)
# Import gtf annotation & setWD, source functions
gencode_path <- "~/dato-activo/reference.genomes_kike/GRCh38/gencode/gencode.v37.annotation.gtf"
hg38_ensembl_gtf <- as.data.frame(rtracklayer::import(gencode_path))
hg38_ensembl_gtf$gene_id <- gsub("_","-",hg38_ensembl_gtf$gene_id)

mitochondrial_ens_ids <- unique(hg38_ensembl_gtf$gene_id[grep("^MT-",hg38_ensembl_gtf$gene_name)])
lncrna_ens_ids <- unique(c(hg38_ensembl_gtf$gene_id[grep("lncRNA",hg38_ensembl_gtf$gene_type)]))
protein_coding_ens_ids <- unique(c(hg38_ensembl_gtf$gene_id[hg38_ensembl_gtf$gene_type=="protein_coding"]))

lncrna_names <- unique(hg38_ensembl_gtf$gene_name[hg38_ensembl_gtf$gene_id %in% lncrna_ens_ids])
protein_coding_names <-  unique(hg38_ensembl_gtf$gene_name[hg38_ensembl_gtf$gene_id %in% protein_coding_ens_ids])

source("/home/egonie/kike/phd/git_rep_hpclogin/manuscript_scRNAseq_lncRNAs/manuscript_functions.r")
setwd("/home/egonie/kike/phd/test_data/paper_figures/figure2/")

#mouse
gencode_path_mouse <- "~/dato-activo/reference.genomes_kike/GRCm39/gencode/gencode.vM27.annotation.gtf"
gtf_backup_mouse <- rtracklayer::import(gencode_path_mouse)
gtf_mouse <- as.data.frame(rtracklayer::import(gencode_path_mouse))
gtf_mouse$gene_id <- gsub("_","-",gtf_mouse$gene_id)

mitochondrial_ens_ids_mouse <- unique(gtf_mouse$gene_id[grep("^mt-",gtf_mouse$gene_name)])
lncrna_ens_ids_mouse <- unique(c(gtf_mouse$gene_id[grep("lncRNA",gtf_mouse$gene_type)]))
protein_coding_ens_ids_mouse <- unique(c(gtf_mouse$gene_id[gtf_mouse$gene_type=="protein_coding"]))

lncrna_names_mouse <- unique(gtf_mouse$gene_name[gtf_mouse$gene_id %in% lncrna_ens_ids_mouse])
protein_coding_names_mouse <-  unique(gtf_mouse$gene_name[gtf_mouse$gene_id %in% protein_coding_ens_ids_mouse])

#####################################################################################################################################################################
#####################################################################################################################################################################
# Extend the benchmarking to include more tissues and organisms: Only focus on comparing Kallisto vs cellRanger
#1. Processing

kallisto_intestine_pool1 <- kallisto_processing(kallisto_dir = "/home/egonie/kike/phd/test_data/GSE158702_intestine/01.Kallisto/output_bus_transcriptome_150bp_onlyR2/bustools_results_no_multimappers", name = "cells_genes_NO_multimapping", mito_gene = mitochondrial_ens_ids)
kallisto_intestine_pool2 <- kallisto_processing(kallisto_dir = "/home/egonie/kike/phd/test_data/GSE158702_intestine/GSM4808348/01.Kallisto/output_bus_transcriptome_91bp_reads/bustools_results_no_multimappers", name = "cells_genes_NO_multimapping", mito_gene = mitochondrial_ens_ids)
kallisto_healthy_lung <- kallisto_processing(kallisto_dir = "/home/egonie/kike/phd/test_data/GSE164829_healthy_lung/01.Kallisto/output_bus_transcriptome/bustools_results_no_multimappers", name = "cells_genes_NO_multimapping", mito_gene = mitochondrial_ens_ids)
kallisto_healthy_lung_GSM4037316 <- kallisto_processing(kallisto_dir = "/home/egonie/kike/datos/Marta_GSE135893_PulmonaryFibrosis/GSM4037316/01.Kallisto/output_bus_transcriptome_split_150bp/bustools_results_no_multimappers", name = "cells_genes_NO_multimapping", mito_gene = mitochondrial_ens_ids)
kallisto_pulmonary_fibrosis <- kallisto_processing(kallisto_dir = "/home/egonie/kike/datos/Marta_GSE135893_PulmonaryFibrosis/03.Kallisto/output_bus_transcriptome_split_cDNA/bustools_results_no_multimappers", name = "cells_genes_NO_multimapping", mito_gene = mitochondrial_ens_ids)
kallisto_PBMCs_5K <- kallisto_processing(kallisto_dir = "/home/egonie/kike/phd/test_data/10X/PBMCS_5K_healthy_donor/01.Kallisto/output_bus_transcriptome/bustools_results_no_multimappers", name = "cells_genes_NO_multimapping", mito_gene = mitochondrial_ens_ids)
kallisto_PBMCs_mouse <- kallisto_processing(kallisto_dir = "/home/egonie/kike/phd/test_data/10X/Mouse_PBMCs_Multiplexed_fastqs/SC3_v3_NextGem_DI_CellPlex_Mouse_PBMC_10K/SC3_v3_NextGem_DI_CellPlex_Mouse_PBMC_10K_1_gex/03.Kallisto/output_bus_transcriptome/bustools_results_no_multimappers", name = "cells_genes_NO_multimapping", mito_gene = mitochondrial_ens_ids_mouse)


cellRanger_intestine_pool1 <- cellRanger_processing(cell_ranger_dir = "/home/egonie/kike/phd/test_data/GSE158702_intestine/01.CellRanger/FES_GSM4037320/outs/raw_feature_bc_matrix", mito_gene = mitochondrial_ens_ids)
cellRanger_intestine_pool2 <- cellRanger_processing(cell_ranger_dir = "/home/egonie/kike/phd/test_data/GSE158702_intestine/GSM4808348/01.CellRanger/FES_GSM4808348/outs/raw_feature_bc_matrix", mito_gene = mitochondrial_ens_ids)
cellRanger_healthy_lung <- cellRanger_processing(cell_ranger_dir = "/home/egonie/kike/phd/test_data/GSE164829_healthy_lung/01.CellRanger/FES_GSM5020383/outs/raw_feature_bc_matrix",  mito_gene = mitochondrial_ens_ids)
cellRanger_healthy_lung_GSM4037316 <-  cellRanger_processing(cell_ranger_dir = "/home/egonie/kike/datos/Marta_GSE135893_PulmonaryFibrosis/GSM4037316/01.CellRanger/VUHD70/outs/raw_feature_bc_matrix", mito_gene = mitochondrial_ens_ids)
cellRanger_pulmonary_fibrosis <- cellRanger_processing(cell_ranger_dir = "/home/egonie/kike/datos/Marta_GSE135893_PulmonaryFibrosis/03.CellRanger/VUILD63/outs/raw_feature_bc_matrix", mito_gene = mitochondrial_ens_ids)
cellRanger_PBMCs_5K <- cellRanger_processing(cell_ranger_dir = "/home/egonie/kike/phd/test_data/10X/PBMCS_5K_healthy_donor/01.CellRanger/5k_pbmc_v3_nextgem_alldata/outs/raw_feature_bc_matrix", mito_gene = mitochondrial_ens_ids)
cellRanger_PBMCs_mouse <- cellRanger_processing(cell_ranger_dir = "/home/egonie/kike/phd/test_data/10X/Mouse_PBMCs_Multiplexed_fastqs/SC3_v3_NextGem_DI_CellPlex_Mouse_PBMC_10K/SC3_v3_NextGem_DI_CellPlex_Mouse_PBMC_10K_1_gex/03.CellRanger/raw_feature_bc_matrix", mito_gene = mitochondrial_ens_ids_mouse, multiplexing = T)




