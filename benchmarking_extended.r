########################################################### Generic #################################################################################################
#####################################################################################################################################################################
libraries <- c("Seurat","patchwork", "tximport", "ggplot2", "scran" , "scater", "SingleCellExperiment", "BUSpaRse", "DropletUtils", "RColorBrewer", "cluster", "UpSetR","scRNAseq", "ggpubr","scDblFinder","reshape","ggupset","dplyr")

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
# Extend the benchmarking to include more tissues and organisms: Only focus on comparing Kallisto vs Cell Ranger
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

# EmptyDrops filtering
kallisto_intestine_pool1_ed <- emptydrops_filt(kallisto_intestine_pool1, lower_ED = 100, EmptyDrops_FDR_thres = 0.001 )
kallisto_intestine_pool2_ed <- emptydrops_filt(kallisto_intestine_pool2, lower_ED = 250, EmptyDrops_FDR_thres = 0.001 )
kallisto_healthy_lung_ed <- emptydrops_filt(kallisto_healthy_lung, lower_ED = 100, EmptyDrops_FDR_thres = 0.001 )
kallisto_healthy_lung_GSM4037316_ed <- emptydrops_filt(kallisto_healthy_lung_GSM4037316, lower_ED = 100, EmptyDrops_FDR_thres = 0.01 )
kallisto_pulmonary_fibrosis_ed <- emptydrops_filt(kallisto_pulmonary_fibrosis, lower_ED = 500, EmptyDrops_FDR_thres = 0.001 )
kallisto_PBMCs_5K_ed <- emptydrops_filt(kallisto_PBMCs_5K, lower_ED = 100, EmptyDrops_FDR_thres = 0.001 )
kallisto_PBMCs_mouse_ed <- emptydrops_filt(kallisto_PBMCs_mouse, lower_ED = 100, EmptyDrops_FDR_thres = 0.001 )

cellRanger_intestine_pool1_ed <- emptydrops_filt(cellRanger_intestine_pool1, lower_ED = 100, EmptyDrops_FDR_thres = 0.001 )
cellRanger_intestine_pool2_ed <- emptydrops_filt(cellRanger_intestine_pool2, lower_ED = 250, EmptyDrops_FDR_thres = 0.001 )
cellRanger_healthy_lung_ed <- emptydrops_filt(cellRanger_healthy_lung, lower_ED = 100, EmptyDrops_FDR_thres = 0.001 )
cellRanger_healthy_lung_GSM4037316_ed <- emptydrops_filt(cellRanger_healthy_lung_GSM4037316, lower_ED = 100, EmptyDrops_FDR_thres = 0.01 )
cellRanger_pulmonary_fibrosis_ed <- emptydrops_filt(cellRanger_pulmonary_fibrosis, lower_ED = 500, EmptyDrops_FDR_thres = 0.001 )
cellRanger_PBMCs_5K_ed <- emptydrops_filt(cellRanger_PBMCs_5K, lower_ED = 100, EmptyDrops_FDR_thres = 0.001 )
cellRanger_PBMCs_mouse_ed <- emptydrops_filt(cellRanger_PBMCs_mouse, lower_ED = 100, EmptyDrops_FDR_thres = 0.001 )

#save RDS
kallisto_objects <- list("intestine1" = kallisto_intestine_pool1_ed, "intestine2" = kallisto_intestine_pool2_ed, "healthy_lung1" = kallisto_healthy_lung_ed, "healthy_lung2" = kallisto_healthy_lung_GSM4037316_ed, "pulmonary_fibrosis" = kallisto_pulmonary_fibrosis_ed, "PBMCs_5K" = kallisto_PBMCs_5K_ed, "PBMCs_mouse" = kallisto_PBMCs_mouse_ed)
saveRDS(kallisto_objects, file = "kallisto_datasets.rds")
kallisto_objects <- readRDS("kallisto_datasets.rds")
kallisto_intestine_pool1_ed <- kallisto_objects[["intestine1"]]
kallisto_intestine_pool2_ed <- kallisto_objects[["intestine2"]]
kallisto_healthy_lung_ed <- kallisto_objects[["healthy_lung1"]]
kallisto_healthy_lung_GSM4037316_ed <- kallisto_objects[["healthy_lung2"]]
kallisto_pulmonary_fibrosis_ed <- kallisto_objects[["pulmonary_fibrosis"]]
kallisto_PBMCs_5K_ed <- kallisto_objects[["PBMCs_5K"]]
kallisto_PBMCs_mouse_ed <- kallisto_objects[["PBMCs_mouse"]]

cellRanger_objects <- list("intestine1" = cellRanger_intestine_pool1_ed, "intestine2" = cellRanger_intestine_pool2_ed, "healthy_lung1" = cellRanger_healthy_lung_ed, "healthy_lung2" = cellRanger_healthy_lung_GSM4037316_ed, "pulmonary_fibrosis" = cellRanger_pulmonary_fibrosis_ed, "PBMCs_5K" = cellRanger_PBMCs_5K_ed, "PBMCs_mouse" = cellRanger_PBMCs_mouse_ed)
saveRDS(cellRanger_objects, file = "cellRanger_datasets.rds")
cellRanger_objects <- readRDS("cellRanger_datasets.rds")
cellRanger_intestine_pool1_ed <- cellRanger_objects[["intestine1"]]
cellRanger_intestine_pool2_ed <- cellRanger_objects[["intestine2"]]
cellRanger_healthy_lung_ed <- cellRanger_objects[["healthy_lung1"]]
cellRanger_healthy_lung_GSM4037316_ed <- cellRanger_objects[["healthy_lung2"]]
cellRanger_pulmonary_fibrosis_ed <- cellRanger_objects[["pulmonary_fibrosis"]]
cellRanger_PBMCs_5K_ed <- cellRanger_objects[["PBMCs_5K"]]
cellRanger_PBMCs_mouse_ed <- cellRanger_objects[["PBMCs_mouse"]]

# Doublet analysis
kallisto_intestine_pool1_ed_nodoubs <- doublet_analysis(kallisto_intestine_pool1_ed)
kallisto_intestine_pool1_ed_nodoubs <- kallisto_intestine_pool1_ed_nodoubs[,kallisto_intestine_pool1_ed_nodoubs$isDoublet == F]
cellRanger_intestine_pool1_ed_nodoubs <- doublet_analysis(cellRanger_intestine_pool1_ed)
cellRanger_intestine_pool1_ed_nodoubs <- cellRanger_intestine_pool1_ed_nodoubs[,cellRanger_intestine_pool1_ed_nodoubs$isDoublet == F]

kallisto_intestine_pool2_ed_nodoubs <- doublet_analysis(kallisto_intestine_pool2_ed)
kallisto_intestine_pool2_ed_nodoubs <- kallisto_intestine_pool2_ed_nodoubs[,kallisto_intestine_pool2_ed_nodoubs$isDoublet == F]
cellRanger_intestine_pool2_ed_nodoubs <- doublet_analysis(cellRanger_intestine_pool2_ed)
cellRanger_intestine_pool2_ed_nodoubs <- cellRanger_intestine_pool2_ed_nodoubs[,cellRanger_intestine_pool2_ed_nodoubs$isDoublet == F]

kallisto_healthy_lung_ed_nodoubs <- doublet_analysis(kallisto_healthy_lung_ed)
kallisto_healthy_lung_ed_nodoubs <- kallisto_healthy_lung_ed_nodoubs[,kallisto_healthy_lung_ed_nodoubs$isDoublet == F]
cellRanger_healthy_lung_ed_nodoubs <- doublet_analysis(cellRanger_healthy_lung_ed)
cellRanger_healthy_lung_ed_nodoubs <- cellRanger_healthy_lung_ed_nodoubs[,cellRanger_healthy_lung_ed_nodoubs$isDoublet == F]

kallisto_healthy_lung_GSM4037316_ed_nodoubs <- doublet_analysis(kallisto_healthy_lung_GSM4037316_ed)
kallisto_healthy_lung_GSM4037316_ed_nodoubs <- kallisto_healthy_lung_GSM4037316_ed_nodoubs[,kallisto_healthy_lung_GSM4037316_ed_nodoubs$isDoublet == F]
cellRanger_healthy_lung_GSM4037316_ed_nodoubs <- doublet_analysis(cellRanger_healthy_lung_GSM4037316_ed)
cellRanger_healthy_lung_GSM4037316_ed_nodoubs <- cellRanger_healthy_lung_GSM4037316_ed_nodoubs[,cellRanger_healthy_lung_GSM4037316_ed_nodoubs$isDoublet == F]

kallisto_pulmonary_fibrosis_ed_nodoubs <- doublet_analysis(kallisto_pulmonary_fibrosis_ed)
kallisto_pulmonary_fibrosis_ed_nodoubs <- kallisto_pulmonary_fibrosis_ed_nodoubs[,kallisto_pulmonary_fibrosis_ed_nodoubs$isDoublet == F]
cellRanger_pulmonary_fibrosis_ed_nodoubs <- doublet_analysis(cellRanger_pulmonary_fibrosis_ed)
cellRanger_pulmonary_fibrosis_ed_nodoubs <- cellRanger_pulmonary_fibrosis_ed_nodoubs[,cellRanger_pulmonary_fibrosis_ed_nodoubs$isDoublet == F]

kallisto_PBMCs_5K_ed_nodoubs <- doublet_analysis(kallisto_PBMCs_5K_ed)
kallisto_PBMCs_5K_ed_nodoubs <- kallisto_PBMCs_5K_ed_nodoubs[,kallisto_PBMCs_5K_ed_nodoubs$isDoublet == F]
cellRanger_PBMCs_5K_ed_nodoubs <- doublet_analysis(cellRanger_PBMCs_5K_ed)
cellRanger_PBMCs_5K_ed_nodoubs <- cellRanger_PBMCs_5K_ed_nodoubs[,cellRanger_PBMCs_5K_ed_nodoubs$isDoublet == F]

kallisto_PBMCs_mouse_ed_nodoubs <- doublet_analysis(kallisto_PBMCs_mouse_ed)
kallisto_PBMCs_mouse_ed_nodoubs <- kallisto_PBMCs_mouse_ed_nodoubs[,kallisto_PBMCs_mouse_ed_nodoubs$isDoublet == F]
cellRanger_PBMCs_mouse_ed_nodoubs <- doublet_analysis(cellRanger_PBMCs_mouse_ed)
cellRanger_PBMCs_mouse_ed_nodoubs <- cellRanger_PBMCs_mouse_ed_nodoubs[,cellRanger_PBMCs_mouse_ed_nodoubs$isDoublet == F]

kallisto_objects <- list("intestine1" = kallisto_intestine_pool1_ed_nodoubs, "intestine2" = kallisto_intestine_pool2_ed_nodoubs, "healthy_lung1" = kallisto_healthy_lung_ed_nodoubs, "healthy_lung2" = kallisto_healthy_lung_GSM4037316_ed_nodoubs, "pulmonary_fibrosis" = kallisto_pulmonary_fibrosis_ed_nodoubs, "PBMCs_5K" = kallisto_PBMCs_5K_ed_nodoubs, "PBMCs_mouse" = kallisto_PBMCs_mouse_ed_nodoubs)
saveRDS(kallisto_objects, file = "kallisto_datasets.rds")
kallisto_objects <- readRDS("kallisto_datasets.rds")
kallisto_intestine_pool1_ed <- kallisto_objects[["intestine1"]]
kallisto_intestine_pool2_ed <- kallisto_objects[["intestine2"]]
kallisto_healthy_lung_ed <- kallisto_objects[["healthy_lung1"]]
kallisto_healthy_lung_GSM4037316_ed <- kallisto_objects[["healthy_lung2"]]
kallisto_pulmonary_fibrosis_ed <- kallisto_objects[["pulmonary_fibrosis"]]
kallisto_PBMCs_5K_ed <- kallisto_objects[["PBMCs_5K"]]
kallisto_PBMCs_mouse_ed <- kallisto_objects[["PBMCs_mouse"]]

cellRanger_objects <- list("intestine1" = cellRanger_intestine_pool1_ed_nodoubs, "intestine2" = cellRanger_intestine_pool2_ed_nodoubs, "healthy_lung1" = cellRanger_healthy_lung_ed_nodoubs, "healthy_lung2" = cellRanger_healthy_lung_GSM4037316_ed_nodoubs, "pulmonary_fibrosis" = cellRanger_pulmonary_fibrosis_ed_nodoubs, "PBMCs_5K" = cellRanger_PBMCs_5K_ed_nodoubs, "PBMCs_mouse" = cellRanger_PBMCs_mouse_ed_nodoubs)
saveRDS(cellRanger_objects, file = "cellRanger_datasets.rds")
cellRanger_objects <- readRDS("cellRanger_datasets.rds")
cellRanger_intestine_pool1_ed <- cellRanger_objects[["intestine1"]]
cellRanger_intestine_pool2_ed <- cellRanger_objects[["intestine2"]]
cellRanger_healthy_lung_ed <- cellRanger_objects[["healthy_lung1"]]
cellRanger_healthy_lung_GSM4037316_ed <- cellRanger_objects[["healthy_lung2"]]
cellRanger_pulmonary_fibrosis_ed <- cellRanger_objects[["pulmonary_fibrosis"]]
cellRanger_PBMCs_5K_ed <- cellRanger_objects[["PBMCs_5K"]]
cellRanger_PBMCs_mouse_ed <- cellRanger_objects[["PBMCs_mouse"]]

#####################################################################################################################################################################
#####################################################################################################################################################################
# 2. Quality control 
# Mito content is the only thing that should be plotted before removing the mitocondrial content
color_palette_iwanthue=c("#9473c6","#cc546d")
kallisto_intestine_pool1_ed$sample <- "Hg_intestine_1"
kallisto_intestine_pool2_ed$sample <- "Hg_intestine_2"
kallisto_healthy_lung_ed$sample <- "Hg_lung_1"
kallisto_healthy_lung_GSM4037316_ed$sample <- "Hg_lung_2"
kallisto_pulmonary_fibrosis_ed$sample <- "Hg_pulm_fibrosis"
kallisto_PBMCs_5K_ed$sample <- "Hg_PBMCs_5k"
kallisto_PBMCs_mouse_ed$sample <- "Mm_PBMCs_10k"

cellRanger_intestine_pool1_ed$sample <- "Hg_intestine_1"
cellRanger_intestine_pool2_ed$sample <- "Hg_intestine_2"
cellRanger_healthy_lung_ed$sample <- "Hg_lung_1"
cellRanger_healthy_lung_GSM4037316_ed$sample <- "Hg_lung_2"
cellRanger_pulmonary_fibrosis_ed$sample <- "Hg_pulm_fibrosis"
cellRanger_PBMCs_5K_ed$sample <- "Hg_PBMCs_5k"
cellRanger_PBMCs_mouse_ed$sample <- "Mm_PBMCs_10k"

# add manually mouse data
all_sce <- cbind(cellRanger_intestine_pool1_ed, cellRanger_intestine_pool2_ed, cellRanger_healthy_lung_ed, cellRanger_healthy_lung_GSM4037316_ed, cellRanger_pulmonary_fibrosis_ed, cellRanger_PBMCs_5K_ed, kallisto_intestine_pool1_ed, kallisto_intestine_pool2_ed, kallisto_healthy_lung_ed, kallisto_healthy_lung_GSM4037316_ed, kallisto_pulmonary_fibrosis_ed, kallisto_PBMCs_5K_ed)
all_sce$ident <- as.character(all_sce$ident)
all_sce$ident <- factor(all_sce$ident, levels = c("CellRanger","kallisto"))

# Mito content/cell
color_palette_iwanthue=c("#9473c6","#cc546d")
dodge <- position_dodge(width = 0.85)

pdf("QC_Mito_cell.pdf")
all_sce_melted <- as.data.frame(cbind(c(as.character(all_sce$orig.ident),as.character(cellRanger_PBMCs_mouse_ed$orig.ident),as.character(kallisto_PBMCs_mouse_ed$orig.ident)), c(as.numeric(all_sce$subsets_Mito_percent),as.numeric(cellRanger_PBMCs_mouse_ed$subsets_Mito_percent),as.numeric(kallisto_PBMCs_mouse_ed$subsets_Mito_percent)),c(as.character(all_sce$sample),as.character(cellRanger_PBMCs_mouse_ed$sample),as.character(kallisto_PBMCs_mouse_ed$sample))))
colnames(all_sce_melted) <- c("ident","Mitochondrial_content", "Sample")
all_sce_melted[,2] <- as.numeric(all_sce_melted[,2])
all_sce_melted$ident <- factor(all_sce_melted$ident, levels = c("CellRanger","kallisto"))
ggplot(all_sce_melted, aes(x=Sample, y=Mitochondrial_content, fill = ident)) + geom_violin(alpha = 1, scale = "width", width = 0.8, position = dodge) + scale_fill_manual(values = color_palette_iwanthue) + theme_classic() + scale_colour_manual(values = color_palette_iwanthue) + ggtitle("Mitochondrial content/cell") + theme(plot.title = element_text(hjust = 0.5, size = 20),axis.text=element_text(size=12), axis.title=element_text(size=14, face = "bold"))  + labs(y= "", x="") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(plot.margin = margin(1,0.5,0.5,2, "cm"))
dev.off()

# Filtering & normalization
threshold_mito_percentage = 15
#delete also cells with >50.000 counts
high_threshold_cell_counts = 50000
cells_min_genes_detected_threshold = 500

kallisto_intestine_pool1_ed_filt <- Filtering(kallisto_intestine_pool1_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
kallisto_intestine_pool2_ed_filt <- Filtering(kallisto_intestine_pool2_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
kallisto_healthy_lung_ed_filt <- Filtering(kallisto_healthy_lung_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
kallisto_healthy_lung_GSM4037316_ed_filt <- Filtering(kallisto_healthy_lung_GSM4037316_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
kallisto_pulmonary_fibrosis_ed_filt <- Filtering(kallisto_pulmonary_fibrosis_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
kallisto_PBMCs_5K_ed_filt <- Filtering(kallisto_PBMCs_5K_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
kallisto_PBMCs_mouse_ed_filt <- Filtering(kallisto_PBMCs_mouse_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)

cellRanger_intestine_pool1_ed_filt <- Filtering(cellRanger_intestine_pool1_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
cellRanger_intestine_pool2_ed_filt <- Filtering(cellRanger_intestine_pool2_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
cellRanger_healthy_lung_ed_filt <- Filtering(cellRanger_healthy_lung_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
cellRanger_healthy_lung_GSM4037316_ed_filt <- Filtering(cellRanger_healthy_lung_GSM4037316_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
cellRanger_pulmonary_fibrosis_ed_filt <- Filtering(cellRanger_pulmonary_fibrosis_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
cellRanger_PBMCs_5K_ed_filt <- Filtering(cellRanger_PBMCs_5K_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
cellRanger_PBMCs_mouse_ed_filt <- Filtering(cellRanger_PBMCs_mouse_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)

all_sce <- cbind(cellRanger_intestine_pool1_ed_filt, cellRanger_intestine_pool2_ed_filt, cellRanger_healthy_lung_ed_filt, cellRanger_healthy_lung_GSM4037316_ed_filt, cellRanger_pulmonary_fibrosis_ed_filt, cellRanger_PBMCs_5K_ed_filt, kallisto_intestine_pool1_ed_filt, kallisto_intestine_pool2_ed_filt, kallisto_healthy_lung_ed_filt, kallisto_healthy_lung_GSM4037316_ed_filt, kallisto_pulmonary_fibrosis_ed_filt, kallisto_PBMCs_5K_ed_filt)
all_sce$ident <- as.character(all_sce$ident)
all_sce$ident <- factor(all_sce$ident, levels = c("CellRanger","kallisto"))

# UMIs/cell
pdf("QC_UMIs_cell.pdf")
all_sce_melted <- as.data.frame(cbind(c(as.character(all_sce$orig.ident),as.character(cellRanger_PBMCs_mouse_ed_filt$orig.ident),as.character(kallisto_PBMCs_mouse_ed_filt$orig.ident)), c(as.numeric(all_sce$nCount_RNA),as.numeric(cellRanger_PBMCs_mouse_ed_filt$nCount_RNA),as.numeric(kallisto_PBMCs_mouse_ed_filt$nCount_RNA)),c(as.character(all_sce$sample),as.character(cellRanger_PBMCs_mouse_ed_filt$sample),as.character(kallisto_PBMCs_mouse_ed_filt$sample))))
colnames(all_sce_melted) <- c("ident","nCount_RNA", "Sample")
all_sce_melted[,2] <- as.numeric(all_sce_melted[,2])
all_sce_melted$ident <- factor(all_sce_melted$ident, levels = c("CellRanger","kallisto"))
ggplot(all_sce_melted, aes(x=Sample, y=nCount_RNA, fill = ident)) + geom_violin(alpha = 1, scale = "width", width = 0.8, position = dodge) + scale_fill_manual(values = color_palette_iwanthue) + theme_classic() + scale_colour_manual(values = color_palette_iwanthue) + ggtitle("UMIs/cell") + theme(plot.title = element_text(hjust = 0.5, size = 20),axis.text=element_text(size=12), axis.title=element_text(size=14, face = "bold"))  + labs(y= "", x="") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(plot.margin = margin(1,0.5,0.5,2, "cm"))
dev.off()

#nLncRNAs detected/cell
pdf("QC_detected_lncRNAs_cell.pdf")
all_sce_melted <- as.data.frame(cbind(c(as.character(all_sce$orig.ident),as.character(cellRanger_PBMCs_mouse_ed_filt$orig.ident),as.character(kallisto_PBMCs_mouse_ed_filt$orig.ident)), c(as.numeric(all_sce$subsets_lncRNA_detected),as.numeric(cellRanger_PBMCs_mouse_ed_filt$subsets_lncRNA_detected),as.numeric(kallisto_PBMCs_mouse_ed_filt$subsets_lncRNA_detected)),c(as.character(all_sce$sample),as.character(cellRanger_PBMCs_mouse_ed_filt$sample),as.character(kallisto_PBMCs_mouse_ed_filt$sample))))
colnames(all_sce_melted) <- c("ident","subsets_lncRNA_detected", "Sample")
all_sce_melted[,2] <- as.numeric(all_sce_melted[,2])
all_sce_melted$ident <- factor(all_sce_melted$ident, levels = c("CellRanger","kallisto"))
ggplot(all_sce_melted, aes(x=Sample, y=subsets_lncRNA_detected, fill = ident)) + geom_violin(alpha = 1, scale = "width", width = 0.8, position = dodge) + scale_fill_manual(values = color_palette_iwanthue) + theme_classic() + scale_colour_manual(values = color_palette_iwanthue) + ggtitle("Detected lncRNAs/cell") + theme(plot.title = element_text(hjust = 0.5, size = 20),axis.text=element_text(size=12), axis.title=element_text(size=14, face = "bold"))  + labs(y= "", x="") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(plot.margin = margin(1,0.5,0.5,2, "cm"))
dev.off()

#nProtein_coding detected/cell
pdf("QC_detected_protein_coding_cell.pdf")
all_sce_melted <- as.data.frame(cbind(c(as.character(all_sce$orig.ident),as.character(cellRanger_PBMCs_mouse_ed_filt$orig.ident),as.character(kallisto_PBMCs_mouse_ed_filt$orig.ident)), c(as.numeric(all_sce$subsets_protien_coding_detected),as.numeric(cellRanger_PBMCs_mouse_ed_filt$subsets_protien_coding_detected),as.numeric(kallisto_PBMCs_mouse_ed_filt$subsets_protien_coding_detected)),c(as.character(all_sce$sample),as.character(cellRanger_PBMCs_mouse_ed_filt$sample),as.character(kallisto_PBMCs_mouse_ed_filt$sample))))
colnames(all_sce_melted) <- c("ident","subsets_protien_coding_detected", "Sample")
all_sce_melted[,2] <- as.numeric(all_sce_melted[,2])
all_sce_melted$ident <- factor(all_sce_melted$ident, levels = c("CellRanger","kallisto"))
ggplot(all_sce_melted, aes(x=Sample, y=subsets_protien_coding_detected, fill = ident)) + geom_violin(alpha = 1, scale = "width", width = 0.8, position = dodge) + scale_fill_manual(values = color_palette_iwanthue) + theme_classic() + scale_colour_manual(values = color_palette_iwanthue) + ggtitle("Detected protein coding genes/cell") + theme(plot.title = element_text(hjust = 0.5, size = 20),axis.text=element_text(size=12), axis.title=element_text(size=14, face = "bold"))  + labs(y= "", x="") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(plot.margin = margin(1,0.5,0.5,2, "cm"))
dev.off()

# UpSet plots
all_sce_melted <- as.data.frame(cbind(c(as.character(all_sce$orig.ident),as.character(cellRanger_PBMCs_mouse_ed_filt$orig.ident),as.character(kallisto_PBMCs_mouse_ed_filt$orig.ident)), c(as.character(colnames(all_sce)),as.character(colnames(cellRanger_PBMCs_mouse_ed_filt)),as.character(colnames(kallisto_PBMCs_mouse_ed_filt))),c(as.character(all_sce$sample),as.character(cellRanger_PBMCs_mouse_ed_filt$sample),as.character(kallisto_PBMCs_mouse_ed_filt$sample))))
colnames(all_sce_melted) <- c("ident","cell_id", "Sample")
color_palette_iwanthue_7 <- c("#cb5b42","#48b1a7","#b45ac2","#67a64e","#6f7ccb","#b69340","#c55d88")

all_sce_melted_tbl=tbl_df(all_sce_melted)
pdf("upset_plot_cells.pdf")
all_sce_melted_tbl %>%
  group_by(cell_id,Sample) %>%
  summarize(ident = list(ident)) %>%
  ggplot(aes(x = ident,fill=Sample)) +
    geom_bar(aes(y = after_stat(comp_pct(count, PANEL, fill))), position="dodge")  +
    scale_x_upset() + theme_classic() + scale_fill_manual(values=color_palette_iwanthue_7, name="") + theme(plot.margin = margin(1,0.5,0.5,2, "cm")) + ylab("Percentage of detected cells for each dataset")
dev.off()

# Gene plots
# Two gradients of defined thresholds: 1) 250 counts and present in more than 25 cells and 2) 50 counts and present in more than 5 cells. Also change the features depending on whether the analysis is regarding lncRNAs or protein-coding genes. Modify the names of the .pdf generated depending on your particular analysis.
threshold_minumun_gene_counts = 250
threshold_cells_detected = 25

features_human <- lncrna_ens_ids
features_mouse <- lncrna_ens_ids_mouse
features_mouse <- gsub("\\..*","",features_mouse)

kallisto_intestine_pool1_ed_filt_top_genes <- top_genes(kallisto_intestine_pool1_ed_filt[features_human,],threshold_minumun_gene_counts,threshold_cells_detected )
kallisto_intestine_pool2_ed_filt_top_genes <- top_genes(kallisto_intestine_pool2_ed_filt[features_human,],threshold_minumun_gene_counts,threshold_cells_detected )
kallisto_healthy_lung_ed_filt_top_genes <- top_genes(kallisto_healthy_lung_ed_filt[features_human,],threshold_minumun_gene_counts,threshold_cells_detected )
kallisto_healthy_lung_GSM4037316_ed_filt_top_genes <- top_genes(kallisto_healthy_lung_GSM4037316_ed_filt[features_human,],threshold_minumun_gene_counts,threshold_cells_detected )
kallisto_pulmonary_fibrosis_ed_filt_top_genes <- top_genes(kallisto_pulmonary_fibrosis_ed_filt[features_human,],threshold_minumun_gene_counts,threshold_cells_detected )
kallisto_PBMCs_5K_ed_filt_top_genes <- top_genes(kallisto_PBMCs_5K_ed_filt[features_human,],threshold_minumun_gene_counts,threshold_cells_detected )
rownames(kallisto_PBMCs_mouse_ed_filt) <- gsub("\\..*","",rownames(kallisto_PBMCs_mouse_ed_filt))
kallisto_PBMCs_mouse_ed_filt_top_genes <- top_genes(kallisto_PBMCs_mouse_ed_filt[features_mouse,],threshold_minumun_gene_counts,threshold_cells_detected )

cellRanger_intestine_pool1_ed_filt_top_genes <- top_genes(cellRanger_intestine_pool1_ed_filt[features_human,],threshold_minumun_gene_counts,threshold_cells_detected )
cellRanger_intestine_pool2_ed_filt_top_genes <- top_genes(cellRanger_intestine_pool2_ed_filt[features_human,],threshold_minumun_gene_counts,threshold_cells_detected )
cellRanger_healthy_lung_ed_filt_top_genes <- top_genes(cellRanger_healthy_lung_ed_filt[features_human,],threshold_minumun_gene_counts,threshold_cells_detected )
cellRanger_healthy_lung_GSM4037316_ed_filt_top_genes <- top_genes(cellRanger_healthy_lung_GSM4037316_ed_filt[features_human,],threshold_minumun_gene_counts,threshold_cells_detected )
cellRanger_pulmonary_fibrosis_ed_filt_top_genes <- top_genes(cellRanger_pulmonary_fibrosis_ed_filt[features_human,],threshold_minumun_gene_counts,threshold_cells_detected )
cellRanger_PBMCs_5K_ed_filt_top_genes <- top_genes(cellRanger_PBMCs_5K_ed_filt[features_human,],threshold_minumun_gene_counts,threshold_cells_detected )
cellRanger_PBMCs_mouse_ed_filt_top_genes <- top_genes(cellRanger_PBMCs_mouse_ed_filt[intersect(features_mouse,rownames(cellRanger_PBMCs_mouse_ed_filt)),],threshold_minumun_gene_counts,threshold_cells_detected )

all_sce_melted <- as.data.frame(cbind(as.character(c(rep("kallisto",nrow(kallisto_intestine_pool1_ed_filt_top_genes)),rep("kallisto",nrow(kallisto_intestine_pool2_ed_filt_top_genes)),rep("kallisto",nrow(kallisto_healthy_lung_ed_filt_top_genes)),rep("kallisto",nrow(kallisto_healthy_lung_GSM4037316_ed_filt_top_genes)),rep("kallisto",nrow(kallisto_pulmonary_fibrosis_ed_filt_top_genes)),rep("kallisto",nrow(kallisto_PBMCs_5K_ed_filt_top_genes)),rep("kallisto",nrow(kallisto_PBMCs_mouse_ed_filt_top_genes)),rep("CellRanger",nrow(cellRanger_intestine_pool1_ed_filt_top_genes)),rep("CellRanger",nrow(cellRanger_intestine_pool2_ed_filt_top_genes)),rep("CellRanger",nrow(cellRanger_healthy_lung_ed_filt_top_genes)),rep("CellRanger",nrow(cellRanger_healthy_lung_GSM4037316_ed_filt_top_genes)),rep("CellRanger",nrow(cellRanger_pulmonary_fibrosis_ed_filt_top_genes)),rep("CellRanger",nrow(cellRanger_PBMCs_5K_ed_filt_top_genes)),rep("CellRanger",nrow(cellRanger_PBMCs_mouse_ed_filt_top_genes)))),as.character(c(rownames(kallisto_intestine_pool1_ed_filt_top_genes),rownames(kallisto_intestine_pool2_ed_filt_top_genes), rownames(kallisto_healthy_lung_ed_filt_top_genes),rownames(kallisto_healthy_lung_GSM4037316_ed_filt_top_genes),rownames(kallisto_pulmonary_fibrosis_ed_filt_top_genes),rownames(kallisto_PBMCs_5K_ed_filt_top_genes),rownames(kallisto_PBMCs_mouse_ed_filt_top_genes), rownames(cellRanger_intestine_pool1_ed_filt_top_genes),rownames(cellRanger_intestine_pool2_ed_filt_top_genes), rownames(cellRanger_healthy_lung_ed_filt_top_genes),rownames(cellRanger_healthy_lung_GSM4037316_ed_filt_top_genes),rownames(cellRanger_pulmonary_fibrosis_ed_filt_top_genes),rownames(cellRanger_PBMCs_5K_ed_filt_top_genes),rownames(cellRanger_PBMCs_mouse_ed_filt_top_genes))),as.character(c(rep("Hg_intestine_1",nrow(kallisto_intestine_pool1_ed_filt_top_genes)),rep("Hg_intestine_2",nrow(kallisto_intestine_pool2_ed_filt_top_genes)),rep("Hg_lung_1",nrow(kallisto_healthy_lung_ed_filt_top_genes)),rep("Hg_lung_2",nrow(kallisto_healthy_lung_GSM4037316_ed_filt_top_genes)),rep("Hg_pulm_fibrosis",nrow(kallisto_pulmonary_fibrosis_ed_filt_top_genes)),rep("Hg_PBMCs_5k",nrow(kallisto_PBMCs_5K_ed_filt_top_genes)),rep("Mm_PBMCs_10k",nrow(kallisto_PBMCs_mouse_ed_filt_top_genes)),rep("Hg_intestine_1",nrow(cellRanger_intestine_pool1_ed_filt_top_genes)),rep("Hg_intestine_2",nrow(cellRanger_intestine_pool2_ed_filt_top_genes)),rep("Hg_lung_1",nrow(cellRanger_healthy_lung_ed_filt_top_genes)),rep("Hg_lung_2",nrow(cellRanger_healthy_lung_GSM4037316_ed_filt_top_genes)),rep("Hg_pulm_fibrosis",nrow(cellRanger_pulmonary_fibrosis_ed_filt_top_genes)),rep("Hg_PBMCs_5k",nrow(cellRanger_PBMCs_5K_ed_filt_top_genes)),rep("Mm_PBMCs_10k",nrow(cellRanger_PBMCs_mouse_ed_filt_top_genes))))))

colnames(all_sce_melted) <- c("ident","gene_id", "Sample")

# Represent as %of the total number of features detected in each SCE object
pdf("upset_plot_lncRNAs_250counts_in_25_genes_percentage.pdf")
all_sce_melted_tbl=tbl_df(all_sce_melted)
all_sce_melted_tbl %>%
  group_by(gene_id,Sample) %>%
  summarize(ident = list(ident)) %>%
  ggplot(aes(x = ident,fill=Sample)) +
    geom_bar(aes(y = after_stat(comp_pct(count, PANEL, fill))), position="dodge")  +
    scale_x_upset(order_by = "degree", reverse = T) + theme_classic() + scale_fill_manual(values=color_palette_iwanthue_7, name="") + theme(plot.margin = margin(1,0.5,0.5,2, "cm")) + ylab("Percentage of detected lncRNAs for each dataset")
dev.off()

# Represent as %of the total counts of features detected in each SCE object
pdf("upset_plot_lncRNAs_250counts_in_25_genes_counts.pdf")
all_sce_melted_tbl=tbl_df(all_sce_melted)
all_sce_melted_tbl %>%
  group_by(gene_id,Sample) %>%
  summarize(ident = list(ident)) %>%
  ggplot(aes(x = ident,fill=Sample)) +
    geom_bar(stat="count", position="dodge")  +
    scale_x_upset(order_by = "degree", reverse = T) + theme_classic() + scale_fill_manual(values=color_palette_iwanthue_7, name="") + theme(plot.margin = margin(1,0.5,0.5,2, "cm")) + ylab("Number of detected lncRNAs for each dataset")
dev.off()

