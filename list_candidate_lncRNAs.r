########################################################### Generic #################################################################################################
#####################################################################################################################################################################
libraries <- c("Seurat","patchwork", "tximport", "ggplot2", "scran" , "scater", "SingleCellExperiment", "BUSpaRse", "DropletUtils", "RColorBrewer", "cluster", "UpSetR","scRNAseq", "ggpubr","scDblFinder","reshape","ggupset","dplyr", "ELATUS")

source("/home/egonie/kike/phd/git_rep_hpclogin/manuscript_scRNAseq_lncRNAs/manuscript_functions.r")
setwd("/home/egonie/kike/phd/test_data/paper_figures/list_candidates_lncRNAs/")
lapply(libraries, require, character.only = TRUE)
# Import gtf annotation & setWD, source functions
gencode_path_human <- "~/dato-activo/reference.genomes_kike/GRCh38/gencode/gencode.v37.annotation.gtf"
hg38_ensembl_gtf <- as.data.frame(rtracklayer::import(gencode_path_human))
hg38_ensembl_gtf$gene_id <- gsub("_","-",hg38_ensembl_gtf$gene_id)

lncrna_ens_ids_human <- unique(c(hg38_ensembl_gtf$gene_id[grep("lncRNA",hg38_ensembl_gtf$gene_type)]))
protein_coding_ens_ids_human <- unique(c(hg38_ensembl_gtf$gene_id[hg38_ensembl_gtf$gene_type=="protein_coding"]))
lncrna_names_human <- unique(hg38_ensembl_gtf$gene_name[hg38_ensembl_gtf$gene_id %in% lncrna_ens_ids_human])
protein_coding_names_human <-  unique(hg38_ensembl_gtf$gene_name[hg38_ensembl_gtf$gene_id %in% protein_coding_ens_ids_human])

#mouse
gencode_path_mouse <- "~/dato-activo/reference.genomes_kike/GRCm39/gencode/gencode.vM27.annotation.gtf"
gtf_backup_mouse <- rtracklayer::import(gencode_path_mouse)
gtf_mouse <- as.data.frame(rtracklayer::import(gencode_path_mouse))
gtf_mouse$gene_id <- gsub("_","-",gtf_mouse$gene_id)

lncrna_ens_ids_mouse <- unique(c(gtf_mouse$gene_id[grep("lncRNA",gtf_mouse$gene_type)]))
protein_coding_ens_ids_mouse <- unique(c(gtf_mouse$gene_id[gtf_mouse$gene_type=="protein_coding"]))

lncrna_names_mouse <- unique(gtf_mouse$gene_name[gtf_mouse$gene_id %in% lncrna_ens_ids_mouse])
protein_coding_names_mouse <-  unique(gtf_mouse$gene_name[gtf_mouse$gene_id %in% protein_coding_ens_ids_mouse])

crispr_data <- readRDS("/home/egonie/kike/databases/hits_info_Liu_science_2015_ensids.rds")

#######################################################################################################################
###################### Load datasets from the benchmarking analysis (only V3 chemistry) ###############################
#######################################################################################################################
kallisto_extended_benchmarking_path <- "/home/egonie/kike/phd/test_data/paper_figures/figure2/kallisto_datasets.rds"
cellRanger_extended_benchmarking_path <- "/home/egonie/kike/phd/test_data/paper_figures/figure2/cellRanger_datasets.rds"
kallisto_objects <- readRDS(kallisto_extended_benchmarking_path)
kallisto_intestine_pool1_ed <- kallisto_objects[["intestine1"]]
kallisto_intestine_pool2_ed <- kallisto_objects[["intestine2"]]
kallisto_healthy_lung_ed <- kallisto_objects[["healthy_lung1"]]
kallisto_PBMCs_5K_ed <- kallisto_objects[["PBMCs_5K"]]
kallisto_PBMCs_mouse_ed <- kallisto_objects[["PBMCs_mouse"]]

cellRanger_objects <- readRDS(cellRanger_extended_benchmarking_path)
cellRanger_intestine_pool1_ed <- cellRanger_objects[["intestine1"]]
cellRanger_intestine_pool2_ed <- cellRanger_objects[["intestine2"]]
cellRanger_healthy_lung_ed <- cellRanger_objects[["healthy_lung1"]]
cellRanger_PBMCs_5K_ed <- cellRanger_objects[["PBMCs_5K"]]
cellRanger_PBMCs_mouse_ed <- cellRanger_objects[["PBMCs_mouse"]]

# Same filterings I applied in the benchmarking
threshold_mito_percentage = 15
high_threshold_cell_counts = 50000
threshold_genes_detected = 500

kallisto_intestine_pool1_ed_filt <- Filtering(kallisto_intestine_pool1_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = threshold_genes_detected)
kallisto_intestine_pool2_ed_filt <- Filtering(kallisto_intestine_pool2_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = threshold_genes_detected)
kallisto_healthy_lung_ed_filt <- Filtering(kallisto_healthy_lung_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = threshold_genes_detected)
kallisto_PBMCs_5K_ed_filt <- Filtering(kallisto_PBMCs_5K_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = threshold_genes_detected)
kallisto_PBMCs_mouse_ed_filt <- Filtering(kallisto_PBMCs_mouse_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = threshold_genes_detected)

cellRanger_intestine_pool1_ed_filt <- Filtering(cellRanger_intestine_pool1_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = threshold_genes_detected)
cellRanger_intestine_pool2_ed_filt <- Filtering(cellRanger_intestine_pool2_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = threshold_genes_detected)
cellRanger_healthy_lung_ed_filt <- Filtering(cellRanger_healthy_lung_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = threshold_genes_detected)
cellRanger_PBMCs_5K_ed_filt <- Filtering(cellRanger_PBMCs_5K_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = threshold_genes_detected)
cellRanger_PBMCs_mouse_ed_filt <- Filtering(cellRanger_PBMCs_mouse_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = threshold_genes_detected)

##############################################################################################################################################
###################### ELATUS to define a list of lncRNAs that exhibit characteristics of biologically relevant lncRNAs ######################
##############################################################################################################################################
# 1. Human PBMCs
human_10k_pbmc_path <- "/home/egonie/kike/phd/test_data/paper_figures/figure3_characterization_ex_kallisto/pbmc_datasets_updated.rds"
pbmc_datasets_updated <- readRDS(human_10k_pbmc_path)
cellRanger_sce_filt_clus <- pbmc_datasets_updated[["cellRanger"]]
kallisto_sce_filt_clus <- pbmc_datasets_updated[["kallisto"]]
candidates_human_pbmcs_10k <- ELATUS::get_candidates(kallisto_sce_filt_clus, cellRanger_sce_filt_clus, threshold_minumun_gene_counts = 250, threshold_cells_detected = 25, lncrna_names = lncrna_names_human,gtf = hg38_ensembl_gtf, exclusive = T)

a <- ELATUS_filtered("Human", kallisto_sce_filt_clus, cellRanger_sce_filt_clus, gene_names = TRUE, threshold_minumun_gene_counts = 250, threshold_cells_detected = 25, dimred_clustering = "PCA", k_neighbors = 5, ratio_threshold = 40, CR_threshold = 10, SI_threshold = 0.15)

