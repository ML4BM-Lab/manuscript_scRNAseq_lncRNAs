########################################################### Generic #################################################################################################
#####################################################################################################################################################################
libraries <- c("Seurat","patchwork", "tximport", "ggplot2", "scran" , "scater", "SingleCellExperiment", "BUSpaRse", "DropletUtils", "RColorBrewer", "cluster", "UpSetR","scRNAseq", "ggpubr","scDblFinder","reshape", "dplyr", "ggupset", "batchelor")

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
mitochondrial_names <- unique(hg38_ensembl_gtf$gene_name[grep("^MT-",hg38_ensembl_gtf$gene_name)])

source("/home/egonie/kike/phd/git_rep_hpclogin/manuscript_scRNAseq_lncRNAs/manuscript_functions.r")
setwd("/home/egonie/kike/datos/119_Amaya_scRNAseq_BLANCA/")

#####################################################################################################################################################################
# 1. Processing: Sample1: AA017
cell_ranger_dir <- "/home/egonie/kike/datos/119_Amaya_scRNAseq_BLANCA/NextSeq2000.RUN62.20220309/01.CellRanger/AA017_alldata/outs/raw_feature_bc_matrix"
cellRanger_data <- Read10X(data.dir = cell_ranger_dir, gene.column = 1)
cellRanger_AA017 <- CreateSeuratObject(cellRanger_data, project = "cellRanger")

kallisto_dir <- "/home/egonie/kike/datos/119_Amaya_scRNAseq_BLANCA/NextSeq2000.RUN62.20220309/output_bus_transcriptome/bustools_results_no_multimappers"
kallisto_data <- BUSpaRse::read_count_output(kallisto_dir, name = "cells_genes_NO_multimapping")
kallisto_AA017 <- CreateSeuratObject(kallisto_data, project = "kallisto")

# Convert to SingleCellExperiment
cellRanger_AA017 <- as.SingleCellExperiment(cellRanger_AA017)
colnames(cellRanger_AA017) <- gsub("\\-.*","",colnames(cellRanger_AA017))
kallisto_AA017 <- as.SingleCellExperiment(kallisto_AA017)

# Quality control
cellRanger_sce_AA017 <- qc_metrics(cellRanger_AA017, mitochondrial_ens_ids)
kallisto_sce_AA017 <- qc_metrics(kallisto_AA017, mitochondrial_ens_ids)

# EmptyDrops filtering
cellRanger_filt_sce_ed_AA017 <- emptydrops_filt(cellRanger_sce_AA017, lower_ED = 500, EmptyDrops_FDR_thres = 0.01 )
kallisto_filt_sce_ed_AA017 <- emptydrops_filt(kallisto_sce_AA017, lower_ED = 500, EmptyDrops_FDR_thres = 0.01 )

threshold_mito_percentage = 20
low_threshold_cell_counts = 3000
high_threshold_cell_counts = 50000
cells_min_genes_detected_threshold = 700
cellRanger_filt_sce_AA017 <- Filtering_TNBC(cellRanger_filt_sce_ed_AA017,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold, cells_min_threshold = low_threshold_cell_counts)
kallisto_filt_sce_AA017 <- Filtering_TNBC(kallisto_filt_sce_ed_AA017,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold, cells_min_threshold = low_threshold_cell_counts)

# uniquifyFeatures
gene_name <- hg38_ensembl_gtf$gene_name[match(rownames(cellRanger_filt_sce_AA017),hg38_ensembl_gtf$gene_id)]
rowData(cellRanger_filt_sce_AA017) <- cbind(ens_id = rownames(cellRanger_filt_sce_AA017),gene_name)
rownames(cellRanger_filt_sce_AA017) <- uniquifyFeatureNames(rownames(cellRanger_filt_sce_AA017), rowData(cellRanger_filt_sce_AA017)$gene_name)

gene_name <- hg38_ensembl_gtf$gene_name[match(rownames(kallisto_filt_sce_AA017),hg38_ensembl_gtf$gene_id)]
rowData(kallisto_filt_sce_AA017) <- cbind(ens_id = rownames(kallisto_filt_sce_AA017),gene_name)
rownames(kallisto_filt_sce_AA017) <- uniquifyFeatureNames(rownames(kallisto_filt_sce_AA017), rowData(kallisto_filt_sce_AA017)$gene_name)

#####################################################################################################################################################################
# 1. Processing: Sample2: AA024
cell_ranger_dir <- "/home/egonie/kike/datos/119_Amaya_scRNAseq_BLANCA/NextSeq2000.RUN74.20220512/AA024_alldata/outs/raw_feature_bc_matrix"
cellRanger_data <- Read10X(data.dir = cell_ranger_dir, gene.column = 1)
cellRanger_AA024 <- CreateSeuratObject(cellRanger_data, project = "cellRanger")

kallisto_dir <- "/home/egonie/kike/datos/119_Amaya_scRNAseq_BLANCA/NextSeq2000.RUN74.20220512/output_bus_transcriptome/bustools_results_no_multimappers"
kallisto_data <- BUSpaRse::read_count_output(kallisto_dir, name = "cells_genes_NO_multimapping")
kallisto_AA024 <- CreateSeuratObject(kallisto_data, project = "kallisto")

# Convert to SingleCellExperiment
cellRanger_AA024 <- as.SingleCellExperiment(cellRanger_AA024)
colnames(cellRanger_AA024) <- gsub("\\-.*","",colnames(cellRanger_AA024))
kallisto_AA024 <- as.SingleCellExperiment(kallisto_AA024)

# Quality control
cellRanger_sce_AA024 <- qc_metrics(cellRanger_AA024, mitochondrial_ens_ids)
kallisto_sce_AA024 <- qc_metrics(kallisto_AA024, mitochondrial_ens_ids)

# EmptyDrops filtering
cellRanger_filt_sce_ed_AA024 <- emptydrops_filt(cellRanger_sce_AA024, lower_ED = 500, EmptyDrops_FDR_thres = 0.01 )
kallisto_filt_sce_ed_AA024 <- emptydrops_filt(kallisto_sce_AA024, lower_ED = 500, EmptyDrops_FDR_thres = 0.01 )

threshold_mito_percentage = 20
low_threshold_cell_counts = 3000
high_threshold_cell_counts = 50000
cells_min_genes_detected_threshold = 700
cellRanger_filt_sce_AA024 <- Filtering_TNBC(cellRanger_filt_sce_ed_AA024,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold, cells_min_threshold = low_threshold_cell_counts)
kallisto_filt_sce_AA024 <- Filtering_TNBC(kallisto_filt_sce_ed_AA024,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold, cells_min_threshold = low_threshold_cell_counts)

# uniquifyFeatures
gene_name <- hg38_ensembl_gtf$gene_name[match(rownames(cellRanger_filt_sce_AA024),hg38_ensembl_gtf$gene_id)]
rowData(cellRanger_filt_sce_AA024) <- cbind(ens_id = rownames(cellRanger_filt_sce_AA024),gene_name)
rownames(cellRanger_filt_sce_AA024) <- uniquifyFeatureNames(rownames(cellRanger_filt_sce_AA024), rowData(cellRanger_filt_sce_AA024)$gene_name)

gene_name <- hg38_ensembl_gtf$gene_name[match(rownames(kallisto_filt_sce_AA024),hg38_ensembl_gtf$gene_id)]
rowData(kallisto_filt_sce_AA024) <- cbind(ens_id = rownames(kallisto_filt_sce_AA024),gene_name)
rownames(kallisto_filt_sce_AA024) <- uniquifyFeatureNames(rownames(kallisto_filt_sce_AA024), rowData(kallisto_filt_sce_AA024)$gene_name)

#####################################################################################################################################################################
# 1. Processing: Sample3: AA025
cell_ranger_dir <- "/home/egonie/kike/datos/119_Amaya_scRNAseq_BLANCA/NextSeq2000.RUN81.aacacho/01.CellRanger/AA025_alldata/outs/raw_feature_bc_matrix"
cellRanger_data <- Read10X(data.dir = cell_ranger_dir, gene.column = 1)
cellRanger_AA025 <- CreateSeuratObject(cellRanger_data, project = "cellRanger")

kallisto_dir <- "/home/egonie/kike/datos/119_Amaya_scRNAseq_BLANCA/NextSeq2000.RUN81.aacacho/01.Kallisto/output_bus_transcriptome/bustools_results_no_multimappers"
kallisto_data <- BUSpaRse::read_count_output(kallisto_dir, name = "cells_genes_NO_multimapping")
kallisto_AA025 <- CreateSeuratObject(kallisto_data, project = "kallisto")

# Convert to SingleCellExperiment
cellRanger_AA025 <- as.SingleCellExperiment(cellRanger_AA025)
colnames(cellRanger_AA025) <- gsub("\\-.*","",colnames(cellRanger_AA025))
kallisto_AA025 <- as.SingleCellExperiment(kallisto_AA025)

# Quality control
cellRanger_sce_AA025 <- qc_metrics(cellRanger_AA025, mitochondrial_ens_ids)
kallisto_sce_AA025 <- qc_metrics(kallisto_AA025, mitochondrial_ens_ids)

# EmptyDrops filtering
cellRanger_filt_sce_ed_AA025 <- emptydrops_filt(cellRanger_sce_AA025, lower_ED = 500, EmptyDrops_FDR_thres = 0.01 )
kallisto_filt_sce_ed_AA025 <- emptydrops_filt(kallisto_sce_AA025, lower_ED = 500, EmptyDrops_FDR_thres = 0.01 )

threshold_mito_percentage = 20
low_threshold_cell_counts = 3000
high_threshold_cell_counts = 50000
cells_min_genes_detected_threshold = 700
cellRanger_filt_sce_AA025 <- Filtering_TNBC(cellRanger_filt_sce_ed_AA025,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold, cells_min_threshold = low_threshold_cell_counts)
kallisto_filt_sce_AA025 <- Filtering_TNBC(kallisto_filt_sce_ed_AA025,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold, cells_min_threshold = low_threshold_cell_counts)

# uniquifyFeatures
gene_name <- hg38_ensembl_gtf$gene_name[match(rownames(cellRanger_filt_sce_AA025),hg38_ensembl_gtf$gene_id)]
rowData(cellRanger_filt_sce_AA025) <- cbind(ens_id = rownames(cellRanger_filt_sce_AA025),gene_name)
rownames(cellRanger_filt_sce_AA025) <- uniquifyFeatureNames(rownames(cellRanger_filt_sce_AA025), rowData(cellRanger_filt_sce_AA025)$gene_name)

gene_name <- hg38_ensembl_gtf$gene_name[match(rownames(kallisto_filt_sce_AA025),hg38_ensembl_gtf$gene_id)]
rowData(kallisto_filt_sce_AA025) <- cbind(ens_id = rownames(kallisto_filt_sce_AA025),gene_name)
rownames(kallisto_filt_sce_AA025) <- uniquifyFeatureNames(rownames(kallisto_filt_sce_AA025), rowData(kallisto_filt_sce_AA025)$gene_name)

#####################################################################################################################################################################
# 1. Processing: Sample4: AA038
cell_ranger_dir <- "/home/egonie/kike/datos/119_Amaya_scRNAseq_BLANCA/NextSeq2000.RUN109.aacacho/01.CellRanger/AA038_alldata/outs/raw_feature_bc_matrix"
cellRanger_data <- Read10X(data.dir = cell_ranger_dir, gene.column = 1)
cellRanger_AA038 <- CreateSeuratObject(cellRanger_data, project = "cellRanger")

kallisto_dir <- "/home/egonie/kike/datos/119_Amaya_scRNAseq_BLANCA/NextSeq2000.RUN109.aacacho/01.Kallisto/output_bus_transcriptome/bustools_results_no_multimappers"
kallisto_data <- BUSpaRse::read_count_output(kallisto_dir, name = "cells_genes_NO_multimapping")
kallisto_AA038 <- CreateSeuratObject(kallisto_data, project = "kallisto")

# Convert to SingleCellExperiment
cellRanger_AA038 <- as.SingleCellExperiment(cellRanger_AA038)
colnames(cellRanger_AA038) <- gsub("\\-.*","",colnames(cellRanger_AA038))
kallisto_AA038 <- as.SingleCellExperiment(kallisto_AA038)

# Quality control
cellRanger_sce_AA038 <- qc_metrics(cellRanger_AA038, mitochondrial_ens_ids)
kallisto_sce_AA038 <- qc_metrics(kallisto_AA038, mitochondrial_ens_ids)

# EmptyDrops filtering
cellRanger_filt_sce_ed_AA038 <- emptydrops_filt(cellRanger_sce_AA038, lower_ED = 500, EmptyDrops_FDR_thres = 0.01 )
kallisto_filt_sce_ed_AA038 <- emptydrops_filt(kallisto_sce_AA038, lower_ED = 500, EmptyDrops_FDR_thres = 0.01 )

threshold_mito_percentage = 20
low_threshold_cell_counts = 3000
high_threshold_cell_counts = 50000
cells_min_genes_detected_threshold = 700
cellRanger_filt_sce_AA038 <- Filtering_TNBC(cellRanger_filt_sce_ed_AA038,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold, cells_min_threshold = low_threshold_cell_counts)
kallisto_filt_sce_AA038 <- Filtering_TNBC(kallisto_filt_sce_ed_AA038,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold, cells_min_threshold = low_threshold_cell_counts)

doublet_analysis_TNBC <- function(data)
{
  data = red_dim(data)
  data = clustering(data, k=8, hg38_ensembl_gtf)
  colLabels(data)=data$louvain_clusters
  dbl.out <- findDoubletClusters(data)
  print(dbl.out)
  chosen.doublet <- rownames(dbl.out)[isOutlier(dbl.out$num.de, type="lower", log=TRUE)]
  print(chosen.doublet)

  dbl.dens <- computeDoubletDensity(data, d=ncol(reducedDim(data)))
  print(summary(dbl.dens))
  data$DoubletScore <- dbl.dens
  dbl.calls <- doubletThresholding(data.frame(score=dbl.dens),method="griffiths", returnType="call")
  print(summary(dbl.calls))

  data <- scDblFinder(data, clusters=colLabels(data))
  print(table(data$scDblFinder.class))
  return(data)
}

# Caution: This sample seem to have a high doublet content
kallisto_filt_sce_AA038 <- doublet_analysis_TNBC(kallisto_filt_sce_AA038)
cellRanger_filt_sce_AA038 <- doublet_analysis_TNBC(cellRanger_filt_sce_AA038)
# Cluster 7 might be containing a lot of doublets
kallisto_filt_sce_AA038_nodoub <- kallisto_filt_sce_AA038[, (kallisto_filt_sce_AA038$louvain_clusters != 7) &(kallisto_filt_sce_AA038$DoubletScore < 3)  & (kallisto_filt_sce_AA038$scDblFinder.score < 0.25)]
# Cluster 13 might be containing a lot of doublets
cellRanger_filt_sce_AA038_nodoub <- cellRanger_filt_sce_AA038[, (cellRanger_filt_sce_AA038$louvain_clusters != 13) & (cellRanger_filt_sce_AA038$DoubletScore < 3)  & (cellRanger_filt_sce_AA038$scDblFinder.score < 0.25)]

#####################################################################################################################################################################
# 1. Processing: Sample5: AA051
cell_ranger_dir <- "/home/egonie/kike/datos/119_Amaya_scRNAseq_BLANCA/NextSeq2000.RUN148.230213/01.CellRanger/AA051_alldata/outs/raw_feature_bc_matrix"
cellRanger_data <- Read10X(data.dir = cell_ranger_dir, gene.column = 1)
cellRanger_AA051 <- CreateSeuratObject(cellRanger_data, project = "cellRanger")

kallisto_dir <- "/home/egonie/kike/datos/119_Amaya_scRNAseq_BLANCA/NextSeq2000.RUN148.230213/01.Kallisto/output_bus_transcriptome/bustools_results_no_multimappers"
kallisto_data <- BUSpaRse::read_count_output(kallisto_dir, name = "cells_genes_NO_multimapping")
kallisto_AA051 <- CreateSeuratObject(kallisto_data, project = "kallisto")

# Convert to SingleCellExperiment
cellRanger_AA051 <- as.SingleCellExperiment(cellRanger_AA051)
colnames(cellRanger_AA051) <- gsub("\\-.*","",colnames(cellRanger_AA051))
kallisto_AA051 <- as.SingleCellExperiment(kallisto_AA051)

# Quality control
cellRanger_sce_AA051 <- qc_metrics(cellRanger_AA051, mitochondrial_ens_ids)
kallisto_sce_AA051 <- qc_metrics(kallisto_AA051, mitochondrial_ens_ids)

# EmptyDrops filtering
cellRanger_filt_sce_ed_AA051 <- emptydrops_filt(cellRanger_sce_AA051, lower_ED = 500, EmptyDrops_FDR_thres = 0.01 )
kallisto_filt_sce_ed_AA051 <- emptydrops_filt(kallisto_sce_AA051, lower_ED = 500, EmptyDrops_FDR_thres = 0.01 )

threshold_mito_percentage = 20
low_threshold_cell_counts = 3000
high_threshold_cell_counts = 50000
cells_min_genes_detected_threshold = 700
cellRanger_filt_sce_AA051 <- Filtering_TNBC(cellRanger_filt_sce_ed_AA051,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold, cells_min_threshold = low_threshold_cell_counts)
kallisto_filt_sce_AA051 <- Filtering_TNBC(kallisto_filt_sce_ed_AA051,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold, cells_min_threshold = low_threshold_cell_counts)

# uniquifyFeatures
gene_name <- hg38_ensembl_gtf$gene_name[match(rownames(cellRanger_filt_sce_AA051),hg38_ensembl_gtf$gene_id)]
rowData(cellRanger_filt_sce_AA051) <- cbind(ens_id = rownames(cellRanger_filt_sce_AA051),gene_name)
rownames(cellRanger_filt_sce_AA051) <- uniquifyFeatureNames(rownames(cellRanger_filt_sce_AA051), rowData(cellRanger_filt_sce_AA051)$gene_name)

gene_name <- hg38_ensembl_gtf$gene_name[match(rownames(kallisto_filt_sce_AA051),hg38_ensembl_gtf$gene_id)]
rowData(kallisto_filt_sce_AA051) <- cbind(ens_id = rownames(kallisto_filt_sce_AA051),gene_name)
rownames(kallisto_filt_sce_AA051) <- uniquifyFeatureNames(rownames(kallisto_filt_sce_AA051), rowData(kallisto_filt_sce_AA051)$gene_name)
# Remove doublets (there are some in this sample also)
cellRanger_filt_sce_AA051 <- doublet_analysis_TNBC(cellRanger_filt_sce_AA051)
kallisto_filt_sce_AA051 <- doublet_analysis_TNBC(kallisto_filt_sce_AA051)
cellRanger_sce_filt_AA051_nodoubs <- cellRanger_filt_sce_AA051[,cellRanger_filt_sce_AA051$DoubletScore< 3]
kallisto_sce_filt_AA051_nodoubs <- kallisto_filt_sce_AA051[,kallisto_filt_sce_AA051$DoubletScore< 3]

#######################################################################################################################################################################
# save RDS 
setwd("/home/egonie/kike/phd/test_data/paper_figures/figure4")
saveRDS(kallisto_filt_sce_AA017,"kallisto_AA017.rds")
saveRDS(kallisto_filt_sce_AA024,"kallisto_AA024.rds")
saveRDS(kallisto_filt_sce_AA025,"kallisto_AA025.rds")
saveRDS(kallisto_filt_sce_AA038_nodoub,"kallisto_AA038.rds")
saveRDS(kallisto_sce_filt_AA051_nodoubs,"kallisto_AA051.rds")

saveRDS(cellRanger_filt_sce_AA017,"cellRanger_AA017.rds")
saveRDS(cellRanger_filt_sce_AA024,"cellRanger_AA024.rds")
saveRDS(cellRanger_filt_sce_AA025,"cellRanger_AA025.rds")
saveRDS(cellRanger_filt_sce_AA038_nodoub,"cellRanger_AA038.rds")
saveRDS(cellRanger_sce_filt_AA051_nodoubs,"cellRanger_AA051.rds")

#######################################################################################################################################################################
######################################################### Integrate AA017, AA024, AA025, AA038 and AA051 ##############################################################
#######################################################################################################################################################################
kallisto_sce_filt_clus_AA017 = readRDS("kallisto_AA017.rds")
kallisto_sce_filt_clus_AA024 = readRDS("kallisto_AA024.rds")
kallisto_sce_filt_clus_AA025 = readRDS("kallisto_AA025.rds")
kallisto_sce_filt_clus_AA038 = readRDS("kallisto_AA038.rds")
kallisto_sce_filt_clus_AA051 = readRDS("kallisto_AA051.rds")

cellRanger_sce_filt_clus_AA017 = readRDS("cellRanger_AA017.rds")
cellRanger_sce_filt_clus_AA024 = readRDS("cellRanger_AA024.rds")
cellRanger_sce_filt_clus_AA025 = readRDS("cellRanger_AA025.rds")
cellRanger_sce_filt_clus_AA038 = readRDS("cellRanger_AA038.rds")
cellRanger_sce_filt_clus_AA051 = readRDS("cellRanger_AA051.rds")

# UpSet plots gene expression
threshold_minumun_gene_counts = 250
threshold_cells_detected = 25

#upset plot lncRNAs
features = lncrna_names
kallisto_AA017_top_genes <- filter_top_genes(top_genes(kallisto_sce_filt_clus_AA017,threshold_minumun_gene_counts,threshold_cells_detected), features)
kallisto_AA024_top_genes <- filter_top_genes(top_genes(kallisto_sce_filt_clus_AA024,threshold_minumun_gene_counts,threshold_cells_detected), features)
kallisto_AA025_top_genes <- filter_top_genes(top_genes(kallisto_sce_filt_clus_AA025,threshold_minumun_gene_counts,threshold_cells_detected), features)
kallisto_AA038_top_genes <- filter_top_genes(top_genes(kallisto_sce_filt_clus_AA038,threshold_minumun_gene_counts,threshold_cells_detected), features)
kallisto_AA051_top_genes <- filter_top_genes(top_genes(kallisto_sce_filt_clus_AA051,threshold_minumun_gene_counts,threshold_cells_detected), features)

cellRanger_AA017_top_genes <- filter_top_genes(top_genes(cellRanger_sce_filt_clus_AA017,threshold_minumun_gene_counts,threshold_cells_detected), features)
cellRanger_AA024_top_genes <- filter_top_genes(top_genes(cellRanger_sce_filt_clus_AA024,threshold_minumun_gene_counts,threshold_cells_detected), features)
cellRanger_AA025_top_genes <- filter_top_genes(top_genes(cellRanger_sce_filt_clus_AA025,threshold_minumun_gene_counts,threshold_cells_detected), features)
cellRanger_AA038_top_genes <- filter_top_genes(top_genes(cellRanger_sce_filt_clus_AA038,threshold_minumun_gene_counts,threshold_cells_detected), features)
cellRanger_AA051_top_genes <- filter_top_genes(top_genes(cellRanger_sce_filt_clus_AA051,threshold_minumun_gene_counts,threshold_cells_detected), features)

kallisto_rows <- sum(nrow(kallisto_AA017_top_genes),nrow(kallisto_AA024_top_genes), nrow(kallisto_AA025_top_genes), nrow(kallisto_AA038_top_genes), nrow(kallisto_AA051_top_genes))
cellRanger_rows <- sum(nrow(cellRanger_AA017_top_genes),nrow(cellRanger_AA024_top_genes), nrow(cellRanger_AA025_top_genes), nrow(cellRanger_AA038_top_genes), nrow(cellRanger_AA051_top_genes))
kallisto_genes <- as.character(c(rownames(kallisto_AA017_top_genes), rownames(kallisto_AA024_top_genes), rownames(kallisto_AA025_top_genes), rownames(kallisto_AA038_top_genes), rownames(kallisto_AA051_top_genes)))
cellRanger_genes <- as.character(c(rownames(cellRanger_AA017_top_genes), rownames(cellRanger_AA024_top_genes), rownames(cellRanger_AA025_top_genes), rownames(cellRanger_AA038_top_genes), rownames(cellRanger_AA051_top_genes)))
sample_names <- as.character(c(rep("Sample_1", nrow(kallisto_AA017_top_genes)), rep("Sample_2", nrow(kallisto_AA024_top_genes)), rep("Sample_3", nrow(kallisto_AA025_top_genes)), rep("Sample_4", nrow(kallisto_AA038_top_genes)), rep("Sample_5", nrow(kallisto_AA051_top_genes)), rep("Sample_1", nrow(cellRanger_AA017_top_genes)), rep("Sample_2", nrow(cellRanger_AA024_top_genes)), rep("Sample_3", nrow(cellRanger_AA025_top_genes)), rep("Sample_4", nrow(cellRanger_AA038_top_genes)), rep("Sample_5", nrow(cellRanger_AA051_top_genes))))

all_sce_melted <- as.data.frame(cbind(as.character(c(rep("kallisto",kallisto_rows),rep("cellRanger",cellRanger_rows))), c(kallisto_genes,cellRanger_genes), rep(sample_names,2)))
colnames(all_sce_melted) <- c("ident","gene_id", "Sample")
color_palette_iwanthue_5 <- c("#ba4a4f","#56ae6c","#a24f99","#af953c","#6971c9")
comp_pct <- function(count, PANEL, cut) {
  data.frame(count, PANEL, cut) %>% 
    group_by(PANEL, cut) %>% 
    mutate(pct = 100*count / sum(count)) %>% 
    pull(pct)
}

pdf("upset_plot_lncRNAs_250counts_in_25_genes.pdf")
all_sce_melted_tbl=tbl_df(all_sce_melted)
all_sce_melted_tbl %>%
  group_by(gene_id,Sample) %>%
  summarize(ident = list(ident)) %>%
  ggplot(aes(x = ident,fill=Sample)) +
    geom_bar(aes(y = after_stat(comp_pct(count, PANEL, fill))), position="dodge")  +
    scale_x_upset(order_by = "degree", reverse = T) + theme_classic() + scale_fill_manual(values=color_palette_iwanthue_5, name="") + theme(plot.margin = margin(1,0.5,0.5,2, "cm")) + ylab("Percentage of detected lncRNAs for each dataset")
dev.off()

#upset plot protein-coding genes
features = protein_coding_names
kallisto_AA017_top_genes <- filter_top_genes(top_genes(kallisto_sce_filt_clus_AA017,threshold_minumun_gene_counts,threshold_cells_detected), features)
kallisto_AA024_top_genes <- filter_top_genes(top_genes(kallisto_sce_filt_clus_AA024,threshold_minumun_gene_counts,threshold_cells_detected), features)
kallisto_AA025_top_genes <- filter_top_genes(top_genes(kallisto_sce_filt_clus_AA025,threshold_minumun_gene_counts,threshold_cells_detected), features)
kallisto_AA038_top_genes <- filter_top_genes(top_genes(kallisto_sce_filt_clus_AA038,threshold_minumun_gene_counts,threshold_cells_detected), features)
kallisto_AA051_top_genes <- filter_top_genes(top_genes(kallisto_sce_filt_clus_AA051,threshold_minumun_gene_counts,threshold_cells_detected), features)

cellRanger_AA017_top_genes <- filter_top_genes(top_genes(cellRanger_sce_filt_clus_AA017,threshold_minumun_gene_counts,threshold_cells_detected), features)
cellRanger_AA024_top_genes <- filter_top_genes(top_genes(cellRanger_sce_filt_clus_AA024,threshold_minumun_gene_counts,threshold_cells_detected), features)
cellRanger_AA025_top_genes <- filter_top_genes(top_genes(cellRanger_sce_filt_clus_AA025,threshold_minumun_gene_counts,threshold_cells_detected), features)
cellRanger_AA038_top_genes <- filter_top_genes(top_genes(cellRanger_sce_filt_clus_AA038,threshold_minumun_gene_counts,threshold_cells_detected), features)
cellRanger_AA051_top_genes <- filter_top_genes(top_genes(cellRanger_sce_filt_clus_AA051,threshold_minumun_gene_counts,threshold_cells_detected), features)

kallisto_rows <- sum(nrow(kallisto_AA017_top_genes),nrow(kallisto_AA024_top_genes), nrow(kallisto_AA025_top_genes), nrow(kallisto_AA038_top_genes), nrow(kallisto_AA051_top_genes))
cellRanger_rows <- sum(nrow(cellRanger_AA017_top_genes),nrow(cellRanger_AA024_top_genes), nrow(cellRanger_AA025_top_genes), nrow(cellRanger_AA038_top_genes), nrow(cellRanger_AA051_top_genes))
kallisto_genes <- as.character(c(rownames(kallisto_AA017_top_genes), rownames(kallisto_AA024_top_genes), rownames(kallisto_AA025_top_genes), rownames(kallisto_AA038_top_genes), rownames(kallisto_AA051_top_genes)))
cellRanger_genes <- as.character(c(rownames(cellRanger_AA017_top_genes), rownames(cellRanger_AA024_top_genes), rownames(cellRanger_AA025_top_genes), rownames(cellRanger_AA038_top_genes), rownames(cellRanger_AA051_top_genes)))
sample_names <- as.character(c(rep("Sample_1", nrow(kallisto_AA017_top_genes)), rep("Sample_2", nrow(kallisto_AA024_top_genes)), rep("Sample_3", nrow(kallisto_AA025_top_genes)), rep("Sample_4", nrow(kallisto_AA038_top_genes)), rep("Sample_5", nrow(kallisto_AA051_top_genes)), rep("Sample_1", nrow(cellRanger_AA017_top_genes)), rep("Sample_2", nrow(cellRanger_AA024_top_genes)), rep("Sample_3", nrow(cellRanger_AA025_top_genes)), rep("Sample_4", nrow(cellRanger_AA038_top_genes)), rep("Sample_5", nrow(cellRanger_AA051_top_genes))))

all_sce_melted <- as.data.frame(cbind(as.character(c(rep("kallisto",kallisto_rows),rep("cellRanger",cellRanger_rows))), c(kallisto_genes,cellRanger_genes), rep(sample_names,2)))
colnames(all_sce_melted) <- c("ident","gene_id", "Sample")
color_palette_iwanthue_5 <- c("#ba4a4f","#56ae6c","#a24f99","#af953c","#6971c9")
comp_pct <- function(count, PANEL, cut) {
  data.frame(count, PANEL, cut) %>% 
    group_by(PANEL, cut) %>% 
    mutate(pct = 100*count / sum(count)) %>% 
    pull(pct)
}

pdf("upset_plot_protein_coding_250counts_in_25_genes.pdf")
all_sce_melted_tbl=tbl_df(all_sce_melted)
all_sce_melted_tbl %>%
  group_by(gene_id,Sample) %>%
  summarize(ident = list(ident)) %>%
  ggplot(aes(x = ident,fill=Sample)) +
    geom_bar(aes(y = after_stat(comp_pct(count, PANEL, fill))), position="dodge")  +
    scale_x_upset(order_by = "degree", reverse = T) + theme_classic() + scale_fill_manual(values=color_palette_iwanthue_5, name="") + theme(plot.margin = margin(1,0.5,0.5,2, "cm")) + ylab("Percentage of detected protein-coding genes for each dataset")
dev.off()

#########################################################################################################################################################################
#########################################################################################################################################################################
# 2. Integration (following OSCA: https://bioconductor.org/books/3.18/OSCA.multisample/integrating-datasets.html)
# 2.1. Kallisto
# Uncorrected
uncorrected <- uncorrected_integration_TNBC(kallisto_sce_filt_clus_AA017, kallisto_sce_filt_clus_AA024, kallisto_sce_filt_clus_AA025, kallisto_sce_filt_clus_AA038, kallisto_sce_filt_clus_AA051)
# FastMNN
# corrected
mnn.out <- mnn_integration_TNBC(kallisto_sce_filt_clus_AA017, kallisto_sce_filt_clus_AA024, kallisto_sce_filt_clus_AA025, kallisto_sce_filt_clus_AA038, kallisto_sce_filt_clus_AA051)
uncorrected$mnn.out_clustering <- mnn.out$louvain_clusters

# cell type classification
canonical_markers <- c("MS4A1", "CD79A","LILRA4", "CLEC4C", "MKI67", "TOP2A", "PECAM1", "EPCAM", "PDGFRB", "CD68", "GNLY", "NKG7", "JCHAIN", "LEF1", "CD3D")
pdf("canonical_markers.pdf", width = 12)
plotUMAP(mnn.out, colour_by = "louvain_clusters", text_by = "louvain_clusters")
plot_markers_integrated(uncorrected, mnn.out, canonical_markers, "mnn.out_clustering",group = "mnn.out_clustering")
dev.off()

# Subcluster analysis: Cluster 1 seems to be a mix of NK and T cells
subcluster <- uncorrected[,uncorrected$mnn.out_clustering==1]
subcluster <- red_dim(subcluster)
subcluster <- clustering(subcluster, k = 5, hg38_ensembl_gtf) 
mnn.out_subcluster <- mnn.out[,mnn.out$louvain_clusters==1]
pdf("subcluster1_kallisto.pdf")
print(plotUMAP(subcluster, colour_by = "louvain_clusters", point_size = 0.5)  + guides(colour = guide_legend(override.aes = list(size=2))))
plot_markers_integrated(subcluster,mnn.out_subcluster , c("CD3D","GNLY","MKI67","LEF1"), "main_cell_types_markers", group = "louvain_clusters")
dev.off()
# Conclusion: Subcluster 4,7 contain NK cells and subclusters 1,2,3,5,6,8,9 contain T cells
NK_cells <- colnames(subcluster)[(subcluster$louvain_clusters==4) | (subcluster$louvain_clusters==7)]
# annotation
kallisto_general_annotation <- c("T cells", "Epithelial", "T cells", "Plasmablasts",rep("T cells",2),"Plasmablasts","T cells","NK","Plasmablasts","Cycling cells","Myeloid","Plasmablasts","B cells","Mesenchymal","DCs","T cells","Endothelial","Plasmablasts")
map <- data.frame(find=names(table(uncorrected$mnn.out_clustering)),replace=kallisto_general_annotation)
uncorrected$cell_type <- as.character(map[match(uncorrected$mnn.out_clustering, map$find), "replace"])
uncorrected[,NK_cells]$cell_type = "NK"
mnn.out$cell_type <- uncorrected$cell_type


# 2.1. Cell Ranger
# Uncorrected
uncorrected_CR <- uncorrected_integration_TNBC(cellRanger_sce_filt_clus_AA017, cellRanger_sce_filt_clus_AA024, cellRanger_sce_filt_clus_AA025, cellRanger_sce_filt_clus_AA038, cellRanger_sce_filt_clus_AA051)
# FastMNN
# corrected
mnn.out_CR <- mnn_integration_TNBC(cellRanger_sce_filt_clus_AA017, cellRanger_sce_filt_clus_AA024, cellRanger_sce_filt_clus_AA025, cellRanger_sce_filt_clus_AA038, cellRanger_sce_filt_clus_AA051)
uncorrected_CR$mnn.out_clustering <- mnn.out_CR$louvain_clusters

# cell type classification
pdf("canonical_markers_CR.pdf", width = 12)
plotUMAP(mnn.out_CR, colour_by = "louvain_clusters", text_by = "louvain_clusters")
plot_markers_integrated(uncorrected_CR, mnn.out_CR, canonical_markers, "mnn.out_clustering",group = "mnn.out_clustering")
dev.off()

# Subcluster analysis: Cluster 1 seems to be a mix of NK and T cells
subcluster <- uncorrected_CR[,uncorrected_CR$mnn.out_clustering==1]
subcluster <- red_dim(subcluster)
subcluster <- clustering(subcluster, k = 5, hg38_ensembl_gtf) 
mnn.out_subcluster <- mnn.out_CR[,mnn.out_CR$louvain_clusters==1]
pdf("subcluster1_CR.pdf")
print(plotUMAP(subcluster, colour_by = "louvain_clusters", point_size = 0.5)  + guides(colour = guide_legend(override.aes = list(size=2))))
plot_markers_integrated(subcluster,mnn.out_subcluster , c("CD3D","GNLY","MKI67","LEF1"), "main_cell_types_markers", group = "louvain_clusters")
dev.off()
# Conclusion: Subcluster 2,8 contain NK cells and remain contain T cells
NK_cells <- colnames(subcluster)[(subcluster$louvain_clusters==2) | (subcluster$louvain_clusters==8)]

#annotation
cellRanger_general_annotation <- c("T cells", "Plasmablasts" ,rep("T cells",2),"Plasmablasts","Epithelial","T cells","NK",rep("Plasmablasts",2),"T cells","Cycling cells","B cells","Myeloid","Plasmablasts","Mesenchymal","DCs","T cells","Plasmablasts","Endothelial","Plasmablasts")
map <- data.frame(find=names(table(uncorrected_CR$mnn.out_clustering)),replace=cellRanger_general_annotation)
uncorrected_CR$cell_type <- as.character(map[match(uncorrected_CR$mnn.out_clustering, map$find), "replace"])
uncorrected_CR[,NK_cells]$cell_type = "NK"
mnn.out_CR$cell_type <- uncorrected_CR$cell_type

# saveRDS
saveRDS(uncorrected_CR,"cellRanger_uncorrected.rds")
saveRDS(mnn.out_CR,"cellRanger_mnnOUT.rds")
saveRDS(uncorrected,"kallisto_uncorrected.rds")
saveRDS(mnn.out,"kallisto_mnnOUT.rds")

#line416
mnn.out_kallisto_seurat <- to_seurat(uncorrected, mnn.out)
mnn.out_CR_seurat <- to_seurat(uncorrected_CR, mnn.out_CR)
#Main cell types
color_palette_iwanthue_10 <- c("#c8723e","#64aa53","#d8f507","#5d3686","#a2b03d","#c26abb","#5794d7","#cb547e","#b84c3f","#4daf95")
pdf("main_cell_types.pdf")
DimPlot(mnn.out_kallisto_seurat, reduction = "UMAP", label = TRUE, pt.size = 0.3, repel = T,label.size = 5, group.by = "cell_type", cols = color_palette_iwanthue_10)
DimPlot(mnn.out_CR_seurat, reduction = "UMAP", label = TRUE, pt.size = 0.3, repel = T,label.size = 5, group.by = "cell_type", cols = color_palette_iwanthue_10)
dev.off()
