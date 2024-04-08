########################################################### Generic #################################################################################################
#####################################################################################################################################################################
libraries <- c("Seurat","patchwork", "tximport", "ggplot2", "scran" , "scater", "SingleCellExperiment", "BUSpaRse", "DropletUtils", "RColorBrewer", "cluster", "UpSetR","scRNAseq", "ggpubr","scDblFinder","reshape", "dplyr", "ggupset", "batchelor","reticulate","argparse","ggbeeswarm","gdata","corrplot","RColorBrewer")

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

uncorrected_CR <- readRDS("cellRanger_uncorrected.rds")
mnn.out_CR <- readRDS("cellRanger_mnnOUT.rds")
uncorrected <- readRDS("kallisto_uncorrected.rds")
mnn.out <- readRDS("kallisto_mnnOUT.rds")

########################################################################################################################################################################
########################################################################################################################################################################
# JIND to annotate subtypes: Using data from https://academic.oup.com/bioinformatics/article/38/9/2488/6543609?login=true & 
# Install it following https://github.com/mohit1997/JIND 

check.integer <- function(x) {
  mean(x == round(x))
}
use_condaenv("jind")
py_config()
parser <- ArgumentParser(description='Run Seurat Classifier')
parser$add_argument('--train_path', default="./data/train.pkl", type="character",
                    help='path to train data frame with labels')
parser$add_argument('--test_path', default="./data/test.pkl", type="character",
                    help='path to test data frame with labels')
parser$add_argument('--column', type="character", default='labels',
                    help='column name for cell types')

args <- parser$parse_args()

jind <- import("jind")

pd <- import("pandas")
lname = args$column

# Use reference annotation from https://www.nature.com/articles/s41588-021-00911-1 (Processed scRNA-seq is available at GSE176078)
f <- Read10X("Wu_etal_2021_BRCA_scRNASeq/", gene.column = 1)
l1 <- read.csv("Wu_etal_2021_BRCA_scRNASeq/metadata.csv")
f_s <- CreateSeuratObject(f, project = "cellRanger")
sce <- as.SingleCellExperiment(f_s)
#get_only_TNBC 
l1 <- l1[l1$subtype=="TNBC",]
sce = sce[,l1$X]
d <- data.frame(t(counts(sce)))

l1_reference <- data.frame(l1$celltype_minor)
colnames(l1_reference) <- c("labels")
rownames(l1_reference) <- l1$X
l1_reference$labels <- as.factor(l1_reference$labels)

k_mat_kallisto <- as.data.frame(t(logcounts(uncorrected)))
k_mat_CellRanger <- as.data.frame(t(logcounts(uncorrected_CR)))

mat1 = d
metadata1 = data.frame(l1_reference)
rownames(metadata1) <- rownames(l1_reference)
colnames(metadata1) <- c("labels")

# kallisto
universe <- intersect(colnames(mat1),colnames(k_mat_kallisto))
mat1 <- mat1[,universe]
mat2 <- k_mat_kallisto[,universe]
obj = jind$JindLib(mat1, as.list(metadata1$labels), path="my_results")
isint = check.integer(mat1[1:100, 1:100]) == 1.

if (isint == TRUE){
  obj$preprocess(count_normalize=TRUE, logt=TRUE)
}

obj$dim_reduction(5000L, 'Var')

train_config = list('val_frac'= 0.2, 'seed' = 0L, 'batch_size'=128L, 'cuda'= FALSE, 'epochs'=15L)
obj$train_classifier(config=train_config, cmat=TRUE)
predicted_label_only_TNBC  = obj$get_filtered_prediction(mat2, frac=0.05, test=FALSE)

#save as .RDS
uncorrected$JIND_raw_predictions <- predicted_label_only_TNBC$raw_predictions
uncorrected$JIND_predictions <- predicted_label_only_TNBC$predictions
mnn.out$JIND_raw_predictions <- predicted_label_only_TNBC$raw_predictions
mnn.out$JIND_predictions <- predicted_label_only_TNBC$predictions

kallisto_integrated_objects <- list("uncorrected" = uncorrected, "mnn.out" = mnn.outj )
saveRDS(kallisto_integrated_objects,"kallisto_5_integrated_objects.rds")

# cellRanger
universe <- intersect(colnames(mat1),colnames(k_mat_CellRanger))
mat1 <- mat1[,universe]
mat2 <- k_mat_kallisto[,universe]
obj = jind$JindLib(mat1, as.list(metadata1$labels), path="my_results")
isint = check.integer(mat1[1:100, 1:100]) == 1.

if (isint == TRUE){
  obj$preprocess(count_normalize=TRUE, logt=TRUE)
}

obj$dim_reduction(5000L, 'Var')

train_config = list('val_frac'= 0.2, 'seed' = 0L, 'batch_size'=128L, 'cuda'= FALSE, 'epochs'=15L)
obj$train_classifier(config=train_config, cmat=TRUE)
predicted_label_only_TNBC  = obj$get_filtered_prediction(mat2, frac=0.05, test=FALSE)

#save as .RDS
uncorrected_CR$JIND_raw_predictions <- predicted_label_only_TNBC$raw_predictions
uncorrected_CR$JIND_predictions <- predicted_label_only_TNBC$predictions
mnn.out_CR$JIND_raw_predictions <- predicted_label_only_TNBC$raw_predictions
mnn.out_CR$JIND_predictions <- predicted_label_only_TNBC$predictions

cellRanger_integrated_objects <- list("uncorrected" = uncorrected_CR, "mnn.out" = mnn.out_CR)
saveRDS(cellRanger_integrated_objects,"cellRanger_5_integrated_objects.rds")

########################################################################################################################################################################
########################################################################################################################################################################
kallisto_integrated_objects <- readRDS("/home/egonie/kike/datos/119_Amaya_scRNAseq_BLANCA/integrated_5_samples/kallisto_5_integrated_objects.rds")
uncorrected_kallisto <- kallisto_integrated_objects[[1]]
mnn.out_kallisto <- kallisto_integrated_objects[[2]]

cellRanger_integrated_objects <- readRDS("/home/egonie/kike/datos/119_Amaya_scRNAseq_BLANCA/integrated_5_samples/cellRanger_5_integrated_objects.rds")
uncorrected_CR <- cellRanger_integrated_objects[[1]]
mnn.out_CR <- cellRanger_integrated_objects[[2]]

uncorrected_kallisto$JIND_simple_predictions <- JIND_simplificate_cell_types(uncorrected_kallisto$JIND_raw_predictions)
uncorrected_CR$JIND_simple_predictions <- JIND_simplificate_cell_types(uncorrected_CR$JIND_raw_predictions)
mnn.out_kallisto$JIND_simple_predictions <- JIND_simplificate_cell_types(uncorrected_kallisto$JIND_raw_predictions)
mnn.out_CR$JIND_simple_predictions <- JIND_simplificate_cell_types(uncorrected_CR$JIND_raw_predictions)

mnn.out_kallisto_seurat <- to_seurat(uncorrected_kallisto, mnn.out_kallisto)
mnn.out_CR_seurat <- to_seurat(uncorrected_CR, mnn.out_CR)

pdf("FeaturePlots_markers.pdf")
markers <- c("CD3E","JCHAIN","GNLY","CD14","PDGFRB","EPCAM","PECAM1","LILRA4","MKI67","KRT19","PDGFRA","MS4A1")
for (i in markers)
{
  print(i)
  print(FeaturePlot(object = mnn.out_kallisto_seurat, features = intersect(i,rownames(mnn.out_kallisto_seurat)), pt.size = 0.25, order = F))
  print(FeaturePlot(object = mnn.out_CR_seurat, features = intersect(i,rownames(mnn.out_kallisto_seurat)), pt.size = 0.25, order = F))
}
dev.off()

pdf("Dotplot_main_markers.pdf")
colnames(uncorrected_CR) <- paste(colnames(uncorrected_CR), uncorrected_CR$batch, sep = "_")
colnames(uncorrected_kallisto) <- paste(colnames(uncorrected_kallisto), uncorrected_kallisto$batch, sep = "_")
DotPlot(as.Seurat(uncorrected_kallisto),scale.max=100,scale.min=0,dot.min=0,assay = "RNA", features = intersect(markers,rownames(uncorrected_kallisto)),group.by = c("JIND_simple_predictions"),scale.by="size",col.min=0,cols=c("white","darkred")) +geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + theme(axis.text.x = element_text(angle = 45, hjust=1,size=10)) +ggtitle("Expression in kallisto")
DotPlot(as.Seurat(uncorrected_CR),scale.max=100,scale.min=0,dot.min=0,assay = "RNA", features = intersect(markers,rownames(uncorrected_kallisto)),group.by = c("JIND_simple_predictions"),scale.by="size",col.min=0,cols=c("white","darkred")) +geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + theme(axis.text.x = element_text(angle = 45, hjust=1,size=10)) +ggtitle("Expression in cellRanger") 
dev.off()

#Sub-cell types
color_palette_iwanthue_15 <- c("#c8723e","#42173f","#305402","#64aa53","#022a57","#9071c7","#5d3686","#a2b03d","#c26abb","#5794d7","#cb547e","#bd5db0","#b84c3f","#05f2b3","#2d9c7e")
pdf("subtypes_kallisto.pdf")
DimPlot(mnn.out_kallisto_seurat, reduction = "UMAP", label = TRUE, pt.size = 0.3, repel = T,label.size = 5, group.by = "JIND_simple_predictions", cols = color_palette_iwanthue_15)
dev.off()
pdf("subtypes_cellRanger.pdf")
DimPlot(mnn.out_CR_seurat, reduction = "UMAP", label = TRUE, pt.size = 0.3, repel = T,label.size = 5, group.by = "JIND_simple_predictions", cols = color_palette_iwanthue_15)
dev.off()

# Cell type specific LncRNAs that are only identified with Kallisto (some of them where found with ELATUS as functionally/biologically plausible lncRNAs)
candidate_lncRNAs <- c("AL136979.1","WT1-AS","LINC02444","AC243960.3","AC022217.2","EGFLAM-AS4","SETBP1-DT","AC009312.1","LINC02065","AL121895.1","AC005972.3")

t=t(logcounts(uncorrected_kallisto[candidate_lncRNAs,]))
t_melted <- as.data.frame(cbind(c(t[,1],t[,2],t[,3],t[,4],t[,5],t[,6],t[,7],t[,8],t[,9],t[,10],t[,11]), c(rep(colnames(t)[1],nrow(t)),rep(colnames(t)[2],nrow(t)),rep(colnames(t)[3],nrow(t)),rep(colnames(t)[4],nrow(t)),rep(colnames(t)[5],nrow(t)),rep(colnames(t)[6],nrow(t)),rep(colnames(t)[7],nrow(t)),rep(colnames(t)[8],nrow(t)),rep(colnames(t)[9],nrow(t)), rep(colnames(t)[10],nrow(t)),rep(colnames(t)[11],nrow(t))),"kallisto"))
colnames(t_melted) <- c("Expression","LncRNA_candidate","ident")
t_melted[,1] <- as.numeric(t_melted[,1])
t_melted$LncRNA_candidate <- factor(t_melted$LncRNA_candidate, levels = candidate_lncRNAs)
t_melted_kallisto <- t_melted

t=t(logcounts(uncorrected_CR[candidate_lncRNAs,]))
t_melted <- as.data.frame(cbind(c(t[,1],t[,2],t[,3],t[,4],t[,5],t[,6],t[,7],t[,8],t[,9],t[,10],t[,11]), c(rep(colnames(t)[1],nrow(t)),rep(colnames(t)[2],nrow(t)),rep(colnames(t)[3],nrow(t)),rep(colnames(t)[4],nrow(t)),rep(colnames(t)[5],nrow(t)),rep(colnames(t)[6],nrow(t)),rep(colnames(t)[7],nrow(t)),rep(colnames(t)[8],nrow(t)),rep(colnames(t)[9],nrow(t)), rep(colnames(t)[10],nrow(t)),rep(colnames(t)[11],nrow(t))),"cellRanger"))
colnames(t_melted) <- c("Expression","LncRNA_candidate","ident")
t_melted[,1] <- as.numeric(t_melted[,1])
t_melted$LncRNA_candidate <- factor(t_melted$LncRNA_candidate, levels = candidate_lncRNAs)

df_melted <- rbind(t_melted_kallisto,t_melted)
color_palette_iwanthue_11 <- c("#b0457b","#92b440","#5a3c90","#5aa554","#c670c3","#45c097","#b6414e","#6d85db","#d19d3a","#be6039","#998c3d")
pdf("violinplots_lncRNAas_candidates.pdf")  #To ease visualization, only plot cells with > 0  expression
p <- ggplot(df_melted[df_melted$ident=="kallisto",], aes(x=LncRNA_candidate, y=Expression, fill = LncRNA_candidate)) + geom_violin(alpha = 1, scale = "width", width = 0.8)  + labs(x="") + theme_classic() + ylim(0,5)+ theme(axis.text.x = element_text(angle = 45, hjust=1,size=10)) 
p2 <- geom_quasirandom(data = df_melted[(df_melted$Expression>0) & (df_melted$ident=="kallisto"),], aes(x=LncRNA_candidate, y=Expression, colour = LncRNA_candidate),width=0.4, groupOnX=TRUE, bandwidth=1,size=0.5,position = position_jitterdodge(seed = 1, dodge.width = 0.8)) 
p + p2 + scale_fill_manual(values = color_palette_iwanthue_11) + scale_colour_manual(values = color_palette_iwanthue_11) +ggtitle("Expression in Kallisto")
p <- ggplot(df_melted[df_melted$ident=="cellRanger",], aes(x=LncRNA_candidate, y=Expression, fill = LncRNA_candidate)) + geom_violin(alpha = 1, scale = "width", width = 0.8)  + labs(x="") + theme_classic() + ylim(0,5)+ theme(axis.text.x = element_text(angle = 45, hjust=1,size=10))
p2 <- geom_quasirandom(data = df_melted[(df_melted$Expression>0) & (df_melted$ident=="cellRanger"),], aes(x=LncRNA_candidate, y=Expression, colour = LncRNA_candidate),width=0.4, groupOnX=TRUE, bandwidth=1,size=0.5,position = position_jitterdodge(seed = 1, dodge.width = 0.8))
p + p2 + scale_fill_manual(values = color_palette_iwanthue_11) + scale_colour_manual(values = color_palette_iwanthue_11) +ggtitle("Expression in cellRanger")
dev.off()

pdf("Dotplots_lncRNAs_candidates.pdf", width = 9)
DotPlot(as.Seurat(uncorrected_kallisto),scale.max=50,scale.min=0,dot.min=0,assay = "RNA", features = intersect(candidate_lncRNAs,rownames(uncorrected_kallisto)),group.by = c("JIND_simple_predictions"),scale.by="size",col.min=0,cols=c("white","darkred")) +geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + theme(axis.text.x = element_text(angle = 45, hjust=1,size=10)) +ggtitle("Expression in kallisto") +  theme(legend.position="bottom", legend.direction="horizontal")
DotPlot(as.Seurat(uncorrected_CR),scale.max=50,scale.min=0,dot.min=0,assay = "RNA", features = intersect(candidate_lncRNAs,rownames(uncorrected_kallisto)),group.by = c("JIND_simple_predictions"),scale.by="size",col.min=0,cols=c("white","darkred")) +geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + theme(axis.text.x = element_text(angle = 45, hjust=1,size=10)) +ggtitle("Expression in cellRanger") + theme(legend.position="bottom", legend.direction="horizontal")
dev.off()

##################################################################################################################################################################
##################################################################################################################################################################
# Correlate expression of AL121895.1 and EPB41L1 in any cell type. 
# Since scRNA-seq is very expression, we are not considering cells that are not expressed in both AL121895.1 and EPB41L1 to avoid incresing the correlation due to this high sparsity
AL121895_expression = matrix(,)
EPB41L1_expression = matrix(,)
for (ct in names(table(uncorrected_kallisto$JIND_simple_predictions)))
{
  print(paste("kallisto",ct))
  A = as.numeric(logcounts(uncorrected_kallisto["AL121895.1",uncorrected_kallisto$JIND_simple_predictions==ct]))
  B = as.numeric(logcounts(uncorrected_kallisto["EPB41L1",uncorrected_kallisto$JIND_simple_predictions==ct],))
  my_data = as.data.frame(cbind("AL121895.1"=A[!(A == 0 & B == 0)],"EPB41L1"=B[!(A ==0 & B == 0)]))
  AL121895_expression <- cbindX(AL121895_expression,as.data.frame(my_data[,1]))
  EPB41L1_expression <- cbindX(EPB41L1_expression,as.data.frame(my_data[,2]))
}
AL121895_expression <- AL121895_expression[,-1]
colnames(AL121895_expression) <- names(table(uncorrected_kallisto$JIND_simple_predictions))
EPB41L1_expression <- EPB41L1_expression[,-1]
colnames(EPB41L1_expression) <- names(table(uncorrected_kallisto$JIND_simple_predictions))

d_cor <- cor(as.data.frame(AL121895_expression), as.data.frame(EPB41L1_expression),method = "spearman",use="complete.obs")
diag(d_cor) <- correlations
for (i in 1:nrow(d_cor))
{
  for (j in 1:ncol(d_cor))
  {
    if (i==j)
    {
    }
    else
    {
      d_cor[i,j]=0
    }
  }
}
pdf("correlations_AL121895.1_EPB41L1.pdf")
corrplot(as.matrix(d_cor),col = COL2('RdBu', 10),tl.col = 'black')
dev.off()

