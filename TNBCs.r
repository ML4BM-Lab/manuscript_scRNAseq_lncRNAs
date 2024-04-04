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

#saveRDS for integration
actual_rds_objects_AA017 <- list("kallisto" = kallisto_sce_filt_clus_AA017, "cellRanger" = cellRanger_sce_filt_clus_AA017)
saveRDS(actual_rds_objects_AA017,"actual_rds_objects_AA017.rds")

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

#saveRDS for integration
actual_rds_objects_AA024 <- list("kallisto" = kallisto_sce_filt_clus_AA024, "cellRanger" = cellRanger_sce_filt_clus_AA024)
saveRDS(actual_rds_objects_AA024,"actual_rds_objects_AA024.rds")


#####################################################################################################################################################################
# 1. Processing: Sample3: AA025
cell_ranger_dir <- "/home/egonie/kike/datos/119_Amaya_scRNAseq_BLANCA/NextSeq2000.RUN81.aacacho/AA025_alldata/outs/raw_feature_bc_matrix"
cellRanger_data <- Read10X(data.dir = cell_ranger_dir, gene.column = 1)
cellRanger_AA025 <- CreateSeuratObject(cellRanger_data, project = "cellRanger")

kallisto_dir <- "/home/egonie/kike/datos/119_Amaya_scRNAseq_BLANCA/NextSeq2000.RUN81.aacacho/output_bus_transcriptome/bustools_results_no_multimappers"
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
cellRanger_filt_sce_ed_AA025 <- emptydrops_filt(cellRanger_sce_AA025, lower_ED = 50, EmptyDrops_FDR_thres = 0.01 )
kallisto_filt_sce_ed_AA025 <- emptydrops_filt(kallisto_sce_AA025, lower_ED = 50, EmptyDrops_FDR_thres = 0.01 )

threshold_mito_percentage = 20
low_threshold_cell_counts = 3000
high_threshold_cell_counts = 50000
cells_min_genes_detected_threshold = 200
cellRanger_filt_sce_AA025 <- Filtering_TNBC(cellRanger_filt_sce_ed_AA025,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold, cells_min_threshold = low_threshold_cell_counts)
kallisto_filt_sce_AA025 <- Filtering_TNBC(kallisto_filt_sce_ed_AA025,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold, cells_min_threshold = low_threshold_cell_counts)

# uniquifyFeatures
gene_name <- hg38_ensembl_gtf$gene_name[match(rownames(cellRanger_filt_sce_AA025),hg38_ensembl_gtf$gene_id)]
rowData(cellRanger_filt_sce_AA025) <- cbind(ens_id = rownames(cellRanger_filt_sce_AA025),gene_name)
rownames(cellRanger_filt_sce_AA025) <- uniquifyFeatureNames(rownames(cellRanger_filt_sce_AA025), rowData(cellRanger_filt_sce_AA025)$gene_name)

gene_name <- hg38_ensembl_gtf$gene_name[match(rownames(kallisto_filt_sce_AA025),hg38_ensembl_gtf$gene_id)]
rowData(kallisto_filt_sce_AA025) <- cbind(ens_id = rownames(kallisto_filt_sce_AA025),gene_name)
rownames(kallisto_filt_sce_AA025) <- uniquifyFeatureNames(rownames(kallisto_filt_sce_AA025), rowData(kallisto_filt_sce_AA025)$gene_name)

#saveRDS for integration
actual_rds_objects_AA025 <- list("kallisto" = kallisto_sce_filt_clus_AA025, "cellRanger" = cellRanger_sce_filt_clus_AA025)
saveRDS(actual_rds_objects_AA025,"actual_rds_objects_AA025.rds")

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

cellRanger_filt_sce_AA038 <- red_dim(cellRanger_filt_sce_AA038)
kallisto_filt_sce_AA038 <- red_dim(kallisto_filt_sce_AA038)
cellRanger_filt_sce_AA038 <- clustering(cellRanger_filt_sce_AA038, k = 10, hg38_ensembl_gtf)
kallisto_filt_sce_AA038 <- clustering(kallisto_filt_sce_AA038, k = 10, hg38_ensembl_gtf)

new_doublet_analysis <- function(data)
{
    dbl.dens <- computeDoubletDensity(data, d=ncol(reducedDim(data)))
    print(summary(dbl.dens))
    data$DoubletScore <- dbl.dens
    dbl.calls <- doubletThresholding(data.frame(score=dbl.dens),method="griffiths", returnType="call")
    data$first_doublet_threshold <- dbl.calls

    data <- scDblFinder(data, clusters=TRUE)
    data$second_doublet_threshold <- data$scDblFinder.class == "doublet"
    return(data)
}
# Caution: This sample seem to have a high doublet content
kallisto_filt_sce_AA038 <- new_doublet_analysis(kallisto_filt_sce_AA038)
cellRanger_filt_sce_AA038 <- new_doublet_analysis(cellRanger_filt_sce_AA038)
kallisto_filt_sce_AA038_nodoub=kallisto_filt_sce_AA038[, (kallisto_filt_sce_AA038$DoubletScore < 3)  & (kallisto_filt_sce_AA038$scDblFinder.score < 0.25)]
cellRanger_filt_sce_AA038_nodoub=cellRanger_filt_sce_AA038[, (cellRanger_filt_sce_AA038$DoubletScore < 3)  & (cellRanger_filt_sce_AA038$scDblFinder.score < 0.25)]

#saveRDS for integration
actual_rds_objects_AA038 <- list("kallisto" = kallisto_filt_sce_AA038_nodoub, "cellRanger" = cellRanger_filt_sce_AA038_nodoub)
saveRDS(actual_rds_objects_AA038,"actual_rds_objects_AA038.rds")

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
cellRanger_filt_sce_AA051 <- red_dim(cellRanger_filt_sce_AA051)
kallisto_filt_sce_AA051 <- red_dim(kallisto_filt_sce_AA051)

dbl.dens <- computeDoubletDensity(cellRanger_filt_sce_AA051, d=ncol(reducedDim(cellRanger_filt_sce_AA051)))
cellRanger_filt_sce_AA051$DoubletScore <- dbl.dens
cellRanger_sce_filt_AA051_nodoubs <- cellRanger_filt_sce_AA051[,cellRanger_filt_sce_AA051$DoubletScore< 3]

dbl.dens <- computeDoubletDensity(kallisto_filt_sce_AA051, d=ncol(reducedDim(kallisto_filt_sce_AA051)))
kallisto_filt_sce_AA051$DoubletScore <- dbl.dens
kallisto_sce_filt_AA051_nodoubs <- kallisto_filt_sce_AA051[,kallisto_filt_sce_AA051$DoubletScore< 3]

#saveRDS for integration
actual_rds_objects_AA051 <- list("kallisto" = kallisto_sce_filt_AA051_nodoubs, "cellRanger" = cellRanger_sce_filt_AA051_nodoubs)
saveRDS(actual_rds_objects_AA051,"actual_rds_objects_AA051.rds")



