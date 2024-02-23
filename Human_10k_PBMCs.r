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
setwd("/home/egonie/kike/phd/test_data/paper_figures/figure1/PBMCs")

#####################################################################################################################################################################
#####################################################################################################################################################################
# 1. Running metrics
color_palette=c("#ACD1C9","#DCB04A","#98B4D0","#DFDCE5")
no_theme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Mapping rate
pdf("mapping_rate.pdf")
df <- data.frame(pipeline = c("CellRanger", "STARsolo","Kallisto","Salmon"),Mapping_Rate = c(50.7, 51.07, 57.3, 49.35))
ggplot(df, aes(x=pipeline, y=Mapping_Rate, fill=pipeline)) + 
  geom_bar(stat="identity",color="black",show.legend = FALSE) +no_theme + theme(axis.title.y = element_text(size=26, hjust = 0.5, vjust = -1),plot.title = element_text(hjust = 0.5, size = 25, face="bold"),axis.text=element_text(size=17, color="black")) + labs(title="Mapping rate", x="Pipeline", y = "Mapping (%)") + ylim(0, 100)+scale_fill_manual(values=color_palette, name="") + scale_x_discrete(limits=c("CellRanger", "STARsolo", "Kallisto","Salmon"))
dev.off()

# Running time (seconds)
pdf("running_time.pdf")
df <- data.frame(pipeline = c("CellRanger", "STARsolo","Kallisto","Salmon"),Running_time = c(92362, 28618, 2235, 17301))
ggplot(df, aes(x=pipeline, y=Running_time, fill=pipeline)) + 
  geom_bar(stat="identity",color="black", show.legend = FALSE) +no_theme + theme(axis.title.y = element_text(size=26, hjust = 0.5, vjust = -1),plot.title = element_text(hjust = 0.5, size = 25, face="bold"),axis.text=element_text(size=17, color="black")) + labs(title="Running time (s)", x="Pipeline", y = "Running time")+scale_fill_manual(values=color_palette, name="") +scale_x_discrete(limits=c("CellRanger", "STARsolo", "Kallisto","Salmon"))
dev.off()

# Computational_load (Gbs)
pdf("memory_used.pdf")
df <- data.frame(pipeline = c("CellRanger", "STARsolo","Kallisto","Salmon"),mem_used = c(26.77, 28.48, 11.61, 29.16))
ggplot(df, aes(x=pipeline, y=mem_used, fill=pipeline)) + 
  geom_bar(stat="identity",color="black", show.legend = FALSE) +no_theme + theme(axis.title.y = element_text(size=26, hjust = 0.5, vjust = -1),plot.title = element_text(hjust = 0.5, size = 25, face="bold"),axis.text=element_text(size=17, color="black")) + labs(title="Max memory used (Gb)", x="Pipeline", y = "Memory used")+scale_fill_manual(values=color_palette, name="") +scale_x_discrete(limits=c("CellRanger", "STARsolo", "Kallisto","Salmon"))
dev.off()

#####################################################################################################################################################################
#####################################################################################################################################################################
# 2.Processing
cell_ranger_dir <- "/home/egonie/kike/phd/test_data/10X/Parent_NGSC3_DI_PBMC_fastqs/03.CellRanger/Parent_NGSC3_DI_PBMC/outs/raw_feature_bc_matrix"
cellRanger_data <- Read10X(data.dir = cell_ranger_dir, gene.column = 1)
cellRanger <- CreateSeuratObject(cellRanger_data, project = "cellRanger")

STARsolo_dir <- "/home/egonie/kike/phd/test_data/10X/Parent_NGSC3_DI_PBMC_fastqs/03.STARsolo/Solo.out/Gene/raw"
STARsolo_data <- DropletUtils::read10xCounts(samples = STARsolo_dir, col.names = TRUE)
STARsolo <- CreateSeuratObject(counts(STARsolo_data),project = "STARsolo")

kallisto_dir <- "/home/egonie/kike/phd/test_data/10X/Parent_NGSC3_DI_PBMC_fastqs/03.Kallisto/output_bus_transcriptome/bustools_results_no_multimappers"
kallisto_data <- BUSpaRse::read_count_output(kallisto_dir, name = "cells_genes_NO_multimapping")
kallisto <- CreateSeuratObject(kallisto_data, project = "kallisto")

alevin_counts <- "/home/egonie/kike/phd/test_data/10X/Parent_NGSC3_DI_PBMC_fastqs/03.ALEVIN_pseudoalignment/alevin_selective_alignment_output/alevin/quants_mat.gz"
txi <- tximport(alevin_counts, type = "alevin")
alevin <- CreateSeuratObject(counts = txi$counts ,  project = "alevin")

# Convert to SingleCellExperiment
cellRanger <- as.SingleCellExperiment(cellRanger)
colnames(cellRanger) <- gsub("\\-.*","",colnames(cellRanger))
kallisto <- as.SingleCellExperiment(kallisto)
alevin <- as.SingleCellExperiment(alevin)
STARsolo <- as.SingleCellExperiment(STARsolo)

# Quality control
cellRanger_sce <- qc_metrics(cellRanger, mitochondrial_ens_ids)
kallisto_sce <- qc_metrics(kallisto, mitochondrial_ens_ids)
alevin_sce <- qc_metrics(alevin, mitochondrial_ens_ids)
STARsolo_sce <- qc_metrics(STARsolo, mitochondrial_ens_ids)

# EmptyDrops filtering
cellRanger_filt_sce_ed <- emptydrops_filt(cellRanger_sce, lower_ED = 100, EmptyDrops_FDR_thres = 0.0001 )
STARsolo_filt_sce_ed <- emptydrops_filt(STARsolo_sce, lower_ED = 100, EmptyDrops_FDR_thres = 0.0001 )
kallisto_filt_sce_ed <- emptydrops_filt(kallisto_sce, lower_ED = 100, EmptyDrops_FDR_thres = 0.0001 )
# alevin is already filtered (for alevin subset all cells with less than the min of the threshold in cellRanger, STARsolo & Kallisto)
alevin_min_counts_filtering <- min(c(min(cellRanger_filt_sce_ed$sum),min(STARsolo_filt_sce_ed$sum),min(kallisto_filt_sce_ed$sum)))
alevin_filt_sce_ed <- alevin_sce[,alevin_sce$sum >= alevin_min_counts_filtering]

# Doublet analysis
cellRanger_filt_sce_ed_nodoubs <- doublet_analysis(cellRanger_filt_sce_ed)
cellRanger_filt_sce_ed_nodoubs <- cellRanger_filt_sce_ed_nodoubs[,cellRanger_filt_sce_ed_nodoubs$isDoublet == F]
STARsolo_filt_sce_ed_nodoubs <- doublet_analysis(STARsolo_filt_sce_ed)
STARsolo_filt_sce_ed_nodoubs <- STARsolo_filt_sce_ed_nodoubs[,STARsolo_filt_sce_ed_nodoubs$isDoublet == F]
kallisto_filt_sce_ed_nodoubs <- doublet_analysis(kallisto_filt_sce_ed)
kallisto_filt_sce_ed_nodoubs <- kallisto_filt_sce_ed_nodoubs[,kallisto_filt_sce_ed_nodoubs$isDoublet == F]
alevin_filt_sce_ed_nodoubs <- doublet_analysis(alevin_filt_sce_ed)
alevin_filt_sce_ed_nodoubs <- alevin_filt_sce_ed_nodoubs[,alevin_filt_sce_ed_nodoubs$isDoublet == F]

# saveRDS
pbmc_datasets <- list("cellRanger" = cellRanger_filt_sce_ed_nodoubs, "STARsolo" = STARsolo_filt_sce_ed_nodoubs,  "kallisto" = kallisto_filt_sce_ed_nodoubs, "Alevin" = alevin_filt_sce_ed_nodoubs)
saveRDS(pbmc_datasets, file = "/home/egonie/kike/phd/test_data/paper_figures/figure1/pbmc_datasets.rds")
pbmc_datasets <- readRDS("/home/egonie/kike/phd/test_data/paper_figures/figure1/pbmc_datasets.rds")
cellRanger_filt_sce_ed_nodoubs <- pbmc_datasets[["cellRanger"]]
STARsolo_filt_sce_ed_nodoubs <- pbmc_datasets[["STARsolo"]]
kallisto_filt_sce_ed_nodoubs <- pbmc_datasets[["kallisto"]]
alevin_filt_sce_ed_nodoubs <- pbmc_datasets[["Alevin"]] 

#####################################################################################################################################################################
#####################################################################################################################################################################
color_palette=c("#ACD1C9","#DFDCE5","#DCB04A","#98B4D0")

# 3.Quality control filtering
universe <- intersect(rownames(cellRanger_filt_sce_ed_nodoubs), intersect(rownames(STARsolo_filt_sce_ed_nodoubs), intersect(rownames(kallisto_filt_sce_ed_nodoubs), rownames(alevin_filt_sce_ed_nodoubs))))
all_sce <- cbind(cellRanger_filt_sce_ed_nodoubs[universe,], STARsolo_filt_sce_ed_nodoubs[universe,], kallisto_filt_sce_ed_nodoubs[universe,], alevin_filt_sce_ed_nodoubs[universe,])
all_sce$ident <- as.character(all_sce$ident)
all_sce$ident[all_sce$ident=="alevin"] = "Salmon"
all_sce$ident[all_sce$ident=="cellRanger"] = "CellRanger"
all_sce$ident[all_sce$ident=="kallisto"] = "Kallisto"
all_sce$ident <- factor(all_sce$ident, levels = c("CellRanger","STARsolo","Kallisto","Salmon"))

# Mitochondrial content
pdf("QC_Mito_cell.pdf")
all_sce_melted <- as.data.frame(cbind(as.character(all_sce$orig.ident), as.numeric(all_sce$subsets_Mito_percent)))
colnames(all_sce_melted) <- c("ident","Mitochondrial_content")
all_sce_melted[,2] <- as.numeric(all_sce_melted[,2])
all_sce_melted$ident <- factor(all_sce$ident, levels = c("CellRanger","STARsolo","Kallisto","Salmon"))
ggplot(all_sce_melted, aes(x=ident, y=Mitochondrial_content, fill = ident)) + geom_violin(alpha = 1, scale = "width", width = 0.8) + scale_fill_manual(values = color_palette) + theme_classic() + scale_colour_manual(values = color_palette) + ggtitle("Mitochondrial content/cell") + theme(plot.title = element_text(hjust = 0.5, size = 28),axis.text=element_text(size=17), axis.title=element_text(size=17, face = "bold")) + theme(legend.position = "none") + labs(y= "", x="")
dev.off()

# Filtering & normalization
threshold_mito_percentage = 15
#delete also cells with >50.000 counts
high_threshold_cell_counts = 50000
cells_min_genes_detected_threshold = 500
cellRanger_filt_sce <- Filtering(cellRanger_filt_sce_ed_nodoubs,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
STARsolo_filt_sce <- Filtering(STARsolo_filt_sce_ed_nodoubs,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
kallisto_filt_sce <- Filtering(kallisto_filt_sce_ed_nodoubs,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
alevin_filt_sce <- Filtering(alevin_filt_sce_ed_nodoubs, cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)

universe <- intersect(rownames(cellRanger_filt_sce), intersect(rownames(STARsolo_filt_sce), intersect(rownames(kallisto_filt_sce), rownames(alevin_filt_sce))))
all_sce <- cbind(cellRanger_filt_sce[universe,], STARsolo_filt_sce[universe,], kallisto_filt_sce[universe,], alevin_filt_sce[universe,])
all_sce$ident <- as.character(all_sce$ident)
all_sce$ident[all_sce$ident=="alevin"] = "Salmon"
all_sce$ident[all_sce$ident=="cellRanger"] = "CellRanger"
all_sce$ident[all_sce$ident=="kallisto"] = "Kallisto"
all_sce$ident <- factor(all_sce$ident, levels = c("CellRanger","STARsolo","Kallisto","Salmon"))

# UMIs/cell
pdf("QC_UMIs_cell.pdf")
all_sce_melted <- as.data.frame(cbind(as.character(all_sce$orig.ident), as.numeric(all_sce$nCount_RNA)))
colnames(all_sce_melted) <- c("ident","nCount_RNA")
all_sce_melted[,2] <- as.numeric(all_sce_melted[,2])
all_sce_melted$ident <- factor(all_sce$ident, levels = c("CellRanger","STARsolo","Kallisto","Salmon"))
ggplot(all_sce_melted, aes(x=ident, y=nCount_RNA, fill = ident)) + geom_violin(alpha = 1, scale = "width", width = 0.8) + scale_fill_manual(values = color_palette) + theme_classic() + scale_colour_manual(values = color_palette) + ggtitle("UMIs/cell") + theme(plot.title = element_text(hjust = 0.5, size = 28),axis.text=element_text(size=17), axis.title=element_text(size=17, face = "bold")) + theme(legend.position = "none") + labs(y= "", x="") + ylim(0,50000)
dev.off()

#nLncRNAs detected/cell
pdf("QC_detected_lncRNAs_cell.pdf")
all_sce_melted <- as.data.frame(cbind(as.character(all_sce$orig.ident), as.numeric(all_sce$subsets_lncRNA_detected)))
colnames(all_sce_melted) <- c("ident","subsets_lncRNA_detected")
all_sce_melted[,2] <- as.numeric(all_sce_melted[,2])
all_sce_melted$ident <- factor(all_sce$ident, levels = c("CellRanger","STARsolo","Kallisto","Salmon"))
ggplot(all_sce_melted, aes(x=ident, y=subsets_lncRNA_detected, fill = ident)) + geom_violin(alpha = 1, scale = "width", width = 0.8) + scale_fill_manual(values = color_palette) + theme_classic() + scale_colour_manual(values = color_palette) + ggtitle("Detected lncRNAs/cell") + theme(plot.title = element_text(hjust = 0.5, size = 28),axis.text=element_text(size=17), axis.title=element_text(size=17, face = "bold")) + theme(legend.position = "none") + labs(y= "", x="")
dev.off()

#nProtein_coding detected/cell
pdf("QC_detected_protein_coding_cell.pdf")
all_sce_melted <- as.data.frame(cbind(as.character(all_sce$orig.ident), as.numeric(all_sce$subsets_protien_coding_detected)))
colnames(all_sce_melted) <- c("ident","subsets_protien_coding_detected")
all_sce_melted[,2] <- as.numeric(all_sce_melted[,2])
all_sce_melted$ident <- factor(all_sce$ident, levels = c("CellRanger","STARsolo","Kallisto","Salmon"))
ggplot(all_sce_melted, aes(x=ident, y=subsets_protien_coding_detected, fill = ident)) + geom_violin(alpha = 1, scale = "width", width = 0.8) + scale_fill_manual(values = color_palette) + theme_classic() + scale_colour_manual(values = color_palette) + ggtitle("Detected protein coding genes/cell") + theme(plot.title = element_text(hjust = 0.5, size = 28),axis.text=element_text(size=17), axis.title=element_text(size=17, face = "bold")) + theme(legend.position = "none") + labs(y= "", x="")
dev.off()

# Overlap of detected cell per pipeline
input_list <- list(CellRanger = colnames(cellRanger_filt_sce), STARsolo = colnames(STARsolo_filt_sce), Kallisto = colnames(kallisto_filt_sce), Salmon = colnames(alevin_filt_sce))
pdf("upset_plot_cells_test.pdf", width = 8, height = 8)
upset(fromList(input_list), order.by = "degree", mainbar.y.label = "", keep.order = TRUE, text.scale = c(2.5,1.7,1.7,1.6,2.0,1.7), sets.x.label = "Number of cells")
dev.off()

# Gene plots
# Gradient of defined thresholds: 1) 250 counts and present in more than 25 cells, 2) 100 counts and present in more than 10 cells, 3) 50 counts and present in more than 5 cells and 4) 25 counts and present in more than 3 cells
threshold_minumun_gene_counts = 250
threshold_cells_detected = 25

kallisto_top_genes <- top_genes(kallisto_filt_sce,threshold_minumun_gene_counts,threshold_cells_detected )
cellRanger_top_genes <- top_genes(cellRanger_filt_sce,threshold_minumun_gene_counts,threshold_cells_detected )
STARsolo_top_genes <- top_genes(STARsolo_filt_sce,threshold_minumun_gene_counts,threshold_cells_detected )
alevin_top_genes <- top_genes(alevin_filt_sce,threshold_minumun_gene_counts,threshold_cells_detected )

input_list <- list(CellRanger = intersect(rownames(cellRanger_top_genes),lncrna_ens_ids), STARsolo = intersect(rownames(STARsolo_top_genes),lncrna_ens_ids), Kallisto = intersect(rownames(kallisto_top_genes),lncrna_ens_ids), Salmon = intersect(rownames(alevin_top_genes),lncrna_ens_ids))
pdf("upset_plot_lncRNAs_250counts_in_25_genes.pdf", width = 8.75, height = 8)
upset(fromList(input_list), order.by = "degree", mainbar.y.label = "", keep.order = TRUE, text.scale = c(2.5,1.7,1.7,1.6,2.0,1.7), sets.x.label = "Number of lncRNAs",
queries =list(list(query = intersects, params = list("Kallisto"), active = T, color = color_palette[3])))
dev.off()

input_list2 <- list(CellRanger = intersect(rownames(cellRanger_top_genes),protein_coding_ens_ids), STARsolo = intersect(rownames(STARsolo_top_genes),protein_coding_ens_ids), Kallisto = intersect(rownames(kallisto_top_genes),protein_coding_ens_ids), Salmon = intersect(rownames(alevin_top_genes),protein_coding_ens_ids))
pdf("upset_plot_protein_coding_names_250counts_in_25_genes.pdf", width = 8, height = 8)
upset(fromList(input_list2), order.by = "degree", mainbar.y.label = "", keep.order = TRUE, text.scale = c(2.5,1.7,1.7,1.6,2.0,1.7), sets.x.label = "Protein coding genes") 
#grid.text("Protein coding genes intersection",x = 0.65, y=0.97, gp=gpar(fontsize=12, fontface="bold"))
dev.off()

#####################################################################################################################################################################
#####################################################################################################################################################################
# 4. Defined cell types

# Dimensionality reduction and clustering
cellRanger_filt_sce <- red_dim(cellRanger_filt_sce)
cellRanger_sce_filt_clus <- clustering(cellRanger_filt_sce, k = 10, hg38_ensembl_gtf)
STARsolo_filt_sce <- red_dim(STARsolo_filt_sce)
STARsolo_sce_filt_clus <- clustering(STARsolo_filt_sce, k = 10, hg38_ensembl_gtf)
kallisto_filt_sce <- red_dim(kallisto_filt_sce)
kallisto_sce_filt_clus <- clustering(kallisto_filt_sce, k = 10, hg38_ensembl_gtf)
alevin_filt_sce <- red_dim(alevin_filt_sce)
alevin_sce_filt_clus <- clustering(alevin_filt_sce, k = 10, hg38_ensembl_gtf)

pdf("red_dim.pdf")
red_dim_plots(cellRanger_sce_filt_clus)
red_dim_plots(STARsolo_sce_filt_clus)
red_dim_plots(kallisto_sce_filt_clus)
red_dim_plots(alevin_sce_filt_clus)
dev.off()

# Annotate cell types using canonical markers
canonical_markers <- c("MS4A1","LILRA4","CD14","GNLY","LEF1")
pdf("canonical_markers.pdf", width = 10)
plot_markers(cellRanger_sce_filt_clus, canonical_markers, title = "CellRanger", group = "louvain_clusters")
plot_markers(STARsolo_sce_filt_clus, canonical_markers, title = "STARsolo", group = "louvain_clusters")
plot_markers(kallisto_sce_filt_clus, canonical_markers, title = "Kallisto", group = "louvain_clusters")
plot_markers(alevin_sce_filt_clus, canonical_markers, title = "Salmon", group = "louvain_clusters")
dev.off()

cellRanger_general_annotation <- c("B cells","Monocytes","T cells", "DCs", "Monocytes","T cells","T cells","NK", "T cells","Monocytes")
STARsolo_general_annotation <- c("B cells","Monocytes","T cells", "DCs", "Monocytes", "T cells", "T cells","NK", "T cells", "Monocytes","DCs")
kallisto_general_annotation <- c("B cells","Monocytes","T cells", "Monocytes","DCs", "T cells", "T cells","NK","T cells", "Monocytes", "DCs")
alevin_general_annotation <- c("B cells","T cells","Monocytes", "DCs","Monocytes","T cells", "NK","DCs" ,"Monocytes","T cells", "T cells")
pdf("cell_types_annotations.pdf")
cellRanger_sce_filt_clus$cell_type = marker_umaps(cellRanger_sce_filt_clus,clusters = cellRanger_sce_filt_clus$louvain_clusters,cellRanger_general_annotation, title = "cellRanger", final_labels_order = c("B cells","DCs","Monocytes","NK","T cells"))
STARsolo_sce_filt_clus$cell_type = marker_umaps(STARsolo_sce_filt_clus,clusters = STARsolo_sce_filt_clus$louvain_clusters,STARsolo_general_annotation, title = "STARsolo", final_labels_order = c("B cells","DCs","Monocytes","NK","T cells"))
kallisto_sce_filt_clus$cell_type = marker_umaps(kallisto_sce_filt_clus,clusters = kallisto_sce_filt_clus$louvain_clusters,kallisto_general_annotation, title = "kallisto", final_labels_order = c("B cells","DCs","Monocytes","NK","T cells"))
alevin_sce_filt_clus$cell_type = marker_umaps(alevin_sce_filt_clus,clusters = alevin_sce_filt_clus$louvain_clusters,alevin_general_annotation, title = "Salmon", final_labels_order = c("B cells","DCs","Monocytes","NK","T cells"))
dev.off()

pdf("markers_per_celltype.pdf", width = 10)
cellRanger_sce_filt_clus$cell_type <- factor(cellRanger_sce_filt_clus$cell_type, levels = c("B cells","DCs","Monocytes","NK","T cells"))
STARsolo_sce_filt_clus$cell_type <- factor(STARsolo_sce_filt_clus$cell_type, levels = c("B cells","DCs","Monocytes","NK","T cells"))
kallisto_sce_filt_clus$cell_type <- factor(kallisto_sce_filt_clus$cell_type, levels = c("B cells","DCs","Monocytes","NK","T cells"))
alevin_sce_filt_clus$cell_type <- factor(alevin_sce_filt_clus$cell_type, levels = c("B cells","DCs","Monocytes","NK","T cells"))

plot_markers(cellRanger_sce_filt_clus, canonical_markers, title = "cellRanger", group = "cell_type")
plot_markers(STARsolo_sce_filt_clus, canonical_markers, title = "STARsolo", group = "cell_type")
plot_markers(kallisto_sce_filt_clus, canonical_markers, title = "kallisto", group = "cell_type")
plot_markers(alevin_sce_filt_clus, canonical_markers, title = "Salmon", group = "cell_type")
dev.off()

#saveRDS
pbmc_datasets_completed <- list("cellRanger" = cellRanger_sce_filt_clus, "STARsolo" = STARsolo_sce_filt_clus,  "kallisto" = kallisto_sce_filt_clus, "Alevin" = alevin_sce_filt_clus)
saveRDS(pbmc_datasets_completed, file = "/home/egonie/kike/phd/test_data/paper_figures/figure1/pbmc_datasets_updated.rds")
pbmc_datasets_updated <- readRDS("/home/egonie/kike/phd/test_data/paper_figures/figure1/pbmc_datasets_updated.rds")
cellRanger_sce_filt_clus <- pbmc_datasets_updated[["cellRanger"]]
STARsolo_sce_filt_clus <- pbmc_datasets_updated[["STARsolo"]]
kallisto_sce_filt_clus <- pbmc_datasets_updated[["kallisto"]]
alevin_sce_filt_clus <- pbmc_datasets_updated[["Alevin"]] 





