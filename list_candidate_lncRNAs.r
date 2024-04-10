########################################################### Generic #################################################################################################
#####################################################################################################################################################################
libraries <- c("Seurat","patchwork", "tximport", "ggplot2", "scran" , "scater", "SingleCellExperiment", "BUSpaRse", "DropletUtils", "RColorBrewer", "cluster", "UpSetR","scRNAseq", "ggpubr","scDblFinder","reshape","ggupset","dplyr", "ELATUS", "batchelor")

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
# 1.1. 10k
human_10k_pbmc_path <- "/home/egonie/kike/phd/test_data/paper_figures/figure3_characterization_ex_kallisto/pbmc_datasets_updated.rds"
pbmc_datasets_updated <- readRDS(human_10k_pbmc_path)
cellRanger_sce_filt_clus <- pbmc_datasets_updated[["cellRanger"]]
kallisto_sce_filt_clus <- pbmc_datasets_updated[["kallisto"]]
candidates_human_pbmcs_10k <- ELATUS::get_candidates(kallisto_sce_filt_clus, cellRanger_sce_filt_clus, threshold_minumun_gene_counts = 250, threshold_cells_detected = 25, lncrna_names = lncrna_names_human,gtf = hg38_ensembl_gtf, exclusive = T)

# 1.2. 5k
kallisto_PBMCs_5K_ed_filt <- uniquifyFeatures(kallisto_PBMCs_5K_ed_filt, hg38_ensembl_gtf)
kallisto_PBMCs_5K_ed_filt <- red_dim(kallisto_PBMCs_5K_ed_filt)
cellRanger_PBMCs_5K_ed_filt <- uniquifyFeatures(cellRanger_PBMCs_5K_ed_filt, hg38_ensembl_gtf)
cellRanger_PBMCs_5K_ed_filt <- red_dim(cellRanger_PBMCs_5K_ed_filt)
candidates_human_pbmcs_5k <- ELATUS::get_candidates(kallisto_PBMCs_5K_ed_filt, cellRanger_PBMCs_5K_ed_filt, threshold_minumun_gene_counts = 250, threshold_cells_detected = 25, lncrna_names_human, gtf = hg38_ensembl_gtf, exclusive = T)

# Integrate samples
uncorrected_PBMCs_kallisto <- uncorrected_integration_2_samples(kallisto_sce_filt_clus, kallisto_PBMCs_5K_ed_filt, label1="10k",label2="5k")
uncorrected_PBMCs_cellRanger <- uncorrected_integration_2_samples(cellRanger_sce_filt_clus, cellRanger_PBMCs_5K_ed_filt, label1="10k",label2="5k")
#mnn
mnn.out_PBMCs_kallisto <- mnn_integration_2_samples(kallisto_sce_filt_clus, kallisto_PBMCs_5K_ed_filt, label1="10k",label2="5k")
mnn.out_PBMCs_cellRanger <- mnn_integration_2_samples(cellRanger_sce_filt_clus, cellRanger_PBMCs_5K_ed_filt, label1="10k",label2="5k")

pdf("corrected_integration_PBMCs.pdf")
plotUMAP(mnn.out_PBMCs_kallisto, colour_by="identity") + ggtitle("Kallisto: Integrated datasets")
plotUMAP(mnn.out_PBMCs_cellRanger, colour_by="identity")+ ggtitle("CellRanger: Integrated datasets")
plotUMAP(mnn.out_PBMCs_kallisto, text_by = "louvain_clusters",colour_by = "louvain_clusters", point_size = 0.5) + ggtitle("Kallisto: Integrated datasets") + guides(colour = guide_legend(override.aes = list(size=2)))
plotUMAP(mnn.out_PBMCs_cellRanger, text_by = "louvain_clusters",colour_by = "louvain_clusters", point_size = 0.5) + ggtitle("CellRanger: Integrated datasets") + guides(colour = guide_legend(override.aes = list(size=2)))
dev.off()

uncorrected_PBMCs_kallisto$mnn.out_clustering <- mnn.out_PBMCs_kallisto$louvain_clusters
colnames(uncorrected_PBMCs_kallisto) <- colnames(mnn.out_PBMCs_kallisto)
uncorrected_PBMCs_cellRanger$mnn.out_clustering <- mnn.out_PBMCs_cellRanger$louvain_clusters
colnames(uncorrected_PBMCs_cellRanger) <- colnames(mnn.out_PBMCs_cellRanger)

canonical_markers <- c("MS4A1","LILRA4","CD14","GNLY","LEF1")
pdf("Generic_markers_PBMCs.pdf",width = 10)
plot_markers(uncorrected_PBMCs_kallisto, canonical_markers, title = "Kallisto: Generic markers", group = "mnn.out_clustering")
DotPlot(as.Seurat(uncorrected_PBMCs_kallisto),scale.max=100,scale.min=0,dot.min=0,assay = "RNA", features = intersect(canonical_markers,rownames(uncorrected_PBMCs_kallisto)),group.by = c("mnn.out_clustering"),scale.by="size",col.min=0,cols=c("white","darkred")) +geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + theme(axis.text.x = element_text(angle = 45, hjust=1,size=10)) +ggtitle("Kallisto: Generic markers")

plot_markers(uncorrected_PBMCs_cellRanger, canonical_markers, title = "CellRanger: Generic markers", group = "mnn.out_clustering") 
DotPlot(as.Seurat(uncorrected_PBMCs_cellRanger),scale.max=100,scale.min=0,dot.min=0,assay = "RNA", features = intersect(canonical_markers,rownames(uncorrected_PBMCs_cellRanger)),group.by = c("mnn.out_clustering"),scale.by="size",col.min=0,cols=c("white","darkred")) +geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + theme(axis.text.x = element_text(angle = 45, hjust=1,size=10)) +ggtitle("CellRanger: Generic markers")
dev.off()

kallisto_general_annotation <- c("B cells","Monocytes","T cells", "Monocytes","Monocytes","T cells","T cells","NK", "T cells","Monocytes","NK","T cells","B cells", "DCs","T cells")
cellRanger_general_annotation <- c("B cells","Monocytes","T cells", "Monocytes","T cells","T cells","NK", "T cells","Monocytes","Monocytes","NK","T cells","B cells","DCs","B cells","Monocytes")

counts(mnn.out_PBMCs_kallisto) <- counts(uncorrected_PBMCs_kallisto)[rownames(mnn.out_PBMCs_kallisto),]
logcounts(mnn.out_PBMCs_kallisto) <- counts(mnn.out_PBMCs_kallisto)
counts(mnn.out_PBMCs_cellRanger) <- counts(uncorrected_PBMCs_cellRanger)[rownames(mnn.out_PBMCs_cellRanger),]
logcounts(mnn.out_PBMCs_cellRanger) <- counts(mnn.out_PBMCs_cellRanger)
pdf("cell_types_annotations_PBMCs.pdf")
mnn.out_PBMCs_kallisto$cell_type = marker_umaps(mnn.out_PBMCs_kallisto,clusters = mnn.out_PBMCs_kallisto$louvain_clusters,kallisto_general_annotation, title = "PBMCs: Integrated Kallisto")
mnn.out_PBMCs_cellRanger$cell_type = marker_umaps(mnn.out_PBMCs_cellRanger,clusters = mnn.out_PBMCs_cellRanger$louvain_clusters,cellRanger_general_annotation, title = "PBMCs: Integrated Cell Ranger")
dev.off()
uncorrected_PBMCs_kallisto$cell_type <- mnn.out_PBMCs_kallisto$cell_type
uncorrected_PBMCs_cellRanger$cell_type <- mnn.out_PBMCs_cellRanger$cell_type

pbmc_objects <- list("kallisto_uncorrected" = uncorrected_PBMCs_kallisto, "cellRanger_uncorrected" = uncorrected_PBMCs_cellRanger, "kallisto_mnn_out" = mnn.out_PBMCs_kallisto, "cellRanger_mnn_out" = mnn.out_PBMCs_kallisto)
saveRDS(pbmc_objects,"pbmc_objects.rds")
pbmc_objects = readRDS("pbmc_objects.rds")
uncorrected_PBMCs_kallisto <- pbmc_objects[["kallisto_uncorrected"]]
uncorrected_PBMCs_cellRanger <- pbmc_objects[["cellRanger_uncorrected"]]

uncorrected_PBMCs_kallisto$louvain_clusters <- uncorrected_PBMCs_kallisto$mnn.out_clustering
uncorrected_PBMCs_cellRanger$louvain_clusters <- uncorrected_PBMCs_cellRanger$mnn.out_clustering

biologically_relevant_lncRNAs_human_PBMCs <- ELATUS_filtered("Human", uncorrected_PBMCs_kallisto, uncorrected_PBMCs_cellRanger, gene_names = TRUE, threshold_minumun_gene_counts = 250, threshold_cells_detected = 25, dimred_clustering = "PCA", k_neighbors = NA, ratio_threshold = 40, CR_threshold = 10, SI_threshold = 0.15)

# Add cell type information
dff <- table(uncorrected_PBMCs_kallisto$cell_type, uncorrected_PBMCs_kallisto$louvain_clusters)
dff_done = get_index_df(dff)
map <- data.frame(find=colnames(dff), replace = dff_done)
biologically_relevant_lncRNAs_human_PBMCs$main_cell_type <- as.character(map[match(biologically_relevant_lncRNAs_human_PBMCs$cell_type_SI, map$find), "replace"])
saveRDS(biologically_relevant_lncRNAs_human_PBMCs, "biologically_relevant_lncRNAs_human_PBMCs.rds")


#################################################################################################################################
# 2. Human intestine
# 2.1. Intestine 1 (GSM4808339)
kallisto_intestine_pool1_ed_filt <- uniquifyFeatures(kallisto_intestine_pool1_ed_filt, hg38_ensembl_gtf)
kallisto_intestine_pool1_ed_filt <- red_dim(kallisto_intestine_pool1_ed_filt)
cellRanger_intestine_pool1_ed_filt <- uniquifyFeatures(cellRanger_intestine_pool1_ed_filt, hg38_ensembl_gtf)
cellRanger_intestine_pool1_ed_filt <- red_dim(cellRanger_intestine_pool1_ed_filt)
candidates_human_intestine_pool1 <- ELATUS::get_candidates(kallisto_intestine_pool1_ed_filt, cellRanger_intestine_pool1_ed_filt, threshold_minumun_gene_counts = 250, threshold_cells_detected = 25, lncrna_names_human, gtf = hg38_ensembl_gtf, exclusive = T)

# 2.2. Intestine 1 (GSM4808348)
kallisto_intestine_pool2_ed_filt <- uniquifyFeatures(kallisto_intestine_pool2_ed_filt, hg38_ensembl_gtf)
kallisto_intestine_pool2_ed_filt <- red_dim(kallisto_intestine_pool2_ed_filt)
cellRanger_intestine_pool2_ed_filt <- uniquifyFeatures(cellRanger_intestine_pool2_ed_filt, hg38_ensembl_gtf)
cellRanger_intestine_pool2_ed_filt <- red_dim(cellRanger_intestine_pool2_ed_filt)
candidates_human_intestine_pool2 <- ELATUS::get_candidates(kallisto_intestine_pool2_ed_filt, cellRanger_intestine_pool2_ed_filt, threshold_minumun_gene_counts = 250, threshold_cells_detected = 25, lncrna_names_human, gtf = hg38_ensembl_gtf, exclusive = T)

# Integrate samples
uncorrected_intestine_kallisto <- uncorrected_integration_2_samples(kallisto_intestine_pool1_ed_filt, kallisto_intestine_pool2_ed_filt,label1="Pool1",label2="Pool2")
uncorrected_intestine_cellRanger <- uncorrected_integration_2_samples(cellRanger_intestine_pool1_ed_filt, cellRanger_intestine_pool2_ed_filt, label1="Pool1",label2="Pool2")
#mnn
mnn.out_intestine_kallisto <- mnn_integration_2_samples(kallisto_intestine_pool1_ed_filt, kallisto_intestine_pool2_ed_filt,label1="Pool1",label2="Pool2")
mnn.out_intestine_cellRanger <- mnn_integration_2_samples(cellRanger_intestine_pool1_ed_filt, cellRanger_intestine_pool2_ed_filt, label1="Pool1",label2="Pool2")

pdf("corrected_integration_intestine.pdf")
plotUMAP(mnn.out_intestine_kallisto, colour_by="identity") + ggtitle("Kallisto: Integrated datasets")
plotUMAP(mnn.out_intestine_cellRanger, colour_by="identity")+ ggtitle("CellRanger: Integrated datasets")
plotUMAP(mnn.out_intestine_kallisto, text_by = "louvain_clusters",colour_by = "louvain_clusters", point_size = 0.5) + ggtitle("Kallisto: Integrated datasets") + guides(colour = guide_legend(override.aes = list(size=2)))
plotUMAP(mnn.out_intestine_cellRanger, text_by = "louvain_clusters",colour_by = "louvain_clusters", point_size = 0.5) + ggtitle("CellRanger: Integrated datasets") + guides(colour = guide_legend(override.aes = list(size=2)))
dev.off()

uncorrected_intestine_kallisto$mnn.out_clustering <- mnn.out_intestine_kallisto$louvain_clusters
colnames(uncorrected_intestine_kallisto) <- colnames(mnn.out_intestine_kallisto)
uncorrected_intestine_cellRanger$mnn.out_clustering <- mnn.out_intestine_cellRanger$louvain_clusters
colnames(uncorrected_intestine_cellRanger) <- colnames(mnn.out_intestine_cellRanger)

canonical_markers <- c("EPCAM","APOA1","MT1G","FAPB1","PECAM1","RAMP2","CLDN5","CDH5","PDFRA","VCAN","PTPRC","CCL3","HLA-DRA","CD74","KLK11","UPK3B","WT1","DES","CNN1","ACTG2","MYH11","TAGLN","PLP1","PHOX2B","S100B","HAND2")
pdf("Canonical_markers_intestine.pdf",width = 10)
plot_markers(uncorrected_intestine_kallisto, canonical_markers, title = "Kallisto: Generic markers", group = "mnn.out_clustering")
DotPlot(as.Seurat(uncorrected_intestine_kallisto),scale.max=100,scale.min=0,dot.min=0,assay = "RNA", features = intersect(canonical_markers,rownames(uncorrected_intestine_kallisto)),group.by = c("mnn.out_clustering"),scale.by="size",col.min=0,cols=c("white","darkred")) +geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + theme(axis.text.x = element_text(angle = 45, hjust=1,size=10)) +ggtitle("Kallisto: Generic markers")

plot_markers(uncorrected_intestine_cellRanger, canonical_markers, title = "CellRanger: Generic markers", group = "mnn.out_clustering") 
DotPlot(as.Seurat(uncorrected_intestine_cellRanger),scale.max=100,scale.min=0,dot.min=0,assay = "RNA", features = intersect(canonical_markers,rownames(uncorrected_intestine_cellRanger)),group.by = c("mnn.out_clustering"),scale.by="size",col.min=0,cols=c("white","darkred")) +geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + theme(axis.text.x = element_text(angle = 45, hjust=1,size=10)) +ggtitle("CellRanger: Generic markers")
dev.off()

kallisto_general_annotation <- c("Muscularis",rep("Epithelial",4),"Neural","Immune",rep("Epithelial",2),"Neural",rep("Epithelial",3),"Fibroblasts","Mesothelial","Immune", "Fibroblasts","Endothelial")

counts(mnn.out_intestine_kallisto) <- counts(uncorrected_intestine_kallisto)[rownames(mnn.out_intestine_kallisto),]
logcounts(mnn.out_intestine_kallisto) <- counts(mnn.out_intestine_kallisto)
counts(mnn.out_intestine_cellRanger) <- counts(uncorrected_intestine_cellRanger)[rownames(mnn.out_intestine_cellRanger),]
logcounts(mnn.out_intestine_cellRanger) <- counts(mnn.out_intestine_cellRanger)
pdf("cell_types_annotations_intestine.pdf")
mnn.out_intestine_kallisto$cell_type = marker_umaps(mnn.out_intestine_kallisto,clusters = mnn.out_intestine_kallisto$louvain_clusters,kallisto_general_annotation, title = "Intestine: Integrated Kallisto")
dev.off()

uncorrected_intestine_kallisto$cell_type <- mnn.out_intestine_kallisto$cell_type
uncorrected_intestine_cellRanger$cell_type <- mnn.out_intestine_cellRanger$cell_type

intestine_objects <- list("kallisto_uncorrected" = uncorrected_intestine_kallisto, "cellRanger_uncorrected" = uncorrected_intestine_cellRanger, "kallisto_mnn_out" = mnn.out_intestine_kallisto, "cellRanger_mnn_out" = mnn.out_intestine_kallisto)
saveRDS(intestine_objects,"intestine_objects.rds")
intestine_objects = readRDS("intestine_objects.rds")
uncorrected_intestine_kallisto <- intestine_objects[["kallisto_uncorrected"]]
uncorrected_intestine_cellRanger <- intestine_objects[["cellRanger_uncorrected"]]

uncorrected_intestine_kallisto$louvain_clusters <- uncorrected_intestine_kallisto$mnn.out_clustering
uncorrected_intestine_cellRanger$louvain_clusters <- uncorrected_intestine_cellRanger$mnn.out_clustering

biologically_relevant_lncRNAs_human_intestine <- ELATUS_filtered("Human", uncorrected_intestine_kallisto, uncorrected_intestine_cellRanger, gene_names = TRUE, threshold_minumun_gene_counts = 250, threshold_cells_detected = 25, dimred_clustering = "PCA", k_neighbors = NA, ratio_threshold = 40, CR_threshold = 10, SI_threshold = 0.15)

# Add cell type information
dff <- table(uncorrected_intestine_kallisto$cell_type, uncorrected_intestine_kallisto$louvain_clusters)
dff_done = get_index_df(dff)
map <- data.frame(find=colnames(dff), replace = dff_done)
biologically_relevant_lncRNAs_human_PBMCs$main_cell_type <- as.character(map[match(biologically_relevant_lncRNAs_human_PBMCs$cell_type_SI, map$find), "replace"])
saveRDS(biologically_relevant_lncRNAs_human_intestine, "biologically_relevant_lncRNAs_human_intestine.rds")

#################################################################################################################################
# 3. Human healthy lung (GSM5020383)
kallisto_healthy_lung_ed_filt <- uniquifyFeatures(kallisto_healthy_lung_ed_filt, hg38_ensembl_gtf)
kallisto_healthy_lung_ed_filt <- red_dim(kallisto_healthy_lung_ed_filt)
cellRanger_healthy_lung_ed_filt <- uniquifyFeatures(cellRanger_healthy_lung_ed_filt, hg38_ensembl_gtf)
cellRanger_healthy_lung_ed_filt <- red_dim(cellRanger_healthy_lung_ed_filt)
candidates_human_healthy_lung <- ELATUS::get_candidates(kallisto_healthy_lung_ed_filt, cellRanger_healthy_lung_ed_filt, threshold_minumun_gene_counts = 250, threshold_cells_detected = 25, lncrna_names_human, gtf = hg38_ensembl_gtf, exclusive = T)
kallisto_healthy_lung_ed_filt <- clustering(kallisto_healthy_lung_ed_filt, 4, hg38_ensembl_gtf, change_rownames = F)

canonical_markers <- c("EPCAM","CLDN5")
pdf("Canonical_markers_healthy_lung.pdf",width = 10)
plot_markers(kallisto_healthy_lung_ed_filt, canonical_markers, title = "Kallisto: Generic markers", group = "louvain_clusters")
DotPlot(as.Seurat(kallisto_healthy_lung_ed_filt),scale.max=100,scale.min=0,dot.min=0,assay = "RNA", features = intersect(canonical_markers,rownames(kallisto_healthy_lung_ed_filt)),group.by = c("louvain_clusters"),scale.by="size",col.min=0,cols=c("white","darkred")) +geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + theme(axis.text.x = element_text(angle = 45, hjust=1,size=10)) +ggtitle("Kallisto: Generic markers")
dev.off()

kallisto_general_annotation <- c(rep("Endothelial",7),"Unassigned", rep("Endothelial",2))
map <- data.frame(find=names(table(kallisto_healthy_lung_ed_filt$louvain_clusters)),replace=kallisto_general_annotation)
kallisto_healthy_lung_ed_filt$cell_type =  as.character(map[match(kallisto_healthy_lung_ed_filt$louvain_clusters, map$find), "replace"])

biologically_relevant_lncRNAs_human_healthy_lung <- ELATUS_filtered("Human", kallisto_healthy_lung_ed_filt, cellRanger_healthy_lung_ed_filt, gene_names = TRUE, threshold_minumun_gene_counts = 250, threshold_cells_detected = 25, dimred_clustering = "PCA", k_neighbors = NA, ratio_threshold = 40, CR_threshold = 10, SI_threshold = 0.15)

# Add cell type information
dff <- table(kallisto_healthy_lung_ed_filt$cell_type, kallisto_healthy_lung_ed_filt$louvain_clusters)
dff_done = get_index_df(dff)
map <- data.frame(find=colnames(dff), replace = dff_done)
biologically_relevant_lncRNAs_human_healthy_lung$main_cell_type <- as.character(map[match(biologically_relevant_lncRNAs_human_healthy_lung$cell_type_SI, map$find), "replace"])
saveRDS(biologically_relevant_lncRNAs_human_healthy_lung, "biologically_relevant_lncRNAs_human_healthy_lung.rds")

#################################################################################################################################
# 4. Mouse PBMCs
kallisto_PBMCs_mouse_ed_filt <- uniquifyFeatures(kallisto_PBMCs_mouse_ed_filt, gtf_mouse)
kallisto_PBMCs_mouse_ed_filt <- red_dim(kallisto_PBMCs_mouse_ed_filt)
kallisto_PBMCs_mouse_ed_filt <- clustering(kallisto_PBMCs_mouse_ed_filt, 10, gtf_mouse, change_rownames = F)

gene_name <- gtf_mouse$gene_name[match(rownames(cellRanger_PBMCs_mouse_ed_filt),gsub("\\..*","",gtf_mouse$gene_id))]
rowData(cellRanger_PBMCs_mouse_ed_filt) <- cbind(ens_id = rownames(cellRanger_PBMCs_mouse_ed_filt),gene_name)
rownames(cellRanger_PBMCs_mouse_ed_filt) <- uniquifyFeatureNames(rownames(cellRanger_PBMCs_mouse_ed_filt), rowData(cellRanger_PBMCs_mouse_ed_filt)$gene_name)
cellRanger_PBMCs_mouse_ed_filt <- red_dim(cellRanger_PBMCs_mouse_ed_filt)

candidates_PBMCs_mouse <- ELATUS::get_candidates(kallisto_PBMCs_mouse_ed_filt, cellRanger_PBMCs_mouse_ed_filt, threshold_minumun_gene_counts = 250, threshold_cells_detected = 25, lncrna_names_mouse, gtf = gtf_mouse, exclusive = T)

canonical_markers <- c("Ms4a1","Cd79a","Cd14","Nkg7","Lef1","Cd3e")
pdf("Canonical_markers_PBMCs_mouse.pdf",width = 10)
plot_markers(kallisto_PBMCs_mouse_ed_filt, canonical_markers, title = "Kallisto: Generic markers", group = "louvain_clusters")
DotPlot(as.Seurat(kallisto_PBMCs_mouse_ed_filt),scale.max=100,scale.min=0,dot.min=0,assay = "RNA", features = intersect(canonical_markers,rownames(kallisto_PBMCs_mouse_ed_filt)),group.by = c("louvain_clusters"),scale.by="size",col.min=0,cols=c("white","darkred")) +geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) + theme(axis.text.x = element_text(angle = 45, hjust=1,size=10)) +ggtitle("Kallisto: Generic markers")
dev.off()

kallisto_general_annotation <- c("B cells","Monocytes","T cells","Unassigned", "B cells", "T cells","Monocytes","NK","B cells","Monocytes")
map <- data.frame(find=names(table(kallisto_PBMCs_mouse_ed_filt$louvain_clusters)),replace=kallisto_general_annotation)
kallisto_PBMCs_mouse_ed_filt$cell_type =  as.character(map[match(kallisto_PBMCs_mouse_ed_filt$louvain_clusters, map$find), "replace"])

kallisto_PBMCs_mouse_ed_filt <- clustering(kallisto_PBMCs_mouse_ed_filt, 5, gtf_mouse, change_rownames = F)

biologically_relevant_lncRNAs_mouse_PBMCs <- ELATUS_filtered("Mouse", kallisto_PBMCs_mouse_ed_filt, cellRanger_PBMCs_mouse_ed_filt, gene_names = TRUE, threshold_minumun_gene_counts = 250, threshold_cells_detected = 25, dimred_clustering = "PCA", k_neighbors = NA, ratio_threshold = 40, CR_threshold = 10, SI_threshold = 0.15)

# Add cell type information
dff <- table(kallisto_PBMCs_mouse_ed_filt$cell_type, kallisto_PBMCs_mouse_ed_filt$louvain_clusters)
dff_done = get_index_df(dff)
map <- data.frame(find=colnames(dff), replace = dff_done)
biologically_relevant_lncRNAs_mouse_PBMCs$main_cell_type <- as.character(map[match(biologically_relevant_lncRNAs_mouse_PBMCs$cell_type_SI, map$find), "replace"])
saveRDS(biologically_relevant_lncRNAs_mouse_PBMCs, "biologically_relevant_lncRNAs_mouse_PBMCs.rds")

#################################################################################################################################
# 5. Mouse brain 1k cells
mouse_1k_brain_path <- "/home/egonie/kike/phd/test_data/paper_figures/figure3_characterization_ex_kallisto/mouse_datasets_completed.rds"
mouse_brain_datasets_updated <- readRDS(mouse_1k_brain_path)
cellRanger_sce_filt_clus <- mouse_brain_datasets_updated[["cellRanger"]]
kallisto_sce_filt_clus <- mouse_brain_datasets_updated[["kallisto"]]

candidates_mouse_brain <- ELATUS::get_candidates(kallisto_sce_filt_clus, cellRanger_sce_filt_clus, threshold_minumun_gene_counts = 250, threshold_cells_detected = 25, lncrna_names_mouse, gtf = gtf_mouse, exclusive = T)

kallisto_sce_filt_clus <- clustering(kallisto_sce_filt_clus,5, gtf_mouse, change_rownames = F)

biologically_relevant_lncRNAs_mouse_brain <- ELATUS_filtered("Mouse", kallisto_sce_filt_clus, cellRanger_sce_filt_clus, gene_names = TRUE, threshold_minumun_gene_counts = 250, threshold_cells_detected = 25, dimred_clustering = "PCA", k_neighbors = NA, ratio_threshold = 40, CR_threshold = 10, SI_threshold = 0.15)

# Add cell type information
dff <- table(kallisto_sce_filt_clus$cell_type, kallisto_sce_filt_clus$louvain_clusters)
dff_done = get_index_df(dff)
map <- data.frame(find=colnames(dff), replace = dff_done)
biologically_relevant_lncRNAs_mouse_brain$main_cell_type <- as.character(map[match(biologically_relevant_lncRNAs_mouse_brain$cell_type_SI, map$find), "replace"])
saveRDS(biologically_relevant_lncRNAs_mouse_brain, "biologically_relevant_lncRNAs_mouse_brain.rds")

#################################################################################################################################
# 6. Intersect candidates to know how many are tissue-specific
# Individual samples without the TNBC data
input_list <- list(PBMCs_10k = rownames(candidates_human_pbmcs_10k), PBMCs_5k = rownames(candidates_human_pbmcs_5k), Intestine_pool1 = rownames(candidates_human_intestine_pool1) , Intestine_pool2 = rownames(candidates_human_intestine_pool2), Healthy_lung = rownames(candidates_human_healthy_lung))
pdf("upset_plot_candidates_human_individual_samples_test_noTNBC.pdf", width = 10, height = 8)
upset(fromList(input_list), order.by = "degree", mainbar.y.label = "", keep.order = TRUE, text.scale = c(2.5,1.7,1.7,1.6,2.0,1.7), sets.x.label = "Number of lncRNAs") 
dev.off()

#mouse
input_list <- list(PBMCs = rownames(candidates_PBMCs_mouse), Brain = rownames(candidates_mouse_brain))
pdf("upset_plot_candidates_mouse.pdf", width = 8, height = 8)
upset(fromList(input_list), order.by = "degree", mainbar.y.label = "", keep.order = TRUE, text.scale = c(2.5,1.7,1.7,1.6,2.0,1.7), sets.x.label = "Number of lncRNAs") 
dev.off()

#################################################################################################################################
# 7. Human TNBC
kallisto_integrated_objects <- readRDS("/home/egonie/kike/datos/119_Amaya_scRNAseq_BLANCA/integrated_5_samples/kallisto_5_integrated_objects.rds")
uncorrected_kallisto <- kallisto_integrated_objects[[1]]
mnn.out_kallisto <- kallisto_integrated_objects[[2]]
uncorrected_kallisto$louvain_clusters = mnn.out_kallisto$louvain_clusters

cellRanger_integrated_objects <- readRDS("/home/egonie/kike/datos/119_Amaya_scRNAseq_BLANCA/integrated_5_samples/cellRanger_5_integrated_objects.rds")
uncorrected_CR <- cellRanger_integrated_objects[[1]]
mnn.out_CR <- cellRanger_integrated_objects[[2]]

biologically_relevant_lncRNAs_human_TNBC <- ELATUS_filtered("Human", uncorrected_kallisto, uncorrected_CR, gene_names = TRUE, threshold_minumun_gene_counts = 250, threshold_cells_detected = 25, dimred_clustering = "PCA", k_neighbors = NA, ratio_threshold = 40, CR_threshold = 10, SI_threshold = 0.15)

# Add cell type information
dff <- table(uncorrected_PBMCs_kallisto$cell_type, uncorrected_PBMCs_kallisto$louvain_clusters)
dff_done = get_index_df(dff)
map <- data.frame(find=colnames(dff), replace = dff_done)
biologically_relevant_lncRNAs_human_PBMCs$main_cell_type <- as.character(map[match(biologically_relevant_lncRNAs_human_PBMCs$cell_type_SI, map$find), "replace"])
saveRDS(biologically_relevant_lncRNAs_human_PBMCs, "biologically_relevant_lncRNAs_human_PBMCs.rds")





