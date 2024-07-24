########################################################### Generic #################################################################################################
#####################################################################################################################################################################
libraries <- c("Seurat","patchwork", "tximport", "ggplot2", "scran" , "scater", "SingleCellExperiment", "BUSpaRse", "DropletUtils", "RColorBrewer", "cluster", "UpSetR","scRNAseq", "ggpubr","scDblFinder","reshape","Signac")

lapply(libraries, require, character.only = TRUE)
# Import gtf annotation & setWD, source functions
gencode_path <- "~/dato-activo/reference.genomes_kike/GRCh38/gencode/gencode.v37.annotation.gtf"
gtf <- rtracklayer::import(gencode_path)
gtf$gene_biotype <- gtf$gene_type

mitochondrial_ens_ids <- unique(gtf$gene_id[grep("^MT-",gtf$gene_name)])
lncrna_ens_ids <- unique(c(gtf$gene_id[grep("lncRNA",gtf$gene_type)]))
protein_coding_ens_ids <- unique(c(gtf$gene_id[gtf$gene_type=="protein_coding"]))

lncrna_names <- unique(gtf$gene_name[gtf$gene_id %in% lncrna_ens_ids])
protein_coding_names <-  unique(gtf$gene_name[gtf$gene_id %in% protein_coding_ens_ids])

source("/home/egonie/kike/phd/git_rep_hpclogin/manuscript_scRNAseq_lncRNAs/manuscript_functions.r")
setwd("/home/egonie/kike/phd/test_data/10X/MULTIOME_PBMCS_3K")

#####################################################################################################################################################################
#####################################################################################################################################################################

# 1. Preprocessing
# 1.1. scATAC-seq & Filter low-quality cells & Generate gene activity matrix
counts <- Read10X_h5("/home/egonie/kike/phd/test_data/10X/MULTIOME_PBMCS_3K/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_raw_feature_bc_matrix.h5")
fragpath <- "/home/egonie/kike/phd/test_data/10X/MULTIOME_PBMCS_3K/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz"
pbmc <- CreateSeuratObject(counts = counts$`Peaks`,assay = "ATAC")
pbmc[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks,sep = c(":", "-"),fragments = fragpath,annotation = gtf)
DefaultAssay(pbmc) <- "ATAC"
pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)

pbmc_all <- pbmc
pbmc <- subset(x = pbmc,subset = nCount_ATAC < 100000 & nCount_ATAC > 1000 &nucleosome_signal < 2 &TSS.enrichment > 1)
#QC plots ATAC
pdf("/home/egonie/kike/phd/test_data/paper_figures/figure3/QC_ATAC.pdf")
qc_atac_data <- as.data.frame(cbind(pbmc$nCount_ATAC, pbmc$TSS.enrichment, pbmc$nucleosome_signal))
color_ATAC="#a09344"
colnames(qc_atac_data) <- c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal")
VlnPlot(object = pbmc,features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),ncol = 3,pt.size = 1)
VlnPlot(object = pbmc,features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),ncol = 3,pt.size = 0)
VlnPlot(object = pbmc_all,features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),ncol = 3,pt.size = 0)
#ncounts
ggplot(qc_atac_data, aes(x="",y=nCount_ATAC, fill = color_ATAC)) + geom_violin(alpha = 1, scale = "width", width = 0.8) + scale_fill_manual(values = color_ATAC) + theme_classic() + scale_colour_manual(values = color_ATAC) + theme_classic() + ggtitle("Counts ATAC/cell") + theme(plot.title = element_text(hjust = 0.5, size = 28),axis.text=element_text(size=17), axis.title=element_text(size=17, face = "bold")) + theme(legend.position = "none") + labs(y= "", x="ATAC")
#TSS enrichment: TSS enrichment score as the number of Tn5 integration site around the TSS normalized to the number of Tn5 integration sites in flanking regions
ggplot(qc_atac_data, aes(x="",y=TSS.enrichment, fill = color_ATAC)) + geom_violin(alpha = 1, scale = "width", width = 0.8) + scale_fill_manual(values = color_ATAC) + theme_classic() + scale_colour_manual(values = color_ATAC) + theme_classic() + ggtitle("TSS enrichment/cell") + theme(plot.title = element_text(hjust = 0.5, size = 28),axis.text=element_text(size=17), axis.title=element_text(size=17, face = "bold")) + theme(legend.position = "none") + labs(y= "", x="ATAC")
#Nucleosome signal: approximate ratio of mononucleosomal to nucleosome-free fragments 
ggplot(qc_atac_data, aes(x="",y=nucleosome_signal, fill = color_ATAC)) + geom_violin(alpha = 1, scale = "width", width = 0.8) + scale_fill_manual(values = color_ATAC) + theme_classic() + scale_colour_manual(values = color_ATAC) + theme_classic() + ggtitle("Nucleosome signal/cell") + theme(plot.title = element_text(hjust = 0.5, size = 28),axis.text=element_text(size=17), axis.title=element_text(size=17, face = "bold")) + theme(legend.position = "none") + labs(y= "", x="ATAC")
dev.off()

# Normalization
pbmc <- RunTFIDF(pbmc) 
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
# The first LSI component captures sequencing depth (technical variation) rather than biological variation. The component should be removed from downstream analysis and we perform downstream steps without this component
# Reddim
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)

# Build gene activities considering all genes (https://rdrr.io/cran/Signac/man/GeneActivity.html)
gene.activities <- GeneActivity(pbmc, biotypes = NULL)
colnames(gene.activities) <- gsub("-.*","",(gsub("-.*","",colnames(gene.activities))))

# 1.2. GEX
# 1.2.1. Cell Ranger
cellRanger_data <- Read10X(data.dir = "/home/egonie/kike/phd/test_data/10X/MULTIOME_PBMCS_3K/01.CellRanger/pbmc_granulocyte_sorted_3k_alldata_cellRanger_modified_GEX_test2/outs/raw_feature_bc_matrix", gene.column = 1)
colnames(cellRanger_data) <- gsub("-.*","",colnames(cellRanger_data))
cellRanger <- CreateSeuratObject(cellRanger_data, project = "cellRanger")
cellRanger <- as.SingleCellExperiment(cellRanger)
gene_name <- gtf$gene_name[match(rownames(cellRanger),gtf$gene_id)]
rowData(cellRanger) <- cbind(ens_id = rownames(cellRanger),gene_name)
rownames(cellRanger) <- uniquifyFeatureNames(rownames(cellRanger), rowData(cellRanger)$gene_name)

# 1.2.2. Kallisto
kallisto <- kallisto_processing(kallisto_dir = "/home/egonie/kike/phd/test_data/10X/MULTIOME_PBMCS_3K/01.Kallisto/output_bus_transcriptome_introns/bustools_results_no_multimappers", name = "cells_genes_NO_multimapping", mito_gene = mitochondrial_ens_ids)
gene_name <- gtf$gene_name[match(rownames(kallisto),gtf$gene_id)]
rowData(kallisto) <- cbind(ens_id = rownames(kallisto),gene_name)
rownames(kallisto) <- uniquifyFeatureNames(rownames(kallisto), rowData(kallisto)$gene_name)

# get same universe of cells and genes
universe_cells <- intersect(colnames(cellRanger),intersect(colnames(gene.activities), colnames(kallisto)))
kallisto_GEX_filtered <- kallisto[,universe_cells]
cellRanger_GEX_filtered <- cellRanger[,universe_cells]
gene.activities_filtered <- gene.activities[,universe_cells]

# Plot and filter GEX cells with high mitocondrial content (minimun 1000 counts and <20% of mito content). Then filter also the atac
cellRanger_GEX_filtered <- qc_metrics(pipeline = cellRanger_GEX_filtered, mitochondrial_ens_ids = mitochondrial_names)  
#QC RNA-seq (differences between Kallisto and cellRanger)
all_sce <- as.data.frame(cbind(c(cellRanger_GEX_filtered$subsets_Mito_percent, kallisto_GEX_filtered$subsets_Mito_percent),c(rep("cellRanger",ncol(cellRanger_GEX_filtered)),rep("kallisto",ncol(kallisto_GEX_filtered)))))
colnames(all_sce) <- c("Mitochondrial_content","ident")
all_sce$ident <- as.character(all_sce$ident)
all_sce$ident <- factor(all_sce$ident, levels = c("cellRanger","kallisto"))
all_sce[,1] <- as.numeric(all_sce[,1])
color_palette_iwanthue=c("#9473c6","#cc546d")
dodge <- position_dodge(width = 0.85)
pdf("/home/egonie/kike/phd/test_data/paper_figures/figure3/QC_RNA_Mito_cell.pdf")
ggplot(all_sce, aes(x=ident, y=Mitochondrial_content, fill = ident)) + geom_violin(alpha = 1, scale = "width", width = 0.8, position = dodge) + scale_fill_manual(values = color_palette_iwanthue) + theme_classic() + scale_colour_manual(values = color_palette_iwanthue) + ggtitle("Mitochondrial content/cell") + theme(plot.title = element_text(hjust = 0.5, size = 20),axis.text=element_text(size=12), axis.title=element_text(size=14, face = "bold"))  + labs(y= "", x="") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

# Delete nuclei with high mitochondrial content and very low expression
threshold_mito_percentage = 25
low_threshold_cell_counts = 500
cellRanger_GEX_filtered <- cellRanger_GEX_filtered[,(cellRanger_GEX_filtered$sum > low_threshold_cell_counts) & (cellRanger_GEX_filtered$subsets_Mito_percent< threshold_mito_percentage)]
kallisto_GEX_filtered <- kallisto_GEX_filtered[,(kallisto_GEX_filtered$sum > low_threshold_cell_counts) & (kallisto_GEX_filtered$subsets_Mito_percent< threshold_mito_percentage)]

# Get common cells and genes to be able to compare matrices
universe_cells <- intersect(colnames(cellRanger_GEX_filtered),intersect(colnames(gene.activities_filtered), colnames(kallisto_GEX_filtered)))
universe_genes <- intersect(rownames(cellRanger_GEX_filtered),intersect(rownames(gene.activities_filtered), rownames(kallisto_GEX_filtered)))
kallisto_GEX_filtered <- kallisto_GEX_filtered[universe_genes,universe_cells]
cellRanger_GEX_filtered <- cellRanger_GEX_filtered[universe_genes,universe_cells]
gene.activities_filtered <- gene.activities_filtered[universe_genes,universe_cells]
#QC RNA-seq (differences between Kallisto and cellRanger)
all_sce <- as.data.frame(cbind(c(cellRanger_GEX_filtered$nCount_RNA, kallisto_GEX_filtered$nCount_RNA),c(rep("cellRanger",ncol(cellRanger_GEX_filtered)),rep("kallisto",ncol(kallisto_GEX_filtered)))))
colnames(all_sce) <- c("nCount_RNA","ident")
all_sce$ident <- as.character(all_sce$ident)
all_sce$ident <- factor(all_sce$ident, levels = c("cellRanger","kallisto"))
all_sce[,1] <- as.numeric(all_sce[,1])
pdf("/home/egonie/kike/phd/test_data/paper_figures/figure3/QC_RNA_nCount_RNA_cell.pdf")
ggplot(all_sce, aes(x=ident, y=nCount_RNA, fill = ident)) + geom_violin(alpha = 1, scale = "width", width = 0.8, position = dodge) + scale_fill_manual(values = color_palette_iwanthue) + theme_classic() + scale_colour_manual(values = color_palette_iwanthue) + ggtitle("UMIs/cell") + theme(plot.title = element_text(hjust = 0.5, size = 20),axis.text=element_text(size=12), axis.title=element_text(size=14, face = "bold"))  + labs(y= "", x="") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

#save intermediate RDS objects
objects_list <- list(kallisto_GEX = kallisto_GEX_filtered, cellRanger_GEX = cellRanger_GEX_filtered, gene_activities_ATAC = gene.activities_filtered)
saveRDS(objects_list, "preprocessed_objects.rds")
# load rds
objects_list <- readRDS("preprocessed_objects.rds")
kallisto_GEX <- objects_list[["kallisto_GEX"]]
cellRanger_GEX <- objects_list[["cellRanger_GEX"]]
gene.activities <- objects_list[["gene_activities_ATAC"]]

# 2. Dimensionality reduction (WNN)
cellRanger_GEX <- logNormCounts(cellRanger_GEX)
kallisto_GEX <- logNormCounts(kallisto_GEX)
cellRanger_GEX <- red_dim(cellRanger_GEX)
kallisto_GEX <- red_dim(kallisto_GEX)

# Merge in a single object for WNN
r= CreateChromatinAssay(counts = counts$Peaks[,paste(colnames(cellRanger_GEX),'-1',sep="")],sep = c(":", "-"),fragments = fragpath,annotation = gtf)
cellRanger_GEX_backup <- cellRanger_GEX
colnames(cellRanger_GEX_backup) <- paste(colnames(cellRanger_GEX),'-1',sep="")
pbmc_CR = as.Seurat(cellRanger_GEX_backup)
pbmc_CR[['ATAC']] <- r
DefaultAssay(pbmc_CR) <- "ATAC"

kallisto_GEX_backup <- kallisto_GEX
colnames(kallisto_GEX_backup) <- paste(colnames(kallisto_GEX),'-1',sep="")
pbmc_kallisto = as.Seurat(kallisto_GEX_backup)
pbmc_kallisto[['ATAC']] <- r
DefaultAssay(pbmc_kallisto) <- "ATAC"

pbmc_CR_sce <- wnn_clustering(pbmc_CR)
pbmc_kallisto_sce <- wnn_clustering(pbmc_kallisto)

# Annotate cell types using canonical markers
canonical_markers <- c("MS4A1","LILRA4","CD14","GNLY","NKG7","LEF1","CD3D")
pdf("/home/egonie/kike/phd/test_data/paper_figures/figure3/canonical_markers.pdf", width = 10)
plot(DimPlot(pbmc_CR_sce[[1]], reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 4, repel = TRUE) + ggtitle("WNN: CellRanger"))
plot_markers(pbmc_CR_sce[[2]], canonical_markers, title = "CellRanger",group = "seurat_clusters")
plot(DimPlot(pbmc_kallisto_sce[[1]], reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 4, repel = TRUE) + ggtitle("WNN: Kallisto"))
plot_markers(pbmc_kallisto_sce[[2]], canonical_markers, title = "Kallisto", group = "seurat_clusters")
dev.off()

cellRanger_general_annotation <- c("Monocytes",rep("T cells",4),"NKT","Monocytes","NK",rep("B cells",2),"Monocytes","T cells","DCs")
kallisto_general_annotation <- c("Monocytes",rep("T cells",4),"NKT","Monocytes","NK","T cells",rep("B cells",2),"Monocytes","DCs")
cellRanger_GEX_clus <- pbmc_CR_sce[[2]]
kallisto_GEX_clus <- pbmc_kallisto_sce[[2]]
pdf("/home/egonie/kike/phd/test_data/paper_figures/figure3/scRNAseq_celltypes_WNN_clusters.pdf")
cellRanger_GEX_clus$cell_type_WNN_clusters = marker_umaps(cellRanger_GEX_clus,clusters = cellRanger_GEX_clus$seurat_clusters,cellRanger_general_annotation, title = "cellRanger WNN GEX & ATAC",reduction = "WNN.UMAP", final_labels_order = c("B cells","DCs","Monocytes","NK","NKT","T cells"))
kallisto_GEX_clus$cell_type_WNN_clusters = marker_umaps(kallisto_GEX_clus,clusters = kallisto_GEX_clus$seurat_clusters,kallisto_general_annotation, title = "kallisto WNN GEX & ATAC",reduction = "WNN.UMAP" , final_labels_order = c("B cells","DCs","Monocytes","NK","NKT","T cells"))
dev.off()

#save intermediate RDS objects
objects_list <- list(kallisto_GEX = kallisto_GEX_clus, cellRanger_GEX = cellRanger_GEX_clus, gene_activities_ATAC = gene.activities)
saveRDS(objects_list, "preprocessed_WNN_ANNOTATED_objects.rds")
# load rds
objects_list <- readRDS("preprocessed_WNN_ANNOTATED_objects.rds")
kallisto_GEX <- objects_list[["kallisto_GEX"]]
cellRanger_GEX <- objects_list[["cellRanger_GEX"]]
gene.activities <- objects_list[["gene_activities_ATAC"]]

pdf("/home/egonie/kike/phd/test_data/paper_figures/figure3/markers_figure3_PBMCs.pdf", width = 10)
cellRanger_GEX$cell_type <- factor(cellRanger_GEX$cell_type_WNN_clusters, levels = c("B cells","DCs","Monocytes","NK","NKT","T cells"))
kallisto_GEX$cell_type <- factor(kallisto_GEX$cell_type_WNN_clusters, levels = c("B cells","DCs","Monocytes","NK","NKT","T cells"))
plot_markers(cellRanger_GEX, canonical_markers, title = "Cell Ranger", group = "cell_type")
plot_markers(kallisto_GEX, canonical_markers, title = "Kallisto", group = "cell_type")
dev.off()

# Gene expression detected by Cell Ranger & Kallisto
#nLncRNAs detected/cell
detected_lncRNAs_nuclei_CellRanger <- colSums(counts(cellRanger_GEX[intersect(rownames(cellRanger_GEX),lncrna_names)],)!=0)
detected_lncRNAs_nuclei_Kallisto <- colSums(counts(kallisto_GEX[intersect(rownames(kallisto_GEX),lncrna_names)],)!=0)

all_sce <- as.data.frame(cbind(c(detected_lncRNAs_nuclei_CellRanger, detected_lncRNAs_nuclei_Kallisto),c(rep("cellRanger",length(detected_lncRNAs_nuclei_CellRanger)),rep("kallisto",length(detected_lncRNAs_nuclei_Kallisto)))))
colnames(all_sce) <- c("Detected_lncRNAs","ident")
all_sce$ident <- as.character(all_sce$ident)
all_sce$ident <- factor(all_sce$ident, levels = c("cellRanger","kallisto"))
all_sce[,1] <- as.numeric(all_sce[,1])
color_palette_iwanthue=c("#9473c6","#cc546d")

pdf("/home/egonie/kike/phd/test_data/paper_figures/figure3/Detected_lncRNAs_nuclei.pdf")
ggplot(all_sce, aes(x=ident, y=Detected_lncRNAs, fill = ident)) + geom_violin(alpha = 1, scale = "width", width = 0.8) + scale_fill_manual(values = color_palette_iwanthue) + theme_classic() + scale_colour_manual(values = color_palette_iwanthue) + ggtitle("Detected lncRNAs/cell") + theme(plot.title = element_text(hjust = 0.5, size = 20),axis.text=element_text(size=12), axis.title=element_text(size=14, face = "bold"))  + labs(y= "", x="") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

threshold_minumun_gene_counts = 250
threshold_cells_detected = 25
kallisto_top_genes <- top_genes(kallisto_GEX,threshold_minumun_gene_counts,threshold_cells_detected )
cellRanger_top_genes <- top_genes(cellRanger_GEX,threshold_minumun_gene_counts,threshold_cells_detected )
input_list <- list(CellRanger = intersect(rownames(cellRanger_top_genes),lncrna_names), Kallisto = intersect(rownames(kallisto_top_genes),lncrna_names))
pdf("/home/egonie/kike/phd/test_data/paper_figures/figure3/highly_expressed_lncRNAs.pdf")
print(upset(fromList(input_list), order.by = "degree", mainbar.y.label = "Highly-expressed lncRNAs", keep.order = TRUE, text.scale = c(2.5,1.7,1.7,1.6,2.0,1.7),, sets.x.label = "Number of lncRNAs"))
df <- data.frame(Pipeline=c("Cell Ranger","Kallisto"),exclusive_lncRNAs=c(length(input_list$CellRanger), length(input_list$Kallisto)))
p1 <- ggplot(df,aes(x=Pipeline, y=exclusive_lncRNAs)) +geom_bar(stat="identity", width=0.65) +xlab("Pipeline") + ylab("Highly-expressed lncRNAs") + theme_classic() + theme(axis.text = element_text(size=13), axis.title = element_text(size=14)) + ggtitle("Highly-expressed lncRNAs") + theme(plot.title = element_text(hjust = 0.5, size = 15))
print(p1)
dev.off()

####################################################################################################################################
###################################### Compare which GEX is more similar to ATAC-seq ###############################################
####################################################################################################################################

# test that order is the same
test_match_order(paste(colnames(gene.activities),"-1",sep=""), colnames(cellRanger_GEX))
test_match_order(paste(colnames(gene.activities),"-1",sep=""), colnames(kallisto_GEX))
test_match_order(rownames(gene.activities), rownames(cellRanger_GEX))
test_match_order(rownames(gene.activities), rownames(kallisto_GEX))

# Metric to compare: For each nuclei, count the number of genes that are simultaneously activatated both in the scRNA-seq (for Kallisto and cellRanger) and in the gene.activities matrix of scATAC-seq
# Define "simultaneous activity" thresholds
threshold_RNA_all <- c(0,2,5,10,10)
threshold_ATAC_all <- c(0,2,5,5,10)
lncRNA_positions <- which(rownames(gene.activities)%in% lncrna_names)
protein_coding_positions <- which(rownames(gene.activities)%in% protein_coding_names)

for (j in 2:length(threshold_RNA_all))
{
    threshold_ATAC <- threshold_ATAC_all[j]
    threshold_RNA <- threshold_RNA_all[j]

    kallisto_ATAC_expressed_genes <- c()
    cellRanger_ATAC_expressed_genes <- c()
    kallisto_ATAC_expressed_lncRNAs <- c()
    cellRanger_ATAC_expressed_lncRNAs <- c()
    kallisto_ATAC_expressed_PCs_genes <- c()
    cellRanger_ATAC_expressed_PCs_genes <- c()

    for (i in 1:ncol(gene.activities))
    {
        print(i)
        
        kallisto_ATAC_expressed_genes <- c(kallisto_ATAC_expressed_genes,length(intersect(which(gene.activities[,i]>threshold_ATAC),which(as.data.frame(counts(kallisto_GEX[,i]))>threshold_RNA))))
        cellRanger_ATAC_expressed_genes <- c(cellRanger_ATAC_expressed_genes,length(intersect(which(gene.activities[,i]>threshold_ATAC),which(as.data.frame(counts(cellRanger_GEX[,i]))>threshold_RNA))))

        kallisto_ATAC_expressed_lncRNAs <- c(kallisto_ATAC_expressed_lncRNAs,length(intersect(which(gene.activities[lncRNA_positions,i]>threshold_ATAC),which(as.data.frame(counts(kallisto_GEX[lncRNA_positions,i]))>threshold_RNA))))
        cellRanger_ATAC_expressed_lncRNAs <- c(cellRanger_ATAC_expressed_lncRNAs,length(intersect(which(gene.activities[lncRNA_positions,i]>threshold_ATAC),which(as.data.frame(counts(cellRanger_GEX[lncRNA_positions,i]))>threshold_RNA))))

        kallisto_ATAC_expressed_PCs_genes <- c(kallisto_ATAC_expressed_PCs_genes,length(intersect(which(gene.activities[protein_coding_positions,i]>threshold_ATAC),which(as.data.frame(counts(kallisto_GEX[protein_coding_positions,i]))>threshold_RNA))))
        cellRanger_ATAC_expressed_PCs_genes <- c(cellRanger_ATAC_expressed_PCs_genes,length(intersect(which(gene.activities[protein_coding_positions,i]>threshold_ATAC),which(as.data.frame(counts(cellRanger_GEX[protein_coding_positions,i]))>threshold_RNA))))
    
    }

    names(kallisto_ATAC_expressed_genes) <- colnames(gene.activities)
    names(cellRanger_ATAC_expressed_genes) <- colnames(gene.activities)
    names(kallisto_ATAC_expressed_lncRNAs) <- colnames(gene.activities)
    names(cellRanger_ATAC_expressed_lncRNAs) <- colnames(gene.activities)
    names(kallisto_ATAC_expressed_PCs_genes) <- colnames(gene.activities)
    names(cellRanger_ATAC_expressed_PCs_genes) <- colnames(gene.activities)

    common_expressed_cells <- list(kallisto_nuclei_ATAC = kallisto_ATAC_expressed_genes, cellRanger_ATAC = cellRanger_ATAC_expressed_genes, kallisto_nuclei_ATAC_lncRNAs = kallisto_ATAC_expressed_lncRNAs, cellRanger_nuclei_ATAC_lncRNAs = cellRanger_ATAC_expressed_lncRNAs, kallisto_nuclei_ATAC_PCs_genes = kallisto_ATAC_expressed_PCs_genes, cellRanger_nuclei_ATAC_PCs_genes = cellRanger_ATAC_expressed_PCs_genes)
    filename <- paste("common_expressed_genes_threshold_",threshold_ATAC,"ATAC_",threshold_RNA,"RNA.rds",sep="")

    saveRDS(common_expressed_cells, file = paste("validations_paper/",filename,sep=""))
}

# Boxplots

boxplot_data_threshold_0ATAC_0RNA <- boxplot_nuclei_scMULTIOME("validations_paper/common_expressed_genes_threshold_0ATAC_0RNA.rds")
boxplot_data_threshold_2ATAC_2RNA <- boxplot_nuclei_scMULTIOME("validations_paper/common_expressed_genes_threshold_2ATAC_2RNA.rds")
boxplot_data_threshold_5ATAC_5RNA <- boxplot_nuclei_scMULTIOME("validations_paper/common_expressed_genes_threshold_5ATAC_5RNA.rds")
boxplot_data_threshold_5ATAC_10RNA <- boxplot_nuclei_scMULTIOME("validations_paper/common_expressed_genes_threshold_5ATAC_10RNA.rds")
boxplot_data_threshold_10ATAC_10RNA <- boxplot_nuclei_scMULTIOME("validations_paper/common_expressed_genes_threshold_10ATAC_10RNA.rds")

all_boxplots <- rbind(boxplot_data_threshold_0ATAC_0RNA,boxplot_data_threshold_2ATAC_2RNA,boxplot_data_threshold_5ATAC_5RNA,boxplot_data_threshold_5ATAC_10RNA,boxplot_data_threshold_10ATAC_10RNA)
all_boxplots$threshold <- factor(all_boxplots$threshold, levels=c("t0ATAC_0RNA","t2ATAC_2RNA","t5ATAC_5RNA","t5ATAC_10RNA","t10ATAC_10RNA"))

pdf("/home/egonie/kike/phd/test_data/paper_figures/figure3/nuclei_boxplots.pdf")
ggplot(all_boxplots[(all_boxplots$threshold=="t0ATAC_0RNA") | (all_boxplots$threshold=="t2ATAC_2RNA"),], aes(y = common_expression, x = threshold, fill = pipeline)) + geom_boxplot(outlier.shape=NA,position=position_dodge(width=0.8)) + theme_classic() + ggtitle("All genes")+ stat_compare_means(method = "t.test",aes(label = ..p.signif..),label.x = 1.5) +ylim(0,500)
ggplot(all_boxplots[(all_boxplots$threshold!="t0ATAC_0RNA") & (all_boxplots$threshold!="t2ATAC_2RNA"),], aes(y = common_expression, x = threshold, fill = pipeline)) + geom_boxplot(outlier.shape=NA,position=position_dodge(width=0.8)) + theme_classic() + ggtitle("All genes")+ stat_compare_means(method = "t.test",aes(label = ..p.signif..),label.x = 1.5) +ylim(0,25)
ggplot(all_boxplots[(all_boxplots$threshold=="t0ATAC_0RNA") | (all_boxplots$threshold=="t2ATAC_2RNA"),], aes(y = common_expression, x = threshold, fill = pipeline)) + geom_boxplot(outlier.shape=NA,position=position_dodge(width=0.8)) + theme_classic() + ggtitle("All genes")+ stat_compare_means(method = "t.test",symnum.args = list(cutpoints = c(0, 0.0005, 0.005, 0.05, 0.1, 1), symbols = c("****", "***", "**", "*", "ns")),aes(label = ..p.signif..),label.x = 1.5) +ylim(0,500)
ggplot(all_boxplots[(all_boxplots$threshold!="t0ATAC_0RNA") & (all_boxplots$threshold!="t2ATAC_2RNA"),], aes(y = common_expression, x = threshold, fill = pipeline)) + geom_boxplot(outlier.shape=NA,position=position_dodge(width=0.8)) + theme_classic() + ggtitle("All genes")+ stat_compare_means(method = "t.test",symnum.args = list(cutpoints = c(0, 0.0005, 0.005, 0.05, 0.1, 1), symbols = c("****", "***", "**", "*", "ns")),aes(label = ..p.signif..),label.x = 1.5) +ylim(0,25)
ggplot(all_boxplots[(all_boxplots$threshold=="t0ATAC_0RNA") | (all_boxplots$threshold=="t2ATAC_2RNA"),], aes(y = common_expression, x = threshold, fill = pipeline)) + geom_boxplot(outlier.shape=NA,position=position_dodge(width=0.8)) + theme_classic() + ggtitle("All genes")+ stat_compare_means(method = "t.test",label.x = 1.5) +ylim(0,500)
ggplot(all_boxplots[(all_boxplots$threshold!="t0ATAC_0RNA") & (all_boxplots$threshold!="t2ATAC_2RNA"),], aes(y = common_expression, x = threshold, fill = pipeline)) + geom_boxplot(outlier.shape=NA,position=position_dodge(width=0.8)) + theme_classic() + ggtitle("All genes")+ stat_compare_means(method = "t.test",label.x = 1.5) +ylim(0,25)
dev.off()

####
# Ratio of simultaneous activation
###
thresholds_0ATAC_0RNA <- readRDS("validations_paper/common_expressed_genes_threshold_0ATAC_0RNA.rds")
thresholds_2ATAC_2RNA <- readRDS("validations_paper/common_expressed_genes_threshold_2ATAC_2RNA.rds")
thresholds_5ATAC_5RNA <- readRDS("validations_paper/common_expressed_genes_threshold_5ATAC_5RNA.rds")
thresholds_5ATAC_10RNA <- readRDS("validations_paper/common_expressed_genes_threshold_5ATAC_10RNA.rds")
thresholds_10ATAC_10RNA <- readRDS("validations_paper/common_expressed_genes_threshold_10ATAC_10RNA.rds")

# All genes
kallisto_greater <- c(table(thresholds_0ATAC_0RNA[["kallisto_nuclei_ATAC"]]> thresholds_0ATAC_0RNA[["cellRanger_ATAC"]])[2], table(thresholds_2ATAC_2RNA[["kallisto_nuclei_ATAC"]]> thresholds_2ATAC_2RNA[["cellRanger_ATAC"]])[2], table(thresholds_5ATAC_5RNA[["kallisto_nuclei_ATAC"]]> thresholds_5ATAC_5RNA[["cellRanger_ATAC"]])[2], table(thresholds_5ATAC_10RNA[["kallisto_nuclei_ATAC"]]> thresholds_5ATAC_10RNA[["cellRanger_ATAC"]])[2], table(thresholds_10ATAC_10RNA[["kallisto_nuclei_ATAC"]]> thresholds_10ATAC_10RNA[["cellRanger_ATAC"]])[2])
cellRanger_greater <- c(table(thresholds_0ATAC_0RNA[["kallisto_nuclei_ATAC"]]< thresholds_0ATAC_0RNA[["cellRanger_ATAC"]])[2], table(thresholds_2ATAC_2RNA[["kallisto_nuclei_ATAC"]]< thresholds_2ATAC_2RNA[["cellRanger_ATAC"]])[2], table(thresholds_5ATAC_5RNA[["kallisto_nuclei_ATAC"]]< thresholds_5ATAC_5RNA[["cellRanger_ATAC"]])[2], table(thresholds_5ATAC_10RNA[["kallisto_nuclei_ATAC"]]< thresholds_5ATAC_10RNA[["cellRanger_ATAC"]])[2], table(thresholds_10ATAC_10RNA[["kallisto_nuclei_ATAC"]]< thresholds_10ATAC_10RNA[["cellRanger_ATAC"]])[2])
all_genes <- kallisto_greater/cellRanger_greater

# Protein-coding genes
kallisto_greater_PCs <- c(table(thresholds_0ATAC_0RNA[["kallisto_nuclei_ATAC_PCs_genes"]]> thresholds_0ATAC_0RNA[["cellRanger_nuclei_ATAC_PCs_genes"]])[2], table(thresholds_2ATAC_2RNA[["kallisto_nuclei_ATAC_PCs_genes"]]> thresholds_2ATAC_2RNA[["cellRanger_nuclei_ATAC_PCs_genes"]])[2], table(thresholds_5ATAC_5RNA[["kallisto_nuclei_ATAC_PCs_genes"]]> thresholds_5ATAC_5RNA[["cellRanger_nuclei_ATAC_PCs_genes"]])[2], table(thresholds_5ATAC_10RNA[["kallisto_nuclei_ATAC_PCs_genes"]]> thresholds_5ATAC_10RNA[["cellRanger_nuclei_ATAC_PCs_genes"]])[2], table(thresholds_10ATAC_10RNA[["kallisto_nuclei_ATAC_PCs_genes"]]> thresholds_10ATAC_10RNA[["cellRanger_nuclei_ATAC_PCs_genes"]])[2])
cellRanger_greater_PCs <- c(table(thresholds_0ATAC_0RNA[["kallisto_nuclei_ATAC_PCs_genes"]]< thresholds_0ATAC_0RNA[["cellRanger_nuclei_ATAC_PCs_genes"]])[2], table(thresholds_2ATAC_2RNA[["kallisto_nuclei_ATAC_PCs_genes"]]< thresholds_2ATAC_2RNA[["cellRanger_nuclei_ATAC_PCs_genes"]])[2], table(thresholds_5ATAC_5RNA[["kallisto_nuclei_ATAC_PCs_genes"]]< thresholds_5ATAC_5RNA[["cellRanger_nuclei_ATAC_PCs_genes"]])[2], table(thresholds_5ATAC_10RNA[["kallisto_nuclei_ATAC_PCs_genes"]]< thresholds_5ATAC_10RNA[["cellRanger_nuclei_ATAC_PCs_genes"]])[2], table(thresholds_10ATAC_10RNA[["kallisto_nuclei_ATAC_PCs_genes"]]< thresholds_10ATAC_10RNA[["cellRanger_nuclei_ATAC_PCs_genes"]])[2])
PCs <- kallisto_greater_PCs/cellRanger_greater_PCs

# LncRNAs
kallisto_greater_lncRNAs <- c(table(thresholds_0ATAC_0RNA[["kallisto_nuclei_ATAC_lncRNAs"]]> thresholds_0ATAC_0RNA[["cellRanger_nuclei_ATAC_lncRNAs"]])[2], table(thresholds_2ATAC_2RNA[["kallisto_nuclei_ATAC_lncRNAs"]]> thresholds_2ATAC_2RNA[["cellRanger_nuclei_ATAC_lncRNAs"]])[2], table(thresholds_5ATAC_5RNA[["kallisto_nuclei_ATAC_lncRNAs"]]> thresholds_5ATAC_5RNA[["cellRanger_nuclei_ATAC_lncRNAs"]])[2], table(thresholds_5ATAC_10RNA[["kallisto_nuclei_ATAC_lncRNAs"]]> thresholds_5ATAC_10RNA[["cellRanger_nuclei_ATAC_lncRNAs"]])[2], table(thresholds_10ATAC_10RNA[["kallisto_nuclei_ATAC_lncRNAs"]]> thresholds_10ATAC_10RNA[["cellRanger_nuclei_ATAC_lncRNAs"]])[2])
cellRanger_greater_lncRNAs <- c(table(thresholds_0ATAC_0RNA[["kallisto_nuclei_ATAC_lncRNAs"]]< thresholds_0ATAC_0RNA[["cellRanger_nuclei_ATAC_lncRNAs"]])[2], table(thresholds_2ATAC_2RNA[["kallisto_nuclei_ATAC_lncRNAs"]]< thresholds_2ATAC_2RNA[["cellRanger_nuclei_ATAC_lncRNAs"]])[2], table(thresholds_5ATAC_5RNA[["kallisto_nuclei_ATAC_lncRNAs"]]< thresholds_5ATAC_5RNA[["cellRanger_nuclei_ATAC_lncRNAs"]])[2], table(thresholds_5ATAC_10RNA[["kallisto_nuclei_ATAC_lncRNAs"]]< thresholds_5ATAC_10RNA[["cellRanger_nuclei_ATAC_lncRNAs"]])[2], table(thresholds_10ATAC_10RNA[["kallisto_nuclei_ATAC_lncRNAs"]]< thresholds_10ATAC_10RNA[["cellRanger_nuclei_ATAC_lncRNAs"]])[2])
lncRNAs <- kallisto_greater_lncRNAs/cellRanger_greater_lncRNAs

tests <- rep(c("t > 0", "t > 2", "t > 5","t > 10/5", "t > 10"),3)
proportions <- c(log2(all_genes),log2(PCs),log2(lncRNAs))
gene_sets <- c(rep("All",5),rep("Protein-coding",5),rep("LncRNAs",5))
df <- data.frame(cbind(tests,gene_sets,proportions))
df$proportions <- as.numeric(df$proportions)
df$tests <- as.factor(df$test)
df$gene_sets <- as.factor(df$gene_sets)
# custom order
df$tests <- factor(df$tests, levels = c("t > 0", "t > 2", "t > 5","t > 10/5", "t > 10"))
df$gene_sets <- factor(df$gene_sets, levels = c("All","Protein-coding", "LncRNAs"))
colors <- c("#7f64b9","#9f9244","#c36785")

pdf("/home/egonie/kike/phd/test_data/paper_figures/figure3/nuclei_ratios_genes_simultaneous_activation.pdf")
ggplot(data=df, aes(x=tests, y=proportions, fill = gene_sets)) + geom_bar(stat="identity",position = position_dodge(), width = 0.5, color="black")  +
theme_classic() +scale_fill_manual(values=colors) + ylab("log2(Genes Kallisto/Genes cellRanger)") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 20), axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 17))
dev.off()

####
# Odds ratio simultaneous activation
####
# Odds ratio genes
OR_t_0_genes <- odds_ratio_genes(input_RDS = "validations_paper/common_expressed_genes_threshold_0ATAC_0RNA.rds", gene.activities)
OR_t_2_genes <- odds_ratio_genes(input_RDS = "validations_paper/common_expressed_genes_threshold_2ATAC_2RNA.rds", gene.activities)
OR_t_5_genes <- odds_ratio_genes(input_RDS = "validations_paper/common_expressed_genes_threshold_5ATAC_5RNA.rds", gene.activities)
OR_t_5ATAC_10RNA_genes <- odds_ratio_genes(input_RDS = "validations_paper/common_expressed_genes_threshold_5ATAC_10RNA.rds", gene.activities)
OR_t_10_genes <- odds_ratio_genes(input_RDS = "validations_paper/common_expressed_genes_threshold_10ATAC_10RNA.rds", gene.activities)


df_threshold <- data.frame(pipeline = c(rep("t > 0",length(OR_t_0_genes[[1]])), rep("t > 2",length(OR_t_2_genes[[1]])), rep("t > 5",length(OR_t_5_genes[[1]])), rep("t > 10/5", length(OR_t_5ATAC_10RNA_genes[[1]])), rep("t > 10",length(OR_t_10_genes[[1]]))),OR = c(OR_t_0_genes[[1]], OR_t_2_genes[[1]], OR_t_5_genes[[1]], OR_t_5ATAC_10RNA_genes[[1]], OR_t_10_genes[[1]]))
df_threshold$pipeline <- factor(df_threshold$pipeline, levels = c("t > 0", "t > 2", "t > 5", "t > 10/5", "t > 10"))

df_threshold_lncRNAs <- data.frame(pipeline = c(rep("t > 0",length(OR_t_0_genes[[2]])), rep("t > 2",length(OR_t_2_genes[[2]])), rep("t > 5",length(OR_t_5_genes[[2]])), rep("t > 10/5", length(OR_t_5ATAC_10RNA_genes[[2]])), rep("t > 10",length(OR_t_10_genes[[2]]))),OR = c(OR_t_0_genes[[2]], OR_t_2_genes[[2]], OR_t_5_genes[[2]], OR_t_5ATAC_10RNA_genes[[2]], OR_t_10_genes[[2]]))
df_threshold_lncRNAs$pipeline <- factor(df_threshold_lncRNAs$pipeline, levels = c("t > 0", "t > 2", "t > 5", "t > 10/5", "t > 10"))

df_threshold_PCs <- data.frame(pipeline = c(rep("t > 0",length(OR_t_0_genes[[3]])), rep("t > 2",length(OR_t_2_genes[[3]])), rep("t > 5",length(OR_t_5_genes[[3]])), rep("t > 10/5", length(OR_t_5ATAC_10RNA_genes[[3]])), rep("t > 10",length(OR_t_10_genes[[3]]))),OR = c(OR_t_0_genes[[3]], OR_t_2_genes[[3]], OR_t_5_genes[[3]], OR_t_5ATAC_10RNA_genes[[3]], OR_t_10_genes[[3]]))
df_threshold_PCs$pipeline <- factor(df_threshold_PCs$pipeline, levels = c("t > 0", "t > 2", "t > 5", "t > 10/5", "t > 10"))

#plot all OR together 
df_threshold$genes <- "All"
df_threshold_lncRNAs$genes <- "LncRNAs"
df_threshold_PCs$genes <- "Protein-coding"

all_df <- rbind(df_threshold,df_threshold_lncRNAs,df_threshold_PCs)
all_df$pipeline <- factor(all_df$pipeline, levels = c("t > 0", "t > 2", "t > 5","t > 10/5", "t > 10"))
all_df$genes <- factor(all_df$genes, levels = c("All","Protein-coding", "LncRNAs"))
colors_5 <- c("#ba543d","#56ae6c","#7f63b8","#ac9c3d","#b84c7d")

pdf("/home/egonie/kike/phd/test_data/paper_figures/figure3/nuclei_simultaneous_activation_ODDS_ratio.pdf")
ggplot(df_threshold, aes(x = pipeline, y = OR, fill = pipeline)) + geom_violin(position=position_dodge(width=0.8))  + theme_classic() + geom_hline(aes(yintercept=0),linetype="dashed")  + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 14),axis.text.y = element_text(size = 14))+  scale_fill_manual(values = colors_5) +  ggtitle("Nuclei ODDS RATIO: All genes")
ggplot(df_threshold_PCs, aes(x = pipeline, y = OR, fill = pipeline)) + geom_violin(position=position_dodge(width=0.8))  + theme_classic() + geom_hline(aes(yintercept=0),linetype="dashed")  + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 14),axis.text.y = element_text(size = 14))+ ggtitle("Nuclei ODDS RATIO: Protein-coding genes")+ scale_fill_manual(values = colors_5)
ggplot(df_threshold_lncRNAs, aes(x = pipeline, y = OR, fill = pipeline)) + geom_violin(position=position_dodge(width=0.8))  + theme_classic() + geom_hline(aes(yintercept=0),linetype="dashed")  + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 14),axis.text.y = element_text(size = 14))+ ggtitle("Nuclei ODDS RATIO: lncRNAs")+ scale_fill_manual(values = colors_5)
dev.off()

# Coverage plots in order to see genes whose scATAC-seq measurements are only corresponded with scRNA-seq when preprocessed with Kallisto
# Also showing potential conflicts with the intronic annotation
gtf$tx_id=gtf$transcript_id
r= CreateChromatinAssay(counts = counts$Peaks[,colnames(cellRanger_GEX)],sep = c(":", "-"),fragments = fragpath,annotation = gtf)
pbmc_CR = as.Seurat(cellRanger_GEX)
pbmc_CR[['ATAC']] <- r
DefaultAssay(pbmc_CR) <- "ATAC"

pbmc_kallisto = as.Seurat(kallisto_GEX)
pbmc_kallisto[['ATAC']] <- r
DefaultAssay(pbmc_kallisto) <- "ATAC"

pdf("/home/egonie/kike/phd/test_data/paper_figures/figure3/candidates_coverage_tracks_GEX.pdf",width=18)
candidates <- c("CYP2F1","AC243960.3")
for(i in candidates)
{
    print(i)
    exp_t <- ExpressionPlot(object = pbmc_CR,features = i,assay = "RNA",group.by="cell_type") + xlim(0,2.5)
    cov_t <- CoveragePlot(object = pbmc_CR,region = i,annotation = TRUE,peaks = FALSE,group.by="cell_type",extend.upstream = 2000,extend.downstream = 1000, window = 200)
    p <- cov_t + exp_t + plot_layout(widths = c(10,2))
    print(p & scale_fill_manual(values = c("#c8723e","#9071c7","#84a142","#cb547e","#bd5db0","#4daf95")) & ggtitle("CellRanger GEX & ATAC coverage"))

    exp_t <- ExpressionPlot(object = pbmc_kallisto,features = i,assay = "RNA",group.by="cell_type") + xlim(0,2.5)
    cov_t <- CoveragePlot(object = pbmc_kallisto,region = i,annotation = TRUE,peaks = FALSE,group.by="cell_type",extend.upstream = 2000,extend.downstream = 1000, window = 200)
    p <- cov_t + exp_t + plot_layout(widths = c(10,2))
    print(p & scale_fill_manual(values = c("#c8723e","#9071c7","#84a142","#cb547e","#bd5db0","#4daf95")) & ggtitle("Kallisto GEX & ATAC coverage"))
}
dev.off()

intronic_annotation_conflict_candidates <- c("IGKC","JUNB")
pdf("/home/egonie/kike/phd/test_data/paper_figures/figure3/candidates_coverage_intronic_annotation_conflicts.pdf",width=18)
for(i in intronic_annotation_conflict_candidates)
{
    p <- CoveragePlot(object = pbmc_CR,region = i,annotation = TRUE,peaks = FALSE, features = i,group.by="cell_type",extend.upstream = 2000,extend.downstream = 1000, window = 200, heights = c(10,2),widths=c(10,2))& ggtitle("CellRanger GEX & ATAC coverage")
    print(p & scale_fill_manual(values = c("#c8723e","#9071c7","#84a142","#cb547e","#bd5db0","#4daf95")))
    p <- CoveragePlot(object = pbmc_kallisto,region = i ,annotation = TRUE,peaks = FALSE, features = i,group.by="cell_type",extend.upstream = 2000,extend.downstream = 1000, window = 200,heights = c(10,2),widths=c(10,2)) & ggtitle("Kallisto GEX & ATAC coverage")
    print(p & scale_fill_manual(values = c("#c8723e","#9071c7","#84a142","#cb547e","#bd5db0","#4daf95")))
}
dev.off()
