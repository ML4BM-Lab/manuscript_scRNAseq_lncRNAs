########################################################### Generic #################################################################################################
#####################################################################################################################################################################
libraries <- c("Seurat","patchwork", "tximport", "ggplot2", "scran" , "scater", "SingleCellExperiment", "BUSpaRse", "DropletUtils", "RColorBrewer", "cluster", "UpSetR","scRNAseq", "ggpubr","scDblFinder","reshape","Signac")

lapply(libraries, require, character.only = TRUE)
# Import gtf annotation & setWD, source functions
gencode_path <- "~/dato-activo/reference.genomes_kike/GRCh38/gencode/gencode.v37.annotation.gtf"
gtf <- rtracklayer::import(gencode_path)
gtf$gene_biotype <- gtf$gene_type

mitochondrial_ens_ids <- unique(hg38_ensembl_gtf$gene_id[grep("^MT-",hg38_ensembl_gtf$gene_name)])
lncrna_ens_ids <- unique(c(hg38_ensembl_gtf$gene_id[grep("lncRNA",hg38_ensembl_gtf$gene_type)]))
protein_coding_ens_ids <- unique(c(hg38_ensembl_gtf$gene_id[hg38_ensembl_gtf$gene_type=="protein_coding"]))

lncrna_names <- unique(hg38_ensembl_gtf$gene_name[hg38_ensembl_gtf$gene_id %in% lncrna_ens_ids])
protein_coding_names <-  unique(hg38_ensembl_gtf$gene_name[hg38_ensembl_gtf$gene_id %in% protein_coding_ens_ids])

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


