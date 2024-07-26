########################################################### Generic #################################################################################################
#####################################################################################################################################################################
libraries <- c("Seurat","patchwork", "tximport", "ggplot2", "scran" , "scater", "SingleCellExperiment", "BUSpaRse", "DropletUtils", "RColorBrewer", "cluster", "UpSetR","scRNAseq", "ggpubr","scDblFinder","reshape","dplyr","ggupset")

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

# Test if detection is affected when matching expression levels of protein-coding genes and lncRNAs
all_lncRNAs_expression <- c((rowSums(logcounts(cellRanger_filt_sce))[lncrna_ens_ids[lncrna_ens_ids %in% rownames(cellRanger_filt_sce)==T]]), (rowSums(logcounts(STARsolo_filt_sce))[lncrna_ens_ids[lncrna_ens_ids %in% rownames(STARsolo_filt_sce)==T]]), (rowSums(logcounts(kallisto_filt_sce))[lncrna_ens_ids[lncrna_ens_ids %in% rownames(kallisto_filt_sce)==T]]), (rowSums(logcounts(alevin_filt_sce))[lncrna_ens_ids[lncrna_ens_ids %in% rownames(alevin_filt_sce)==T]]))

for (i in c(0.5,1,2))
{
  max_gene_expression <- mean(all_lncRNAs_expression)+i*sd(all_lncRNAs_expression)

  cellRanger_filt_genes <- filtered_gene_expression(sce=cellRanger_filt_sce, max_gene_expression=max_gene_expression, lncrna_ens_ids,protein_coding_ens_ids )
  STARsolo_filt_genes <- filtered_gene_expression(sce=STARsolo_filt_sce, max_gene_expression=max_gene_expression, lncrna_ens_ids,protein_coding_ens_ids )
  kallisto_filt_genes <- filtered_gene_expression(sce=kallisto_filt_sce, max_gene_expression=max_gene_expression, lncrna_ens_ids,protein_coding_ens_ids )
  alevin_filt_genes <- filtered_gene_expression(sce=alevin_filt_sce, max_gene_expression=max_gene_expression, lncrna_ens_ids,protein_coding_ens_ids )

  detected_lncRNAs_CellRanger <- colSums(logcounts(cellRanger_filt_genes[intersect(rownames(cellRanger_filt_genes),lncrna_ens_ids)],)!=0)
  detected_lncRNAs_STARsolo <- colSums(logcounts(STARsolo_filt_genes[intersect(rownames(STARsolo_filt_genes),lncrna_ens_ids)],)!=0)
  detected_lncRNAs_kallisto <- colSums(logcounts(kallisto_filt_genes[intersect(rownames(kallisto_filt_genes),lncrna_ens_ids)],)!=0)
  detected_lncRNAs_alevin <- colSums(logcounts(alevin_filt_genes[intersect(rownames(alevin_filt_genes),lncrna_ens_ids)],)!=0)

  detected_PCs_CellRanger <- colSums(logcounts(cellRanger_filt_genes[intersect(rownames(cellRanger_filt_genes),protein_coding_ens_ids)],)!=0)
  detected_PCs_STARsolo <- colSums(logcounts(STARsolo_filt_genes[intersect(rownames(STARsolo_filt_genes),protein_coding_ens_ids)],)!=0)
  detected_PCs_kallisto <- colSums(logcounts(kallisto_filt_genes[intersect(rownames(kallisto_filt_genes),protein_coding_ens_ids)],)!=0)
  detected_PCs_alevin <- colSums(logcounts(alevin_filt_genes[intersect(rownames(alevin_filt_genes),protein_coding_ens_ids)],)!=0)

  all_sce_lncRNAs_detected <- combine_4sce(detected_lncRNAs_CellRanger, detected_lncRNAs_STARsolo, detected_lncRNAs_kallisto, detected_lncRNAs_alevin, "detected_filtered_lncRNAs" )
  all_sce_PCs_detected <- combine_4sce(detected_PCs_CellRanger, detected_PCs_STARsolo, detected_PCs_kallisto, detected_PCs_alevin, "detected_filtered_PCGs" )

  pdf(paste("QC_lncRNAs_PCGs_similar_counts_mean_plus_",i,"sd.pdf",sep=""))
  ggplot(all_sce_lncRNAs_detected, aes(x=ident, y=detected_filtered_lncRNAs, fill = ident)) + geom_violin(alpha = 1, scale = "width", width = 0.8) + scale_fill_manual(values = color_palette) + theme_classic() + scale_colour_manual(values = color_palette) + ggtitle("Filtered lncRNAs detected/cell") + theme(plot.title = element_text(hjust = 0.5, size = 28),axis.text=element_text(size=17), axis.title=element_text(size=17, face = "bold")) + theme(legend.position = "none") + labs(y= "", x="")

  ggplot(all_sce_PCs_detected, aes(x=ident, y=detected_filtered_PCGs, fill = ident)) + geom_violin(alpha = 1, scale = "width", width = 0.8) + scale_fill_manual(values = color_palette) + theme_classic() + scale_colour_manual(values = color_palette) + ggtitle("Filtered PCGs detected/cell") + theme(plot.title = element_text(hjust = 0.5, size = 28),axis.text=element_text(size=17), axis.title=element_text(size=17, face = "bold")) + theme(legend.position = "none") + labs(y= "", x="")
  dev.off()
}



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
input_list_lncRNAs <- list(CellRanger = intersect(rownames(cellRanger_top_genes),lncrna_ens_ids), Kallisto = intersect(rownames(kallisto_top_genes),lncrna_ens_ids))
ratio_lncRNAs_v38 <- length(setdiff(input_list_lncRNAs$Kallisto, input_list_lncRNAs$CellRanger))/length(intersect(input_list_lncRNAs$Kallisto, input_list_lncRNAs$CellRanger))

input_list2 <- list(CellRanger = intersect(rownames(cellRanger_top_genes),protein_coding_ens_ids), STARsolo = intersect(rownames(STARsolo_top_genes),protein_coding_ens_ids), Kallisto = intersect(rownames(kallisto_top_genes),protein_coding_ens_ids), Salmon = intersect(rownames(alevin_top_genes),protein_coding_ens_ids))
pdf("upset_plot_protein_coding_names_250counts_in_25_genes.pdf", width = 8, height = 8)
upset(fromList(input_list2), order.by = "degree", mainbar.y.label = "", keep.order = TRUE, text.scale = c(2.5,1.7,1.7,1.6,2.0,1.7), sets.x.label = "Protein coding genes") 
#grid.text("Protein coding genes intersection",x = 0.65, y=0.97, gp=gpar(fontsize=12, fontface="bold"))
dev.off()
input_list_PCs <- list(CellRanger = intersect(rownames(cellRanger_top_genes),protein_coding_ens_ids), Kallisto = intersect(rownames(kallisto_top_genes),protein_coding_ens_ids))
ratio_PCs_v38 <- length(setdiff(input_list_PCs$Kallisto, input_list_PCs$CellRanger))/length(intersect(input_list_PCs$Kallisto, input_list_PCs$CellRanger))



##########################################################################################################################################################################
# Detection of protein-coding genes & lncRNAs affected by more/less precise annotations (compare commonly detected vs exclusively detected)
# 
hg38_ensembl_gtf_ <- as.data.frame(rtracklayer::import("/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode_release_45/gencode..annotation.gtf"))
hg38_ensembl_gtf_$gene_id <- gsub("_","-",hg38_ensembl_gtf_$gene_id)
mitochondrial_ens_ids_ <- unique(hg38_ensembl_gtf_$gene_id[grep("^MT-",hg38_ensembl_gtf_$gene_name)])
lncrna_ens_ids_hg38_ <- unique(c(hg38_ensembl_gtf_$gene_id[grep("lncRNA",hg38_ensembl_gtf_$gene_type)]))
protein_coding_ens_ids_hg38_ <- unique(c(hg38_ensembl_gtf_$gene_id[hg38_ensembl_gtf_$gene_type=="protein_coding"]))

cellRanger_ <- CreateSeuratObject(Read10X(data.dir = "/home/egonie/kike/phd/test_data/10X/Parent_NGSC3_DI_PBMC_fastqs/hg38_gencode_release_/01.CellRanger/Parent_NGSC3_DI_PBMC/outs/raw_feature_bc_matrix", gene.column = 1), project = "cellRanger")
kallisto_ <- CreateSeuratObject(BUSpaRse::read_count_output("/home/egonie/kike/phd/test_data/10X/Parent_NGSC3_DI_PBMC_fastqs/hg38_gencode_release_/01.Kallisto/output_bus_transcriptome/bustools_results", name = "cells_genes_NO_multimapping"), project = "kallisto")

# Convert to SingleCellExperiment
cellRanger_ <- as.SingleCellExperiment(cellRanger_)
colnames(cellRanger_) <- gsub("\\-.*","",colnames(cellRanger_))
kallisto_ <- as.SingleCellExperiment(kallisto_)

# Quality control
cellRanger_sce_ <- qc_metrics(cellRanger_, mitochondrial_ens_ids)
kallisto_sce_ <- qc_metrics(kallisto_, mitochondrial_ens_ids)

# EmptyDrops filtering
cellRanger_filt_sce_ed_ <- emptydrops_filt(cellRanger_sce_, lower_ED = 100, EmptyDrops_FDR_thres = 0.0001 )
kallisto_filt_sce_ed_ <- emptydrops_filt(kallisto_sce_, lower_ED = 100, EmptyDrops_FDR_thres = 0.0001 )

# Doublet analysis
cellRanger_filt_sce_ed_nodoubs_ <- doublet_analysis(cellRanger_filt_sce_ed_)
cellRanger_filt_sce_ed_nodoubs_ <- cellRanger_filt_sce_ed_nodoubs_[,cellRanger_filt_sce_ed_nodoubs_$isDoublet == F]
kallisto_filt_sce_ed_nodoubs_ <- doublet_analysis(kallisto_filt_sce_ed_)
kallisto_filt_sce_ed_nodoubs_ <- kallisto_filt_sce_ed_nodoubs_[,kallisto_filt_sce_ed_nodoubs_$isDoublet == F]

pbmc_datasets_ <- list("cellRanger" = cellRanger_filt_sce_ed_nodoubs_, "kallisto" = kallisto_filt_sce_ed_nodoubs_)
saveRDS(pbmc_datasets_, file = "/home/egonie/kike/phd/test_data/paper_figures/figure1/pbmc_datasets_hg38_.rds")
pbmc_datasets_ <- readRDS("/home/egonie/kike/phd/test_data/paper_figures/figure1/pbmc_datasets_hg38_.rds")
cellRanger_filt_sce_ed_nodoubs_ <- pbmc_datasets_[["cellRanger"]]
kallisto_filt_sce_ed_nodoubs_ <- pbmc_datasets_[["kallisto"]]

cellRanger_filt_sce_ed_ <- Filtering(cellRanger_filt_sce_ed_nodoubs_,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
kallisto_filt_sce_ed_ <- Filtering(kallisto_filt_sce_ed_nodoubs_,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)

kallisto_top_genes_ <- top_genes(kallisto_filt_sce_ed_,threshold_minumun_gene_counts,threshold_cells_detected )
cellRanger_top_genes_ <- top_genes(cellRanger_filt_sce_ed_,threshold_minumun_gene_counts,threshold_cells_detected )

input_list__lncRNAs <- list(CellRanger = intersect(rownames(cellRanger_top_genes_),lncrna_ens_ids_hg38_), Kallisto = intersect(rownames(kallisto_top_genes_),lncrna_ens_ids_hg38_))
ratio_lncRNAs_ <- length(setdiff(input_list__lncRNAs$Kallisto, input_list__lncRNAs$CellRanger))/length(intersect(input_list__lncRNAs$Kallisto, input_list__lncRNAs$CellRanger))

input_list__PCs <- list(CellRanger = intersect(rownames(cellRanger_top_genes_),protein_coding_ens_ids_hg38_), Kallisto = intersect(rownames(kallisto_top_genes_),protein_coding_ens_ids_hg38_))
ratio_PCs_ <- length(setdiff(input_list__PCs$Kallisto, input_list__PCs$CellRanger))/length(intersect(input_list__PCs$Kallisto, input_list__PCs$CellRanger))



# hg19 v19
hg19_ensembl_gtf <- as.data.frame(rtracklayer::import("/home/egonie/dato-activo/reference.genomes_kike/GRCh37/gencode/gencode.v19.annotation.gtf"))
hg19_ensembl_gtf$gene_id <- gsub("_","-",hg19_ensembl_gtf$gene_id)
mitochondrial_ens_ids_v19 <- unique(hg19_ensembl_gtf$gene_id[grep("^MT-",hg19_ensembl_gtf$gene_name)])
lncrna_ens_ids_v19 <- unique(c(hg19_ensembl_gtf$gene_id[(hg19_ensembl_gtf$gene_type=="3prime_overlapping_ncRNA") | (hg19_ensembl_gtf$gene_type=="antisense") | (hg19_ensembl_gtf$gene_type=="bidirectional_promoter_lncRNA") | (hg19_ensembl_gtf$gene_type=="lincRNA") | (hg19_ensembl_gtf$gene_type=="macro_lncRNA") | (hg19_ensembl_gtf$gene_type=="non_coding") | (hg19_ensembl_gtf$gene_type=="processed_transcript") | (hg19_ensembl_gtf$gene_type=="sense_intronic") | (hg19_ensembl_gtf$gene_type=="sense_overlapping")]))
protein_coding_ens_ids_v19 <- unique(c(hg19_ensembl_gtf$gene_id[hg19_ensembl_gtf$gene_type=="protein_coding"]))

cellRanger_v19 <- CreateSeuratObject(Read10X(data.dir = "/home/egonie/kike/phd/test_data/10X/Parent_NGSC3_DI_PBMC_fastqs/hg19_gencode/01.CellRanger/Parent_NGSC3_DI_PBMC/outs/raw_feature_bc_matrix", gene.column = 1), project = "cellRanger")
kallisto_v19 <- CreateSeuratObject(BUSpaRse::read_count_output("/home/egonie/kike/phd/test_data/10X/Parent_NGSC3_DI_PBMC_fastqs/hg19_gencode/01.Kallisto/output_bus_transcriptome/bustools_results_no_multimappers", name = "cells_genes_NO_multimapping"), project = "kallisto")

# Convert to SingleCellExperiment
cellRanger_v19 <- as.SingleCellExperiment(cellRanger_v19)
colnames(cellRanger_v19) <- gsub("\\-.*","",colnames(cellRanger_v19))
kallisto_v19 <- as.SingleCellExperiment(kallisto_v19)

# Quality control
cellRanger_sce_v19 <- qc_metrics(cellRanger_v19, mitochondrial_ens_ids_v19)
kallisto_sce_v19 <- qc_metrics(kallisto_v19, mitochondrial_ens_ids_v19)

# EmptyDrops filtering
cellRanger_filt_sce_ed_v19 <- emptydrops_filt(cellRanger_sce_v19, lower_ED = 100, EmptyDrops_FDR_thres = 0.0001 )
kallisto_filt_sce_ed_v19 <- emptydrops_filt(kallisto_sce_v19, lower_ED = 100, EmptyDrops_FDR_thres = 0.0001 )

# Doublet analysis
cellRanger_filt_sce_ed_nodoubs_v19 <- doublet_analysis(cellRanger_filt_sce_ed_v19)
cellRanger_filt_sce_ed_nodoubs_v19 <- cellRanger_filt_sce_ed_nodoubs_v19[,cellRanger_filt_sce_ed_nodoubs_v19$isDoublet == F]
kallisto_filt_sce_ed_nodoubs_v19 <- doublet_analysis(kallisto_filt_sce_ed_v19)
kallisto_filt_sce_ed_nodoubs_v19 <- kallisto_filt_sce_ed_nodoubs_v19[,kallisto_filt_sce_ed_nodoubs_v19$isDoublet == F]

pbmc_datasets_v19 <- list("cellRanger" = cellRanger_filt_sce_ed_nodoubs_v19, "kallisto" = kallisto_filt_sce_ed_nodoubs_v19)
saveRDS(pbmc_datasets_v19, file = "/home/egonie/kike/phd/test_data/paper_figures/figure1/pbmc_datasets_hg19_v19.rds")
pbmc_datasets_v19 <- readRDS("/home/egonie/kike/phd/test_data/paper_figures/figure1/pbmc_datasets_hg19_v19.rds")
cellRanger_filt_sce_ed_nodoubs_v19 <- pbmc_datasets_v19[["cellRanger"]]
kallisto_filt_sce_ed_nodoubs_v19 <- pbmc_datasets_v19[["kallisto"]]

cellRanger_filt_sce_ed_v19 <- Filtering(cellRanger_filt_sce_ed_nodoubs_v19,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
kallisto_filt_sce_ed_v19 <- Filtering(kallisto_filt_sce_ed_nodoubs_v19,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)

kallisto_top_genes_v19 <- top_genes(kallisto_filt_sce_ed_v19,threshold_minumun_gene_counts,threshold_cells_detected )
cellRanger_top_genes_v19 <- top_genes(cellRanger_filt_sce_ed_v19,threshold_minumun_gene_counts,threshold_cells_detected )

input_list_hg19_lncRNAs <- list(CellRanger = intersect(rownames(cellRanger_top_genes_v19),lncrna_ens_ids_v19), Kallisto = intersect(rownames(kallisto_top_genes_v19),lncrna_ens_ids_v19))
ratio_lncRNAs_v19 <- length(setdiff(input_list_hg19_lncRNAs$Kallisto, input_list_hg19_lncRNAs$CellRanger))/length(intersect(input_list_hg19_lncRNAs$Kallisto, input_list_hg19_lncRNAs$CellRanger))

input_list_hg19_PCs <- list(CellRanger = intersect(rownames(cellRanger_top_genes_v19),protein_coding_ens_ids_v19), Kallisto = intersect(rownames(kallisto_top_genes_v19),protein_coding_ens_ids_v19))
ratio_PCs_v19 <- length(setdiff(input_list_hg19_PCs$Kallisto, input_list_hg19_PCs$CellRanger))/length(intersect(input_list_hg19_PCs$Kallisto, input_list_hg19_PCs$CellRanger))

# NONCODE (v5)
cellRanger_NONCODE <- CreateSeuratObject(Read10X(data.dir = "/home/egonie/kike/phd/test_data/10X/Parent_NGSC3_DI_PBMC_fastqs/NONCODE/01.CellRanger/Parent_NGSC3_DI_PBMC/outs/raw_feature_bc_matrix", gene.column = 1), project = "cellRanger")
kallisto_NONCODE <- CreateSeuratObject(BUSpaRse::read_count_output("/home/egonie/kike/phd/test_data/10X/Parent_NGSC3_DI_PBMC_fastqs/NONCODE/01.Kallisto/output_bus_transcriptome/bustools_results", name = "cells_genes_NO_multimapping"), project = "kallisto")

cellRanger_NONCODE <- as.SingleCellExperiment(cellRanger_NONCODE)
kallisto_NONCODE <- as.SingleCellExperiment(kallisto_NONCODE)

universe_genes <- intersect(rownames(cellRanger_NONCODE), rownames(kallisto_NONCODE))
cellRanger_NONCODE <- cellRanger_NONCODE[universe_genes,]
kallisto_NONCODE <- kallisto_NONCODE[universe_genes,]

qc <- function(pipeline )
{
  sce <- pipeline  
  stats <- perCellQCMetrics(sce)
  sum_log=data.frame(log10(stats$sum))
  detected_log=data.frame(log10(stats$detected))
  colData(sce) <- c(colData(sce),stats,sum_log,detected_log)
  
  return(sce)
}
cellRanger_sce_NONCODE <- qc(cellRanger_NONCODE)
kallisto_sce_NONCODE <- qc(kallisto_NONCODE)

cellRanger_filt_sce_ed_NONCODE <- emptydrops_filt(cellRanger_NONCODE, lower_ED = 100, EmptyDrops_FDR_thres = 0.0001 )
kallisto_filt_sce_ed_NONCODE <- emptydrops_filt(kallisto_NONCODE, lower_ED = 100, EmptyDrops_FDR_thres = 0.0001 )

cellRanger_filt_sce_ed_nodoubs_NONCODE <- doublet_analysis(cellRanger_filt_sce_ed_NONCODE)
cellRanger_filt_sce_ed_nodoubs_NONCODE <- cellRanger_filt_sce_ed_nodoubs_NONCODE[,cellRanger_filt_sce_ed_nodoubs_NONCODE$isDoublet == F]
kallisto_filt_sce_ed_nodoubs_NONCODE <- doublet_analysis(kallisto_filt_sce_ed_NONCODE)
kallisto_filt_sce_ed_nodoubs_NONCODE <- kallisto_filt_sce_ed_nodoubs_NONCODE[,kallisto_filt_sce_ed_nodoubs_NONCODE$isDoublet == F]

f <- function(sce, cells_max_threshold, cells_min_genes_detected_threshold)
{
    total_max_expression <- sce$nCount_RNA > cells_max_threshold
    total_min_genes_detected <- sce$nFeature_RNA < cells_min_genes_detected_threshold

    discard <-  total_max_expression  | total_min_genes_detected
    print(table(discard))
    print(DataFrame(total_max_expression=sum(total_max_expression),total_min_genes_detected=sum(total_min_genes_detected), Total=sum(discard)))
    sce$discard <- discard
    #filter and normalize data
    sce_filt <- sce[,discard==F]

    sce_filt <- logNormCounts(sce_filt)

    return(sce_filt)
}

cellRanger_filt_sce_NONCODE <- f(cellRanger_filt_sce_ed_nodoubs_NONCODE, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
kallisto_filt_sce_NONCODE <- f(kallisto_filt_sce_ed_nodoubs_NONCODE, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)

cellRanger_top_genes_NONCODE <- top_genes(cellRanger_filt_sce_NONCODE,threshold_minumun_gene_counts,threshold_cells_detected )
kallisto_top_genes_NONCODE <- top_genes(kallisto_filt_sce_NONCODE,threshold_minumun_gene_counts,threshold_cells_detected )
input_list_NONCODE <- list(CellRanger = rownames(cellRanger_top_genes_NONCODE), Kallisto =rownames(kallisto_top_genes_NONCODE))
ratio_lncRNAs_NONCODE <- length(setdiff(input_list_NONCODE$Kallisto, input_list_NONCODE$CellRanger))/length(intersect(input_list_NONCODE$Kallisto, input_list_NONCODE$CellRanger))


# Ratio of = Number of Kallisto exclusive lncRNAs / Number of commonly detected lncRNAs
df <- data.frame(assembly = c("GENCODE_hg19", "GENCODE_hg38_v37","GENCODE_hg38_","NONCODE"),Ratio_exclusive_common = c(ratio_lncRNAs_v19/ratio_lncRNAs_v19, ratio_lncRNAs_v38/ratio_lncRNAs_v19, ratio_lncRNAs_/ratio_lncRNAs_v19, ratio_lncRNAs_NONCODE/ratio_lncRNAs_v19))
pdf("ratios_exclusive_common_lncRNAs.pdf", width = 10)
no_theme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggplot(df, aes(x=assembly, y=Ratio_exclusive_common, fill=assembly)) + 
  geom_bar(stat="identity",color="black",show.legend = FALSE) +no_theme + theme(axis.title.y = element_text(size=15, hjust = 0.5, vjust = -1),plot.title = element_text(hjust = 0.5, size = 14, face="bold"),axis.text=element_text(size=12, color="black")) + labs(title="Ratio exclusive to common lncRNAs", x="Pipeline", y = "Ratio (norm to GENCODE_hg19)") + ylim(0, 1) + scale_x_discrete(limits=c("GENCODE_hg19", "GENCODE_hg38_v37","GENCODE_hg38_","NONCODE"))
dev.off()
# Also calculate the ratios for protein-coding genes
df <- data.frame(assembly = c("GENCODE_hg19", "GENCODE_hg38_v37","GENCODE_hg38_"),Ratio_exclusive_common = c(ratio_PCs_v19/ratio_PCs_v19, ratio_PCs_v38/ratio_PCs_v19, ratio_PCs_/ratio_PCs_v19))
pdf("ratios_exclusive_common_PCs.pdf", width = 10)
no_theme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggplot(df, aes(x=assembly, y=Ratio_exclusive_common, fill=assembly)) + 
  geom_bar(stat="identity",color="black",show.legend = FALSE) +no_theme + theme(axis.title.y = element_text(size=15, hjust = 0.5, vjust = -1),plot.title = element_text(hjust = 0.5, size = 14, face="bold"),axis.text=element_text(size=12, color="black")) + labs(title="Ratio exclusive to common protein-coding", x="Pipeline", y = "Ratio (norm to GENCODE_hg19)") + ylim(0, 1) + scale_x_discrete(limits=c("GENCODE_hg19", "GENCODE_hg38_v37","GENCODE_hg38_"))
dev.off()

# Also represent as UpSet plots
#lncRNAs

all_sce_melted <- as.data.frame(cbind(as.character(c(rep("kallisto",length(input_list_hg19_lncRNAs[[2]])),rep("kallisto",length(input_list_lncRNAs[[2]])),rep("kallisto",length(input_list__lncRNAs[[2]])),rep("kallisto",length(input_list_NONCODE[[2]])),rep("CellRanger",length(input_list_hg19_lncRNAs[[1]])),rep("CellRanger",length(input_list_lncRNAs[[1]])),rep("CellRanger",length(input_list__lncRNAs[[1]])),rep("CellRanger",length(input_list_NONCODE[[1]])))),as.character(c(input_list_hg19_lncRNAs[[2]],input_list_lncRNAs[[2]], input_list__lncRNAs[[2]],input_list_NONCODE[[2]],input_list_hg19_lncRNAs[[1]],input_list_lncRNAs[[1]], input_list__lncRNAs[[1]],input_list_NONCODE[[1]])),as.character(c(rep("GENCODE_hg19",length(input_list_hg19_lncRNAs[[2]])),rep("GENCODE_hg38_v37",length(input_list_lncRNAs[[2]])),rep("GENCODE_hg38_",length(input_list__lncRNAs[[2]])),rep("NONCODE",length(input_list_NONCODE[[2]])), rep("GENCODE_hg19",length(input_list_hg19_lncRNAs[[1]])),rep("GENCODE_hg38_v37",length(input_list_lncRNAs[[1]])),rep("GENCODE_hg38_",length(input_list__lncRNAs[[1]])),rep("NONCODE",length(input_list_NONCODE[[1]]))))))
colnames(all_sce_melted) <- c("ident","gene_id", "Assembly")

# Represent as %of the total number of features detected in each SCE object
pdf("upset_plot_lncRNAs_250counts_in_25_cells.pdf")
all_sce_melted_tbl=tbl_df(all_sce_melted)
all_sce_melted_tbl %>%
  group_by(gene_id,Assembly) %>%
  summarize(ident = list(ident)) %>%
  ggplot(aes(x = ident,fill=Assembly)) +
    geom_bar(aes(y = after_stat(comp_pct(count, PANEL, fill))), position="dodge")  +
    scale_x_upset(order_by = "degree", reverse = T) + theme_classic() + theme(plot.margin = margin(1,0.5,0.5,2, "cm")) + ylab("Percentage of detected lncRNAs for each assembly")
dev.off()

all_sce_melted <- as.data.frame(cbind(as.character(c(rep("kallisto",length(input_list_hg19_PCs[[2]])),rep("kallisto",length(input_list_PCs[[2]])),rep("kallisto",length(input_list_v45_PCs[[2]])),rep("CellRanger",length(input_list_hg19_PCs[[1]])),rep("CellRanger",length(input_list_PCs[[1]])),rep("CellRanger",length(input_list_v45_PCs[[1]])))),as.character(c(input_list_hg19_PCs[[2]],input_list_PCs[[2]], input_list_v45_PCs[[2]],input_list_hg19_PCs[[1]],input_list_PCs[[1]], input_list_v45_PCs[[1]])),as.character(c(rep("GENCODE_hg19",length(input_list_hg19_PCs[[2]])),rep("GENCODE_hg38_v37",length(input_list_PCs[[2]])),rep("GENCODE_hg38_",length(input_list_v45_PCs[[2]])), rep("GENCODE_hg19",length(input_list_hg19_PCs[[1]])),rep("GENCODE_hg38_v37",length(input_list_PCs[[1]])),rep("GENCODE_hg38_",length(input_list_v45_PCs[[1]]))))))
colnames(all_sce_melted) <- c("ident","gene_id", "Assembly")

# Represent as %of the total number of features detected in each SCE object
pdf("upset_plot_PCs_250counts_in_25_cells.pdf")
all_sce_melted_tbl=tbl_df(all_sce_melted)
all_sce_melted_tbl %>%
  group_by(gene_id,Assembly) %>%
  summarize(ident = list(ident)) %>%
  ggplot(aes(x = ident,fill=Assembly)) +
    geom_bar(aes(y = after_stat(comp_pct(count, PANEL, fill))), position="dodge")  +
    scale_x_upset(order_by = "degree", reverse = T) + theme_classic() + theme(plot.margin = margin(1,0.5,0.5,2, "cm")) + ylab("Percentage of detected protein-coding genes for each assembly")
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






