# Simulate scRNA-seq expression and then run Kallisto and Cell Ranger. Which pipeline detects more of the simulated lncRNAs? Which pipeline is more similar to the ground truth simulated count matrix
.libPaths( c( "/home/egonie/data/egonie/R/x86_64-pc-linux-gnu-library/4.1", .libPaths()) )

libraries <- c("Seurat","patchwork", "iSEE", "tximport", "ggplot2", "scran" , "scater", "SingleCellExperiment", "BUSpaRse", "DropletUtils", "RColorBrewer", "biomaRt", "bluster", "batchelor", "cluster", "RCurl", "SingleR", "pheatmap", "fossil", "UpSetR", "celldex", "Biostrings", "stringr", "scuttle","scRNAseq", "ggpubr","RBioinf")

lapply(libraries, require, character.only = TRUE)

gencode_path <- "/~/dato-activo/reference.genomes_kike/GRCh38/gencode/gencode.v37.annotation.gtf"
hg38_ensembl_gtf <- rtracklayer::import(gencode_path)
source("/home/egonie/kike/phd/git_rep_hpclogin/manuscript_scRNAseq_lncRNAs/manuscript_functions.r")

# 1. Modify reference
# Create 1 exon in the middle of the largest intron of each transcript.
start(hg38_ensembl_gtf)  <- start(hg38_ensembl_gtf)-100
end(hg38_ensembl_gtf) <- end(hg38_ensembl_gtf) + 100

transcripts_ids <- unique(hg38_ensembl_gtf$transcript_id)
transcripts_ids <- transcripts_ids[!is.na(transcripts_ids)]
exons <- hg38_ensembl_gtf[hg38_ensembl_gtf$type=="exon"]
transcripts <-  hg38_ensembl_gtf[hg38_ensembl_gtf$type=="transcript"]
new_isoforms <- as.data.frame(exons)[1,]
for (i in 1:length(transcripts_ids))
{
    print(i)
    subset_transcript <- transcripts[transcripts$transcript_id==transcripts_ids[i]]
    subset_exons <- exons[exons$transcript_id==transcripts_ids[i]]
    if (length(subset_exons)==1)
    {
        next
    }
    intronic_regions <- GenomicRanges::setdiff(subset_transcript, subset_exons)
    median_exonic_length <- median(width(subset_exons))
    if (median_exonic_length < max(width(intronic_regions)))
    {
        intronic_region_to_convert <- intronic_regions[which.max(width(intronic_regions))]
        aux <- round((width(intronic_region_to_convert)-median_exonic_length)/2)
        start(intronic_region_to_convert) <- start(intronic_region_to_convert) +aux 
        end(intronic_region_to_convert) <- end(intronic_region_to_convert) -aux
    }
    else # in that case we just add an exon which is all the intron minus +- 100 bp
    {
       aux <- round(max(width(intronic_regions))/3)
       intronic_region_to_convert <- intronic_regions[which.max(width(intronic_regions))]
       start(intronic_region_to_convert) <- start(intronic_region_to_convert) +aux 
       end(intronic_region_to_convert) <- end(intronic_region_to_convert) -aux
    }
    intronic_region_to_convert_df <- as.data.frame(intronic_region_to_convert)
    a <- as.data.frame(subset_exons)[1,]
    a[,1:5] <- intronic_region_to_convert_df
    a$exon_id <- paste("ENSE000000_NEW",i,sep="")
    new_isoforms <- rbind(new_isoforms, a)
}

# remove first row
added_isoforms_GR <- new_isoforms[-1,]
added_isoforms <- makeGRangesFromDataFrame(added_isoforms, keep.extra.columns=TRUE)
modified_gtf <- c(hg38_ensembl_gtf, added_isoforms_GR)
# Add +-100 bp
start(modified_gtf) <- start(modified_gtf) - 100
end(modified_gtf) <- end(modified_gtf) + 100
rtracklayer::export(modified_gtf, "modified_gtf_hg38_GENCODE_v37.gtf")

# 2. Get reads from the “fake GTF”. Create the fastq
#       ~3000 genes (50% lncRNAs)
#       ~100 cells (~70.000 reads per cell)
#       ~3000 genes per cell
#       23000 reads for each gene in each cell. In total 2300 reads for each gene. 2300*1500 = 3450000 Total reads


modified_gtf <- rtracklayer::import("/home/egonie/data/egonie/phd/test_data/Minju_ha_nat_comm/modified_gtf_hg38_GENCODE_v37.gtf")
lncrna_ens_ids <- unique(c(modified_gtf$gene_id[grep("lncRNA",modified_gtf$gene_type)]))
protein_coding_ens_ids <- unique(c(modified_gtf$gene_id[modified_gtf$gene_type=="protein_coding"]))
exons <- modified_gtf[modified_gtf$type == "exon",]
CB_10X <- read.table("3M-february-2018.txt")
library(Biostrings)
genome_fasta <- readDNAStringSet("/home/egonie/data/egonie/phd/reference_genomes_kike/GRCh38/gencode/GRCh38.primary_assembly_GENCODE.genome.fa")
names(genome_fasta) <- gsub(" .*","",names(genome_fasta))

simulate_reads(lncrna_ens_ids, modified_gtf, protein_coding_ens_ids, simu1.rds, CB_10X, simu1, genome_fasta)
simulate_reads(lncrna_ens_ids, modified_gtf, protein_coding_ens_ids, simu2.rds, CB_10X, simu2, genome_fasta)
simulate_reads(lncrna_ens_ids, modified_gtf, protein_coding_ens_ids, simu3.rds, CB_10X, simu3, genome_fasta)
simulate_reads(lncrna_ens_ids, modified_gtf, protein_coding_ens_ids, simu4.rds, CB_10X, simu4, genome_fasta)
simulate_reads(lncrna_ens_ids, modified_gtf, protein_coding_ens_ids, simu5.rds, CB_10X, simu5, genome_fasta)
simulate_reads(lncrna_ens_ids, modified_gtf, protein_coding_ens_ids, simu6.rds, CB_10X, simu6, genome_fasta)
simulate_reads(lncrna_ens_ids, modified_gtf, protein_coding_ens_ids, simu7.rds, CB_10X, simu7, genome_fasta)
simulate_reads(lncrna_ens_ids, modified_gtf, protein_coding_ens_ids, simu8.rds, CB_10X, simu8, genome_fasta)
simulate_reads(lncrna_ens_ids, modified_gtf, protein_coding_ens_ids, simu9.rds, CB_10X, simu9, genome_fasta)
simulate_reads(lncrna_ens_ids, modified_gtf, protein_coding_ens_ids, simu10.rds, CB_10X, simu10, genome_fasta)


# Cell Ranger and Kallisto are used to preprocessed the simulated fastq files using the unmodified reference annotation in the scrit simulate_expression.sh
setwd("/home/egonie/kike/phd/test_data/paper_figures")

all_genes_captured_kallisto <- c()
all_lncRNAs_captured_kallisto <- c()
all_PCs_captured_kallisto <- c()
all_lncRNAs_500_captured_kallisto <- c()
all_PCs_500_captured_kallisto <- c()
all_genes_captured_CellRanger <- c()
all_lncRNAs_captured_CellRanger <- c()
all_PCs_captured_CellRanger <- c()
all_lncRNAs_500_captured_CellRanger <- c()
all_PCs_500_captured_CellRanger <- c()

frob_distance_all_genes_kallisto <- c()
frob_distance_real_genes_kallisto <- c()
frob_distance_lncRNAs_kallisto <- c()
frob_distance_PCs_kallisto <- c()
frob_distance_all_genes_CellRanger <- c()
frob_distance_real_genes_CellRanger <- c()
frob_distance_lncRNAs_CellRanger <- c()
frob_distance_PCs_CellRanger <- c()

for (i in 1:10)
{
    print(i)
    cell_ranger_dir <- paste("/home/egonie/kike/phd/test_data/Minju_ha_nat_comm/simulation_scRNAseq/more_simulations/simulation",i,"_CellRanger/outs/filtered_feature_bc_matrix",sep="")
    cellRanger_data <- Read10X(data.dir = cell_ranger_dir, gene.column = 1)
    cellRanger <- CreateSeuratObject(cellRanger_data, project = "cellRanger")

    kallisto_dir <- paste("/home/egonie/kike/phd/test_data/Minju_ha_nat_comm/simulation_scRNAseq/more_simulations/simulation",i,"_kallisto/bustools_results",sep="")
    kallisto_data <- BUSpaRse::read_count_output(kallisto_dir, name = "cells_genes_NO_multimapping")
    kallisto <- CreateSeuratObject(kallisto_data, project = "kallisto")

    cellRanger <- as.SingleCellExperiment(cellRanger)
    colnames(cellRanger) <- gsub("\\-.*","",colnames(cellRanger))
    kallisto <- as.SingleCellExperiment(kallisto)

    # Real matrix
    real_expression <- readRDS(paste("simulation_scRNAseq/more_simulations/simulationsimu",i,"_real_count_matrix.rds",sep=""))
    universe_genes <- intersect(rownames(kallisto), rownames(real_expression))
    real_expression_sorted <- real_expression[universe_genes,colnames(kallisto)]
    kallisto_filtered <- kallisto[universe_genes,]
    cellRanger_filtered <- cellRanger[universe_genes,colnames(kallisto)]
    real_expressed_genes = rownames(real_expression_sorted)[rowSums(real_expression_sorted)>500]
    real_expressed_lncRNAs <- intersect(lncrna_ens_ids, real_expressed_genes)
    real_expressed_PCs <- intersect(protein_coding_ens_ids, real_expressed_genes)

    # we consider highly expressed genes as those with more than 500 counts since we have specifically enriched for the targeted genes
    threshold_minumun_gene_counts = 250 
    threshold_cells_detected = 0

    kallisto_top_genes <- top_genes(kallisto,threshold_minumun_gene_counts,threshold_cells_detected )
    cellRanger_top_genes <- top_genes(cellRanger,threshold_minumun_gene_counts,threshold_cells_detected )

    all_genes_captured_kallisto <- c(all_genes_captured_kallisto, 100*length(intersect(rownames(kallisto_top_genes), real_expressed_genes)) / length(real_expressed_genes))
    all_lncRNAs_captured_kallisto <- c(all_lncRNAs_captured_kallisto, 100 * length(intersect(rownames(kallisto_top_genes), real_expressed_lncRNAs)) / length(real_expressed_lncRNAs))
    all_PCs_captured_kallisto <- c(all_PCs_captured_kallisto, 100 * length(intersect(rownames(kallisto_top_genes), real_expressed_PCs)) / length(real_expressed_PCs))

    all_genes_captured_CellRanger <- c(all_genes_captured_CellRanger, 100 * length(intersect(rownames(cellRanger_top_genes), real_expressed_genes)) / length(real_expressed_genes))
    all_lncRNAs_captured_CellRanger <- c(all_lncRNAs_captured_CellRanger, 100* length(intersect(rownames(cellRanger_top_genes), real_expressed_lncRNAs)) / length(real_expressed_lncRNAs))
    all_PCs_captured_CellRanger <- c(all_PCs_captured_CellRanger, 100 * length(intersect(rownames(cellRanger_top_genes), real_expressed_PCs)) / length(real_expressed_PCs))

    threshold_minumun_gene_counts = 500 
    threshold_cells_detected = 0

    kallisto_top_genes <- top_genes(kallisto,threshold_minumun_gene_counts,threshold_cells_detected )
    cellRanger_top_genes <- top_genes(cellRanger,threshold_minumun_gene_counts,threshold_cells_detected )

    all_lncRNAs_500_captured_kallisto <- c(all_lncRNAs_500_captured_kallisto, 100 * length(intersect(rownames(kallisto_top_genes), real_expressed_lncRNAs)) / length(real_expressed_lncRNAs))
    all_PCs_500_captured_kallisto <- c(all_PCs_500_captured_kallisto, 100 * length(intersect(rownames(kallisto_top_genes), real_expressed_PCs)) / length(real_expressed_PCs))

    all_lncRNAs_500_captured_CellRanger <- c(all_lncRNAs_500_captured_CellRanger, 100* length(intersect(rownames(cellRanger_top_genes), real_expressed_lncRNAs)) / length(real_expressed_lncRNAs))
    all_PCs_500_captured_CellRanger <- c(all_PCs_500_captured_CellRanger, 100 * length(intersect(rownames(cellRanger_top_genes), real_expressed_PCs)) / length(real_expressed_PCs))


    test_match_order(colnames(counts(kallisto_filtered)),colnames(real_expression_sorted))
    test_match_order(colnames(counts(cellRanger_filtered)),colnames(real_expression_sorted))
    test_match_order(rownames(counts(kallisto_filtered)),rownames(real_expression_sorted))
    test_match_order(rownames(counts(cellRanger_filtered)),rownames(real_expression_sorted))

    frob_distance_all_genes_kallisto <- c(frob_distance_all_genes_kallisto, SMFilter::FDist2(as.matrix(counts(kallisto_filtered)), real_expression_sorted))
    # filtered matrix
    frob_distance_real_genes_kallisto <- c(frob_distance_real_genes_kallisto, SMFilter::FDist2(as.matrix(counts(kallisto_filtered)[real_expressed_genes,]), real_expression_sorted[real_expressed_genes,]))
    frob_distance_lncRNAs_kallisto <- c(frob_distance_lncRNAs_kallisto, SMFilter::FDist2(as.matrix(counts(kallisto_filtered)[real_expressed_lncRNAs,]), real_expression_sorted[real_expressed_lncRNAs,]))
    frob_distance_PCs_kallisto <- c(frob_distance_PCs_kallisto, SMFilter::FDist2(as.matrix(counts(kallisto_filtered)[real_expressed_PCs,]), real_expression_sorted[real_expressed_PCs,]))

    # Frobenius norm: CellRanger
    frob_distance_all_genes_CellRanger <- c(frob_distance_all_genes_CellRanger, SMFilter::FDist2(as.matrix(counts(cellRanger_filtered)), real_expression_sorted))
    frob_distance_real_genes_CellRanger <- c(frob_distance_real_genes_CellRanger, SMFilter::FDist2(as.matrix(counts(cellRanger_filtered)[real_expressed_genes,]), real_expression_sorted[real_expressed_genes,]))
    frob_distance_lncRNAs_CellRanger <- c(frob_distance_lncRNAs_CellRanger, SMFilter::FDist2(as.matrix(counts(cellRanger_filtered)[real_expressed_lncRNAs,]), real_expression_sorted[real_expressed_lncRNAs,]))
    frob_distance_PCs_CellRanger <- c(frob_distance_PCs_CellRanger, SMFilter::FDist2(as.matrix(counts(cellRanger_filtered)[real_expressed_PCs,]), real_expression_sorted[real_expressed_PCs,]))

}

pdf("simulations_lncRNAs_expression.pdf", width = 9)
df <- data.frame(Pipeline=c(rep("Cell Ranger",10),rep("Kallisto",10)),frob_distance=c((frob_distance_lncRNAs_CellRanger/mean(frob_distance_lncRNAs_CellRanger)), (frob_distance_lncRNAs_kallisto/mean(frob_distance_lncRNAs_CellRanger))))
p1 <- ggplot(df,aes(x=Pipeline, y=frob_distance)) + geom_boxplot() +xlab("Pipeline") + ylab("Frobenious distance (fold Cell Ranger)") + theme_classic() + theme(axis.text = element_text(size=13), axis.title = element_text(size=14)) + ggtitle("Frobenious distance between REAL count matrix and dataset: LncRNAs") + theme(plot.title = element_text(hjust = 0.5, size = 15)) + stat_compare_means(method = "t.test", label.x = 1.3) +ylim(0.5,1.03)
print(p1)

df <- data.frame(Pipeline=c(rep("Cell Ranger",10),rep("Kallisto",10)),detected_lncRNAs=c(all_lncRNAs_500_captured_CellRanger, all_lncRNAs_500_captured_kallisto))
p1 <- ggplot(df,aes(x=Pipeline, y=detected_lncRNAs)) + geom_boxplot() +xlab("Pipeline") +xlab("Pipeline") + ylab("Detected lncRNAs (%)") + theme_classic() + theme(axis.text = element_text(size=13), axis.title = element_text(size=14)) + ggtitle("Detected lncRNAs from the real ones (>500 counts)") + theme(plot.title = element_text(hjust = 0.5, size = 15)) + stat_compare_means(method = "t.test", label.x = 1.3)  +ylim(70,100)
print(p1)
dev.off()

