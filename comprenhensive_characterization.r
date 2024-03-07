########################################################### Generic #################################################################################################
#####################################################################################################################################################################
libraries <- c("Seurat","patchwork", "tximport", "ggplot2", "scran" , "scater", "SingleCellExperiment", "BUSpaRse", "DropletUtils", "RColorBrewer", "cluster", "UpSetR","scRNAseq", "ggpubr","scDblFinder","reshape","ggupset","dplyr","readxl")

lapply(libraries, require, character.only = TRUE)
# Import gtf annotation & setWD, source functions
human_gencode_path <- "~/dato-activo/reference.genomes_kike/GRCh38/gencode/gencode.v37.annotation.gtf"
hg38_ensembl_gtf <- as.data.frame(rtracklayer::import(human_gencode_path))
hg38_ensembl_gtf$gene_id <- gsub("_","-",hg38_ensembl_gtf$gene_id)
lncrna_ens_ids_human <- unique(c(hg38_ensembl_gtf$gene_id[grep("lncRNA",hg38_ensembl_gtf$gene_type)]))
protein_coding_ens_ids_human <- unique(c(hg38_ensembl_gtf$gene_id[hg38_ensembl_gtf$gene_type=="protein_coding"]))
lncrna_names_human <- unique(hg38_ensembl_gtf$gene_name[hg38_ensembl_gtf$gene_id %in% lncrna_ens_ids_human])
protein_coding_names_human <-  unique(hg38_ensembl_gtf$gene_name[hg38_ensembl_gtf$gene_id %in% protein_coding_ens_ids_human])

mouse_gencode_path <- "/home/egonie/dato-activo/reference.genomes_kike/GRCm39/gencode/gencode.vM27.annotation.gtf"
mouse_gtf <- as.data.frame(rtracklayer::import(mouse_gencode_path))
mouse_gtf$gene_id <- gsub("_","-",mouse_gtf$gene_id)
lncrna_ens_ids_mouse <- unique(c(mouse_gtf$gene_id[grep("lncRNA",mouse_gtf$gene_type)]))
protein_coding_ens_ids_mouse <- unique(c(mouse_gtf$gene_id[mouse_gtf$gene_type=="protein_coding"]))
lncrna_names_mouse <- unique(mouse_gtf$gene_name[mouse_gtf$gene_id %in% lncrna_ens_ids_mouse])
protein_coding_names_mouse <-  unique(mouse_gtf$gene_name[mouse_gtf$gene_id %in% protein_coding_ens_ids_mouse])

source("/home/egonie/kike/phd/git_rep_hpclogin/manuscript_scRNAseq_lncRNAs/manuscript_functions.r")
setwd("/home/egonie/kike/phd/test_data/paper_figures/figure3_characterization_ex_kallisto")

# Datasets
human_10k_pbmc_path <- "/home/egonie/kike/phd/test_data/paper_figures/figure3_characterization_ex_kallisto/pbmc_datasets_updated.rds"
mouse_1k_brain_path <- "/home/egonie/kike/phd/test_data/paper_figures/figure3_characterization_ex_kallisto/mouse_datasets_completed.rds"
kallisto_figure2_path <- "/home/egonie/kike/phd/test_data/paper_figures/figure3_characterization_ex_kallisto/kallisto_datasets.rds"
cellRanger_figure2_path <- "/home/egonie/kike/phd/test_data/paper_figures/figure3_characterization_ex_kallisto/cellRanger_datasets.rds"

###########################################################################################################################################################
###########################################################################################################################################################
# This is generic data that only needs to be run once. After we will just load this information.
all_info_human <- exonic_length_transcript(hg38_ensembl_gtf)
exons = hg38_ensembl_gtf[hg38_ensembl_gtf$type=="exon",]
last_isoform <- exons[nrow(exons),]
all_info_human[nrow(all_info_human),]=c(last_isoform$gene_id,last_isoform$gene_name,last_isoform$transcript_id,last_isoform$width)

all_info_mouse <- exonic_length_transcript(mouse_gtf)
exons = mouse_gtf[mouse_gtf$type=="exon",]
last_isoform <- exons[nrow(exons),]
all_info_mouse[nrow(all_info_mouse),]=c(last_isoform$gene_id,last_isoform$gene_name,last_isoform$transcript_id,last_isoform$width)

saveRDS(all_info_human, "/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/transcrips_lengths.rds")
saveRDS(all_info_mouse, "/home/egonie/dato-activo/reference.genomes_kike/GRCm39/gencode/transcrips_lengths.rds")

all_info_human_path <- "/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/transcrips_lengths.rds"
all_info_mouse_path <- "/home/egonie/dato-activo/reference.genomes_kike/GRCm39/gencode/transcrips_lengths.rds"

longest_transcripts_human <- get_longest_transcript(all_info_human_path)
longest_transcripts_mouse <- get_longest_transcript(all_info_mouse_path)

saveRDS(longest_transcripts_human,"longest_transcripts_human.RDS")
longest_transcripts_human <- readRDS("longest_transcripts_human.RDS")
saveRDS(longest_transcripts_mouse,"longest_transcripts_mouse.RDS")
longest_transcripts_mouse <- readRDS("longest_transcripts_mouse.RDS")

# get for every transcript the list of coordinates of its exons
exons_longest_transcripts_human <- transcripts_exons_coordinates(human_gencode_path, longest_transcripts_human)
saveRDS(exons_longest_transcripts_human, "/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/exons_longest_transcripts.rds")
exons_longest_transcripts_human <- readRDS("/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/exons_longest_transcripts.rds")

exons_longest_transcripts_mouse <- transcripts_exons_coordinates(mouse_gencode_path, longest_transcripts_mouse)
saveRDS(exons_longest_transcripts_mouse, "/home/egonie/dato-activo/reference.genomes_kike/GRCm39/gencode/exons_longest_transcripts.rds")
exons_longest_transcripts_mouse <- readRDS("/home/egonie/dato-activo/reference.genomes_kike/GRCm39/gencode/exons_longest_transcripts.rds")

# Get the number of exons for the longest isoform of each gene
n_exons_human <- number_of_exons_longest_isoform(hg38_ensembl_gtf, longest_transcripts_human)
saveRDS(n_exons_human, "number_exons_human.RDS")
n_exons_human <- readRDS("number_exons_human.RDS")

n_exons_mouse <- number_of_exons_longest_isoform(mouse_gtf, longest_transcripts_mouse)
saveRDS(n_exons_mouse, "number_exons_mouse.RDS")
n_exons_mouse <- readRDS("number_exons_mouse.RDS")

# load crispr data
#CRISPR public data from paper https://www.science.org/doi/10.1126/science.aah7111 (~500 lncRNAs proved to participate in cell growth)
crispr_data <- readRDS("/home/egonie/kike/databases/hits_info_Liu_science_2015_ensids.rds")


############################################################################################################################
############ For every dataset generate the df_vp matrix with the info about the cell specificity and expression ###########
############################################################################################################################
threshold_minumun_gene_counts_v <- c(250,100,50,25)
threshold_cells_detected_v <- c(25,10,5,3)

############################################################################################################################
##################################################### Human 10k PBMCs ######################################################
############################################################################################################################
pbmc_datasets_updated <- readRDS(human_10k_pbmc_path)
cellRanger_sce_filt_clus <- pbmc_datasets_updated[["cellRanger"]]
STARsolo_sce_filt_clus <- pbmc_datasets_updated[["STARsolo"]]
kallisto_sce_filt_clus <- pbmc_datasets_updated[["kallisto"]]
alevin_sce_filt_clus <- pbmc_datasets_updated[["Alevin"]]