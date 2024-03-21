########################################################### Generic #################################################################################################
#####################################################################################################################################################################
libraries <- c("Seurat","patchwork", "tximport", "ggplot2", "scran" , "scater", "SingleCellExperiment", "BUSpaRse", "DropletUtils", "RColorBrewer", "cluster", "UpSetR","scRNAseq", "ggpubr","scDblFinder","reshape","ggupset","dplyr","readxl","stringr")

lapply(libraries, require, character.only = TRUE)
# Import gtf annotation & setWD, source functions
human_gencode_path <- "~/dato-activo/reference.genomes_kike/GRCh38/gencode/gencode.v37.annotation.gtf"
hg38_ensembl_gtf <- as.data.frame(rtracklayer::import(human_gencode_path))
hg38_ensembl_gtf$gene_id <- gsub("_","-",hg38_ensembl_gtf$gene_id)
lncrna_ens_ids_human <- unique(c(hg38_ensembl_gtf$gene_id[grep("lncRNA",hg38_ensembl_gtf$gene_type)]))
protein_coding_ens_ids_human <- unique(c(hg38_ensembl_gtf$gene_id[hg38_ensembl_gtf$gene_type=="protein_coding"]))
lncrna_names_human <- unique(hg38_ensembl_gtf$gene_name[hg38_ensembl_gtf$gene_id %in% lncrna_ens_ids_human])
protein_coding_names_human <-  unique(hg38_ensembl_gtf$gene_name[hg38_ensembl_gtf$gene_id %in% protein_coding_ens_ids_human])
human_repeatMasker_cleaned_path <- "/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/GRCh38.primary_assembly_GENCODE.genome.fa.out_cleaned.gff"
seekr_6_communities_human_path <- "/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/SEEKR_communities_6mers.csv"

mouse_gencode_path <- "/home/egonie/dato-activo/reference.genomes_kike/GRCm39/gencode/gencode.vM27.annotation.gtf"
mouse_gtf <- as.data.frame(rtracklayer::import(mouse_gencode_path))
mouse_gtf$gene_id <- gsub("_","-",mouse_gtf$gene_id)
lncrna_ens_ids_mouse <- unique(c(mouse_gtf$gene_id[grep("lncRNA",mouse_gtf$gene_type)]))
protein_coding_ens_ids_mouse <- unique(c(mouse_gtf$gene_id[mouse_gtf$gene_type=="protein_coding"]))
lncrna_names_mouse <- unique(mouse_gtf$gene_name[mouse_gtf$gene_id %in% lncrna_ens_ids_mouse])
protein_coding_names_mouse <-  unique(mouse_gtf$gene_name[mouse_gtf$gene_id %in% protein_coding_ens_ids_mouse])
mouse_repeatMasker_cleaned_path <- "/home/egonie/dato-activo/reference.genomes_kike/GRCm39/gencode/GRCm39.primary_assembly.genome.fa.out_cleaned.gff"
seekr_6_communities_mouse_path <- "/home/egonie/dato-activo/reference.genomes_kike/GRCm39/gencode/SEEKR_communities_6mers.csv"

source("/home/egonie/kike/phd/git_rep_hpclogin/manuscript_scRNAseq_lncRNAs/manuscript_functions.r")
source("/home/egonie/kike/phd/git_rep_hpclogin/manuscript_scRNAseq_lncRNAs/characterization_function.r")
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

# Load repeats from repeatMasker
human_repeats_cleaned <- rtracklayer::import(human_repeatMasker_cleaned_path)
mouse_repeats_cleaned <- rtracklayer::import(mouse_repeatMasker_cleaned_path)

# load SEEKR data: For analyzing function of lncRNAs according to k-mer content
seekr_6_communities_human <- load_SEEKR_communities(seekr_6_communities_human_path, hg38_ensembl_gtf)
seekr_6_communities_mouse <- load_SEEKR_communities(seekr_6_communities_mouse_path, mouse_gtf)

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

#get universe of genes
universe_genes <- intersect(rownames(kallisto_sce_filt_clus),rownames(alevin_sce_filt_clus))
cellRanger_sce_filt_clus <- cellRanger_sce_filt_clus[universe_genes,]
STARsolo_sce_filt_clus <- STARsolo_sce_filt_clus[universe_genes,]
kallisto_sce_filt_clus <- kallisto_sce_filt_clus[universe_genes,]
alevin_sce_filt_clus <- alevin_sce_filt_clus[universe_genes,]

# Gene length distribution:
length_PBMCs <- length_distributions(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_sce_filt_clus, cellRanger_sce_filt_clus, STARsolo_sce_filt_clus, alevin_sce_filt_clus, lncrna_names_human, protein_coding_names_human,longest_transcripts=longest_transcripts_human,gene_name="gene_name")

# Number of exons
number_exons_PBMCs <- number_exons_distributions(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_sce_filt_clus, cellRanger_sce_filt_clus, STARsolo_sce_filt_clus, alevin_sce_filt_clus, lncrna_names_human, protein_coding_names_human,n_exons_all=n_exons_human,gene_name="gene_name")

# Repeat content
final_repeats_percentage_PBMCs <- repeats_results(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_sce_filt_clus, cellRanger_sce_filt_clus, STARsolo_sce_filt_clus, alevin_sce_filt_clus,lncrna_names_human,protein_coding_names_human,hg38_repeats=human_repeats_cleaned, exons_longest_transcript = exons_longest_transcripts_human)
ratio_repeats_hg_10k_PBMCs <- ratios_repeats(final_repeats_percentage_PBMCs,"hg_10k_PBMCs")

# K-mer analysis (SEEKR)
seekr_results_hg_10k_PBMCs <- SEEKR_results(threshold_minumun_gene_counts_v, threshold_cells_detected_v, kallisto_sce_filt_clus, cellRanger_sce_filt_clus, STARsolo_sce_filt_clus, alevin_sce_filt_clus, seekr_6_communities_human,lncrnas_ids = lncrna_names_human )

# Specificity index (SI)
df_vp <- create_df_vp(kallisto_sce_filt_clus,cellRanger_sce_filt_clus,STARsolo_sce_filt_clus=STARsolo_sce_filt_clus,alevin_sce_filt_clus=alevin_sce_filt_clus, lncrna_names_human, protein_coding_names_human )
saveRDS(df_vp,"df_vp_PBMC.rds")
df_vp_PBMC <- readRDS("df_vp_PBMC.rds")

# Intersect with bibliographically validated lncRNAs by CRISPRi
all_crispr_data_intersection_PBMCs <- crispr_data_intersection(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_sce_filt_clus, cellRanger_sce_filt_clus, STARsolo_sce_filt_clus, alevin_sce_filt_clus, lncrna_names_human,crispr_data,gene_name="gene_name",hg38_ensembl_gtf)


############################################################################################################################
################################################## 1k Mouse Brain cells ####################################################
############################################################################################################################
mouse_datasets_updated <- readRDS(mouse_1k_brain_path)
cellRanger_sce_filt_clus <- mouse_datasets_updated[["cellRanger"]]
STARsolo_sce_filt_clus <- mouse_datasets_updated[["STARsolo"]]
kallisto_sce_filt_clus <- mouse_datasets_updated[["kallisto"]]
alevin_sce_filt_clus <- mouse_datasets_updated[["Alevin"]] 

#get universe of genes
universe_genes <- intersect(rownames(kallisto_sce_filt_clus),rownames(alevin_sce_filt_clus))
cellRanger_sce_filt_clus <- cellRanger_sce_filt_clus[universe_genes,]
STARsolo_sce_filt_clus <- STARsolo_sce_filt_clus[universe_genes,]
kallisto_sce_filt_clus <- kallisto_sce_filt_clus[universe_genes,]
alevin_sce_filt_clus <- alevin_sce_filt_clus[universe_genes,]

# Gene length distribution:
length_Mouse_Brain <- length_distributions(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_sce_filt_clus, cellRanger_sce_filt_clus, STARsolo_sce_filt_clus, alevin_sce_filt_clus, lncrna_names_mouse, protein_coding_names_mouse,longest_transcripts=longest_transcripts_mouse,gene_name="gene_name")

# Number of exons
number_exons_Mouse_Brain <- number_exons_distributions(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_sce_filt_clus, cellRanger_sce_filt_clus, STARsolo_sce_filt_clus, alevin_sce_filt_clus, lncrna_names_mouse, protein_coding_names_mouse,n_exons_all=n_exons_mouse,gene_name="gene_name")

# Repeat content
final_repeats_percentage_Mouse_Brain <- repeats_results(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_sce_filt_clus, cellRanger_sce_filt_clus, STARsolo_sce_filt_clus, alevin_sce_filt_clus,lncrna_names_mouse,protein_coding_names_mouse,hg38_repeats=mouse_repeats_cleaned, exons_longest_transcript = exons_longest_transcripts_mouse)
ratio_repeats_mm_1k_brain <- ratios_repeats(final_repeats_percentage_Mouse_Brain,"mm_1k_brain")

# K-mer analysis (SEEKR)
seekr_results_mm_1k_brain <- SEEKR_results(threshold_minumun_gene_counts_v, threshold_cells_detected_v, kallisto_sce_filt_clus, cellRanger_sce_filt_clus, STARsolo_sce_filt_clus, alevin_sce_filt_clus, seekr_6_communities_mouse,lncrnas_ids = lncrna_names_mouse )

# Specificity index (SI)
df_vp_Mouse_Brain <- create_df_vp(kallisto_sce_filt_clus,cellRanger_sce_filt_clus,STARsolo_sce_filt_clus=STARsolo_sce_filt_clus,alevin_sce_filt_clus=alevin_sce_filt_clus, lncrna_names_mouse, protein_coding_names_mouse )
saveRDS(df_vp_Mouse_Brain,"df_vp_Mouse_Brain.rds")
df_vp_Mouse_Brain <- readRDS("df_vp_Mouse_Brain.rds")


##############################################################################################################################################################
###################################################  Datasets extended benchmark  ############################################################################
##############################################################################################################################################################
kallisto_objects <- readRDS(kallisto_figure2_path)
kallisto_intestine_pool1_ed <- kallisto_objects[["intestine1"]]
kallisto_intestine_pool2_ed <- kallisto_objects[["intestine2"]]
kallisto_healthy_lung_ed <- kallisto_objects[["healthy_lung1"]]
kallisto_healthy_lung_GSM4037316_ed <- kallisto_objects[["healthy_lung2"]]
kallisto_pulmonary_fibrosis_ed <- kallisto_objects[["pulmonary_fibrosis"]]
kallisto_PBMCs_5K_ed <- kallisto_objects[["PBMCs_5K"]]
kallisto_PBMCs_mouse_ed <- kallisto_objects[["PBMCs_mouse"]]

cellRanger_objects <- readRDS(cellRanger_figure2_path)
cellRanger_intestine_pool1_ed <- cellRanger_objects[["intestine1"]]
cellRanger_intestine_pool2_ed <- cellRanger_objects[["intestine2"]]
cellRanger_healthy_lung_ed <- cellRanger_objects[["healthy_lung1"]]
cellRanger_healthy_lung_GSM4037316_ed <- cellRanger_objects[["healthy_lung2"]]
cellRanger_pulmonary_fibrosis_ed <- cellRanger_objects[["pulmonary_fibrosis"]]
cellRanger_PBMCs_5K_ed <- cellRanger_objects[["PBMCs_5K"]]
cellRanger_PBMCs_mouse_ed <- cellRanger_objects[["PBMCs_mouse"]]

# Filtering & normalization
threshold_mito_percentage = 15
high_threshold_cell_counts = 50000
cells_min_genes_detected_threshold = 500

kallisto_intestine_pool1_ed_filt <- Filtering(kallisto_intestine_pool1_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
kallisto_intestine_pool2_ed_filt <- Filtering(kallisto_intestine_pool2_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
kallisto_healthy_lung_ed_filt <- Filtering(kallisto_healthy_lung_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
kallisto_healthy_lung_GSM4037316_ed_filt <- Filtering(kallisto_healthy_lung_GSM4037316_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
kallisto_pulmonary_fibrosis_ed_filt <- Filtering(kallisto_pulmonary_fibrosis_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
kallisto_PBMCs_5K_ed_filt <- Filtering(kallisto_PBMCs_5K_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
kallisto_PBMCs_mouse_ed_filt <- Filtering(kallisto_PBMCs_mouse_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)

cellRanger_intestine_pool1_ed_filt <- Filtering(cellRanger_intestine_pool1_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
cellRanger_intestine_pool2_ed_filt <- Filtering(cellRanger_intestine_pool2_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
cellRanger_healthy_lung_ed_filt <- Filtering(cellRanger_healthy_lung_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
cellRanger_healthy_lung_GSM4037316_ed_filt <- Filtering(cellRanger_healthy_lung_GSM4037316_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
cellRanger_pulmonary_fibrosis_ed_filt <- Filtering(cellRanger_pulmonary_fibrosis_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
cellRanger_PBMCs_5K_ed_filt <- Filtering(cellRanger_PBMCs_5K_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)
cellRanger_PBMCs_mouse_ed_filt <- Filtering(cellRanger_PBMCs_mouse_ed,  cells_mito_threshold = threshold_mito_percentage, cells_max_threshold = high_threshold_cell_counts, cells_min_genes_detected_threshold = cells_min_genes_detected_threshold)

#gene names features
kallisto_intestine_pool1_ed_filt <- gene_names_sce(kallisto_intestine_pool1_ed_filt, hg38_ensembl_gtf)
kallisto_intestine_pool2_ed_filt <- gene_names_sce(kallisto_intestine_pool2_ed_filt, hg38_ensembl_gtf)
kallisto_healthy_lung_ed_filt <- gene_names_sce(kallisto_healthy_lung_ed_filt, hg38_ensembl_gtf)
kallisto_healthy_lung_GSM4037316_ed_filt <- gene_names_sce(kallisto_healthy_lung_GSM4037316_ed_filt, hg38_ensembl_gtf)
kallisto_pulmonary_fibrosis_ed_filt <- gene_names_sce(kallisto_pulmonary_fibrosis_ed_filt, hg38_ensembl_gtf)
kallisto_PBMCs_5K_ed_filt <- gene_names_sce(kallisto_PBMCs_5K_ed_filt, hg38_ensembl_gtf)
kallisto_PBMCs_mouse_ed_filt <- gene_names_sce(kallisto_PBMCs_mouse_ed_filt, mouse_gtf)

cellRanger_intestine_pool1_ed_filt <- gene_names_sce(cellRanger_intestine_pool1_ed_filt, hg38_ensembl_gtf)
cellRanger_intestine_pool2_ed_filt <- gene_names_sce(cellRanger_intestine_pool2_ed_filt, hg38_ensembl_gtf)
cellRanger_healthy_lung_ed_filt <- gene_names_sce(cellRanger_healthy_lung_ed_filt, hg38_ensembl_gtf)
cellRanger_healthy_lung_GSM4037316_ed_filt <- gene_names_sce(cellRanger_healthy_lung_GSM4037316_ed_filt, hg38_ensembl_gtf)
cellRanger_pulmonary_fibrosis_ed_filt <- gene_names_sce(cellRanger_pulmonary_fibrosis_ed_filt, hg38_ensembl_gtf)
cellRanger_PBMCs_5K_ed_filt <- gene_names_sce(cellRanger_PBMCs_5K_ed_filt, hg38_ensembl_gtf)
cellRanger_PBMCs_mouse_ed_filt <- gene_names_sce(cellRanger_PBMCs_mouse_ed_filt, mouse_gtf)

################################################## Intestine pool 1 ##################################################################################
# Gene length distribution:
length_intestine_pool1 <- length_distributions(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_intestine_pool1_ed_filt, cellRanger_intestine_pool1_ed_filt, cellRanger_intestine_pool1_ed_filt, cellRanger_intestine_pool1_ed_filt, lncrna_names_human, protein_coding_names_human,longest_transcripts=longest_transcripts_human,gene_name="gene_name")

# Number of exons
number_exons_intestine_pool1 <- number_exons_distributions(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_intestine_pool1_ed_filt, cellRanger_intestine_pool1_ed_filt, cellRanger_intestine_pool1_ed_filt, cellRanger_intestine_pool1_ed_filt, lncrna_names_human, protein_coding_names_human,n_exons_all=n_exons_human,gene_name="gene_name")

# Repeat content
final_repeats_percentage_intestine_pool1 <- repeats_results(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_intestine_pool1_ed_filt, cellRanger_intestine_pool1_ed_filt, cellRanger_intestine_pool1_ed_filt, cellRanger_intestine_pool1_ed_filt,lncrna_names_human,protein_coding_names_human,hg38_repeats=human_repeats_cleaned, exons_longest_transcript = exons_longest_transcripts_human)
ratio_repeats_intestine_pool1 <- ratios_repeats(final_repeats_percentage_intestine_pool1,"hg_intestine_1")

# K-mer analysis (SEEKR)
seekr_results_intestine_pool1 <- SEEKR_results(threshold_minumun_gene_counts_v, threshold_cells_detected_v, kallisto_intestine_pool1_ed_filt, cellRanger_intestine_pool1_ed_filt, cellRanger_intestine_pool1_ed_filt, cellRanger_intestine_pool1_ed_filt, seekr_6_communities_human,lncrnas_ids = lncrna_names_human )

# Specificity index (SI)
kallisto_intestine_pool1_ed_filt <- red_dim(kallisto_intestine_pool1_ed_filt)
df_vp_intestine_pool1 <- create_df_vp(kallisto_intestine_pool1_ed_filt,cellRanger_intestine_pool1_ed_filt,STARsolo_sce_filt_clus=cellRanger_intestine_pool1_ed_filt,alevin_sce_filt_clus=cellRanger_intestine_pool1_ed_filt, lncrna_names_human, protein_coding_names_human )
saveRDS(df_vp_intestine_pool1,"df_vp_intestine_pool1.rds")
df_vp_intestine_pool1 <- readRDS("df_vp_intestine_pool1.rds")

# Intersect with bibliographically validated lncRNAs by CRISPRi
all_crispr_data_intersection_intestine_pool1 <- crispr_data_intersection(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_intestine_pool1_ed_filt, cellRanger_intestine_pool1_ed_filt, cellRanger_intestine_pool1_ed_filt, cellRanger_intestine_pool1_ed_filt, lncrna_names_human,crispr_data,gene_name="gene_name",hg38_ensembl_gtf)

################################################## Intestine pool 2 ##################################################################################
# Gene length distribution:
length_intestine_pool2 <- length_distributions(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_intestine_pool2_ed_filt, cellRanger_intestine_pool2_ed_filt, cellRanger_intestine_pool2_ed_filt, cellRanger_intestine_pool2_ed_filt, lncrna_names_human, protein_coding_names_human,longest_transcripts=longest_transcripts_human,gene_name="gene_name")

# Number of exons
number_exons_intestine_pool2 <- number_exons_distributions(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_intestine_pool2_ed_filt, cellRanger_intestine_pool2_ed_filt, cellRanger_intestine_pool2_ed_filt, cellRanger_intestine_pool2_ed_filt, lncrna_names_human, protein_coding_names_human,n_exons_all=n_exons_human,gene_name="gene_name")

# Repeat content
final_repeats_percentage_intestine_pool2 <- repeats_results(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_intestine_pool2_ed_filt, cellRanger_intestine_pool2_ed_filt, cellRanger_intestine_pool2_ed_filt, cellRanger_intestine_pool2_ed_filt,lncrna_names_human,protein_coding_names_human,hg38_repeats=human_repeats_cleaned, exons_longest_transcript = exons_longest_transcripts_human)
ratio_repeats_intestine_pool2 <- ratios_repeats(final_repeats_percentage_intestine_pool2,"hg_intestine_2")

# K-mer analysis (SEEKR)
seekr_results_intestine_pool2 <- SEEKR_results(threshold_minumun_gene_counts_v, threshold_cells_detected_v, kallisto_intestine_pool2_ed_filt, cellRanger_intestine_pool2_ed_filt, cellRanger_intestine_pool2_ed_filt, cellRanger_intestine_pool2_ed_filt, seekr_6_communities_human,lncrnas_ids = lncrna_names_human )

# Specificity index (SI)
kallisto_intestine_pool2_ed_filt <- red_dim(kallisto_intestine_pool2_ed_filt)
df_vp_intestine_pool2 <- create_df_vp(kallisto_intestine_pool2_ed_filt,cellRanger_intestine_pool2_ed_filt,STARsolo_sce_filt_clus=cellRanger_intestine_pool2_ed_filt,alevin_sce_filt_clus=cellRanger_intestine_pool2_ed_filt, lncrna_names_human, protein_coding_names_human )
saveRDS(df_vp_intestine_pool2,"df_vp_intestine_pool2.rds")
df_vp_intestine_pool2 <- readRDS("df_vp_intestine_pool2.rds")

# Intersect with bibliographically validated lncRNAs by CRISPRi
all_crispr_data_intersection_intestine_pool2 <- crispr_data_intersection(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_intestine_pool2_ed_filt, cellRanger_intestine_pool2_ed_filt, cellRanger_intestine_pool2_ed_filt, cellRanger_intestine_pool2_ed_filt, lncrna_names_human,crispr_data,gene_name="gene_name",hg38_ensembl_gtf)

################################################## Healthy lung ##################################################################################
# Gene length distribution:
length_healthy_lung <- length_distributions(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_healthy_lung_ed_filt, cellRanger_healthy_lung_ed_filt, cellRanger_healthy_lung_ed_filt, cellRanger_healthy_lung_ed_filt, lncrna_names_human, protein_coding_names_human,longest_transcripts=longest_transcripts_human,gene_name="gene_name")

# Number of exons
number_exons_healthy_lung <- number_exons_distributions(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_healthy_lung_ed_filt, cellRanger_healthy_lung_ed_filt ,cellRanger_healthy_lung_ed_filt, cellRanger_healthy_lung_ed_filt, lncrna_names_human, protein_coding_names_human,n_exons_all=n_exons_human,gene_name="gene_name")

# Repeat content
final_repeats_percentage_healthy_lung <- repeats_results(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_healthy_lung_ed_filt, cellRanger_healthy_lung_ed_filt, cellRanger_healthy_lung_ed_filt, cellRanger_healthy_lung_ed_filt,lncrna_names_human,protein_coding_names_human,hg38_repeats=human_repeats_cleaned, exons_longest_transcript = exons_longest_transcripts_human)
ratio_repeats_healthy_lung <- ratios_repeats(final_repeats_percentage_healthy_lung,"hg_lung_healthy_1")

# K-mer analysis (SEEKR)
seekr_results_healthy_lung <- SEEKR_results(threshold_minumun_gene_counts_v, threshold_cells_detected_v, kallisto_healthy_lung_ed_filt, cellRanger_healthy_lung_ed_filt, cellRanger_healthy_lung_ed_filt, cellRanger_healthy_lung_ed_filt,  seekr_6_communities_human,lncrnas_ids = lncrna_names_human )

# Specificity index (SI)
kallisto_healthy_lung_ed_filt <- red_dim(kallisto_healthy_lung_ed_filt)
df_vp_healthy_lung <- create_df_vp(kallisto_healthy_lung_ed_filt,cellRanger_healthy_lung_ed_filt,STARsolo_sce_filt_clus=cellRanger_healthy_lung_ed_filt,alevin_sce_filt_clus=cellRanger_healthy_lung_ed_filt, lncrna_names_human, protein_coding_names_human )
saveRDS(df_vp_healthy_lung,"df_vp_healthy_lung.rds")
df_vp_healthy_lung <- readRDS("df_vp_healthy_lung.rds")

# Intersect with bibliographically validated lncRNAs by CRISPRi
all_crispr_data_intersection_healthy_lung <- crispr_data_intersection(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_healthy_lung_ed_filt, cellRanger_healthy_lung_ed_filt, cellRanger_healthy_lung_ed_filt, cellRanger_healthy_lung_ed_filt,  lncrna_names_human,crispr_data,gene_name="gene_name",hg38_ensembl_gtf)

################################################## Healthy lung GSM4037316 ##################################################################################
# Gene length distribution:
length_healthy_lung_GSM4037316 <- length_distributions(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_healthy_lung_GSM4037316_ed_filt, cellRanger_healthy_lung_GSM4037316_ed_filt, cellRanger_healthy_lung_GSM4037316_ed_filt, cellRanger_healthy_lung_GSM4037316_ed_filt,  lncrna_names_human, protein_coding_names_human,longest_transcripts=longest_transcripts_human,gene_name="gene_name")

# Number of exons
number_exons_healthy_lung_GSM4037316 <- number_exons_distributions(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_healthy_lung_GSM4037316_ed_filt, cellRanger_healthy_lung_GSM4037316_ed_filt, cellRanger_healthy_lung_GSM4037316_ed_filt, cellRanger_healthy_lung_GSM4037316_ed_filt, lncrna_names_human, protein_coding_names_human,n_exons_all=n_exons_human,gene_name="gene_name")

# Repeat content
final_repeats_percentage_healthy_lung_GSM4037316 <- repeats_results(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_healthy_lung_GSM4037316_ed_filt, cellRanger_healthy_lung_GSM4037316_ed_filt, cellRanger_healthy_lung_GSM4037316_ed_filt, cellRanger_healthy_lung_GSM4037316_ed_filt,lncrna_names_human,protein_coding_names_human,hg38_repeats=human_repeats_cleaned, exons_longest_transcript = exons_longest_transcripts_human)
ratio_repeats_healthy_lung_GSM4037316 <- ratios_repeats(final_repeats_percentage_healthy_lung_GSM4037316,"hg_lung_healthy_2")

# K-mer analysis (SEEKR)
seekr_results_healthy_lung_GSM4037316 <- SEEKR_results(threshold_minumun_gene_counts_v, threshold_cells_detected_v, kallisto_healthy_lung_GSM4037316_ed_filt, cellRanger_healthy_lung_GSM4037316_ed_filt, cellRanger_healthy_lung_GSM4037316_ed_filt, cellRanger_healthy_lung_GSM4037316_ed_filt, seekr_6_communities_human,lncrnas_ids = lncrna_names_human )

# Specificity index (SI)
kallisto_healthy_lung_GSM4037316_ed_filt <- red_dim(kallisto_healthy_lung_GSM4037316_ed_filt)
df_vp_healthy_lung_GSM4037316 <- create_df_vp(kallisto_healthy_lung_GSM4037316_ed_filt,cellRanger_healthy_lung_GSM4037316_ed_filt,STARsolo_sce_filt_clus=cellRanger_healthy_lung_GSM4037316_ed_filt,alevin_sce_filt_clus=cellRanger_healthy_lung_GSM4037316_ed_filt,  lncrna_names_human, protein_coding_names_human )
saveRDS(df_vp_healthy_lung_GSM4037316,"df_vp_healthy_lung_GSM4037316.rds")
df_vp_healthy_lung_GSM4037316 <- readRDS("df_vp_healthy_lung_GSM4037316.rds")

# Intersect with bibliographically validated lncRNAs by CRISPRi
all_crispr_data_intersection_healthy_lung_GSM4037316 <- crispr_data_intersection(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_healthy_lung_GSM4037316_ed_filt, cellRanger_healthy_lung_GSM4037316_ed_filt, cellRanger_healthy_lung_GSM4037316_ed_filt, cellRanger_healthy_lung_GSM4037316_ed_filt, lncrna_names_human,crispr_data,gene_name="gene_name",hg38_ensembl_gtf)

################################################## Pulmonary fibrosis ##################################################################################
# Gene length distribution:
length_pulmonary_fibrosis <- length_distributions(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_pulmonary_fibrosis_ed_filt, cellRanger_pulmonary_fibrosis_ed_filt, cellRanger_pulmonary_fibrosis_ed_filt, cellRanger_pulmonary_fibrosis_ed_filt,  lncrna_names_human, protein_coding_names_human,longest_transcripts=longest_transcripts_human,gene_name="gene_name")

# Number of exons
number_exons_pulmonary_fibrosis <- number_exons_distributions(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_pulmonary_fibrosis_ed_filt, cellRanger_pulmonary_fibrosis_ed_filt, cellRanger_pulmonary_fibrosis_ed_filt, cellRanger_pulmonary_fibrosis_ed_filt, lncrna_names_human, protein_coding_names_human,n_exons_all=n_exons_human,gene_name="gene_name")

# Repeat content
final_repeats_percentage_pulmonary_fibrosis <- repeats_results(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_pulmonary_fibrosis_ed_filt, cellRanger_pulmonary_fibrosis_ed_filt, cellRanger_pulmonary_fibrosis_ed_filt, cellRanger_pulmonary_fibrosis_ed_filt,lncrna_names_human,protein_coding_names_human,hg38_repeats=human_repeats_cleaned, exons_longest_transcript = exons_longest_transcripts_human)
ratio_repeats_pulmonary_fibrosis <- ratios_repeats(final_repeats_percentage_pulmonary_fibrosis,"hg_lung_fibrosis")

# K-mer analysis (SEEKR)
seekr_results_pulmonary_fibrosis <- SEEKR_results(threshold_minumun_gene_counts_v, threshold_cells_detected_v, kallisto_pulmonary_fibrosis_ed_filt, cellRanger_pulmonary_fibrosis_ed_filt, cellRanger_pulmonary_fibrosis_ed_filt, cellRanger_pulmonary_fibrosis_ed_filt, seekr_6_communities_human,lncrnas_ids = lncrna_names_human )

# Specificity index (SI)
kallisto_pulmonary_fibrosis_ed_filt <- red_dim(kallisto_pulmonary_fibrosis_ed_filt)
df_vp_pulmonary_fibrosis <- create_df_vp(kallisto_pulmonary_fibrosis_ed_filt,cellRanger_pulmonary_fibrosis_ed_filt,STARsolo_sce_filt_clus=cellRanger_pulmonary_fibrosis_ed_filt,alevin_sce_filt_clus=cellRanger_pulmonary_fibrosis_ed_filt,  lncrna_names_human, protein_coding_names_human )
saveRDS(df_vp_pulmonary_fibrosis,"df_vp_pulmonary_fibrosis.rds")
df_vp_pulmonary_fibrosis <- readRDS("df_vp_pulmonary_fibrosis.rds")

# Intersect with bibliographically validated lncRNAs by CRISPRi
all_crispr_data_intersection_pulmonary_fibrosis <- crispr_data_intersection(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_pulmonary_fibrosis_ed_filt, cellRanger_pulmonary_fibrosis_ed_filt, cellRanger_pulmonary_fibrosis_ed_filt, cellRanger_pulmonary_fibrosis_ed_filt, lncrna_names_human,crispr_data,gene_name="gene_name",hg38_ensembl_gtf)

################################################## Human 5K PBMCs  ##################################################################################
# Gene length distribution:
length_PBMCs_5K <- length_distributions(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_PBMCs_5K_ed_filt, cellRanger_PBMCs_5K_ed_filt, cellRanger_PBMCs_5K_ed_filt, cellRanger_PBMCs_5K_ed_filt,  lncrna_names_human, protein_coding_names_human,longest_transcripts=longest_transcripts_human,gene_name="gene_name")

# Number of exons
number_exons_PBMCs_5K <- number_exons_distributions(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_PBMCs_5K_ed_filt, cellRanger_PBMCs_5K_ed_filt, cellRanger_PBMCs_5K_ed_filt, cellRanger_PBMCs_5K_ed_filt, lncrna_names_human, protein_coding_names_human,n_exons_all=n_exons_human,gene_name="gene_name")

# Repeat content
final_repeats_percentage_PBMCs_5K <- repeats_results(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_PBMCs_5K_ed_filt, cellRanger_PBMCs_5K_ed_filt, cellRanger_PBMCs_5K_ed_filt, cellRanger_PBMCs_5K_ed_filt,lncrna_names_human,protein_coding_names_human,hg38_repeats=human_repeats_cleaned, exons_longest_transcript = exons_longest_transcripts_human)
ratio_repeats_PBMCs_5K <- ratios_repeats(final_repeats_percentage_PBMCs_5K,"hg_5k_PBMCs")

# K-mer analysis (SEEKR)
seekr_results_PBMCs_5K <- SEEKR_results(threshold_minumun_gene_counts_v, threshold_cells_detected_v, kallisto_PBMCs_5K_ed_filt, cellRanger_PBMCs_5K_ed_filt, cellRanger_PBMCs_5K_ed_filt, cellRanger_PBMCs_5K_ed_filt,  seekr_6_communities_human,lncrnas_ids = lncrna_names_human )

# Specificity index (SI)
kallisto_PBMCs_5K_ed_filt <- red_dim(kallisto_PBMCs_5K_ed_filt)
df_vp_PBMCs_5K <- create_df_vp(kallisto_PBMCs_5K_ed_filt,cellRanger_PBMCs_5K_ed_filt,STARsolo_sce_filt_clus=cellRanger_PBMCs_5K_ed_filt,alevin_sce_filt_clus=cellRanger_PBMCs_5K_ed_filt,  lncrna_names_human, protein_coding_names_human )
saveRDS(df_vp_PBMCs_5K,"df_vp_PBMCs_5K.rds")
df_vp_PBMCs_5K <- readRDS("df_vp_PBMCs_5K.rds")

# Intersect with bibliographically validated lncRNAs by CRISPRi
all_crispr_data_intersection_PBMCs_5K <- crispr_data_intersection(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_PBMCs_5K_ed_filt, cellRanger_PBMCs_5K_ed_filt, cellRanger_PBMCs_5K_ed_filt, cellRanger_PBMCs_5K_ed_filt, lncrna_names_human,crispr_data,gene_name="gene_name",hg38_ensembl_gtf)


################################################## Mouse PBMCs  ##################################################################################
# Gene length distribution:
length_PBMCs_mouse <- length_distributions(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_PBMCs_mouse_ed_filt, cellRanger_PBMCs_mouse_ed_filt, cellRanger_PBMCs_mouse_ed_filt, cellRanger_PBMCs_mouse_ed_filt, lncrna_names_mouse, protein_coding_names_mouse,longest_transcripts=longest_transcripts_mouse,gene_name="gene_name")

# Number of exons
number_exons_PBMCs_mouse <- number_exons_distributions(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_PBMCs_mouse_ed_filt, cellRanger_PBMCs_mouse_ed_filt, cellRanger_PBMCs_mouse_ed_filt, cellRanger_PBMCs_mouse_ed_filt, lncrna_names_mouse, protein_coding_names_mouse,n_exons_all=n_exons_mouse,gene_name="gene_name")

# Repeat content
final_repeats_percentage_PBMCs_mouse <- repeats_results(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_PBMCs_mouse_ed_filt, cellRanger_PBMCs_mouse_ed_filt, cellRanger_PBMCs_mouse_ed_filt, cellRanger_PBMCs_mouse_ed_filt,lncrna_names_mouse,protein_coding_names_mouse,hg38_repeats=mouse_repeats_cleaned, exons_longest_transcript = exons_longest_transcripts_mouse)
ratio_repeats_PBMCs_mouse <- ratios_repeats(final_repeats_percentage_PBMCs_mouse,"mm_10k_PBMCs")

# K-mer analysis (SEEKR)
seekr_results_PBMCs_mouse <- SEEKR_results(threshold_minumun_gene_counts_v, threshold_cells_detected_v, kallisto_PBMCs_mouse_ed_filt, cellRanger_PBMCs_mouse_ed_filt, cellRanger_PBMCs_mouse_ed_filt, cellRanger_PBMCs_mouse_ed_filt, seekr_6_communities_mouse,lncrnas_ids = lncrna_names_mouse )

# Specificity index (SI)
kallisto_PBMCs_mouse_ed_filt <- red_dim(kallisto_PBMCs_mouse_ed_filt)
df_vp_PBMCs_mouse <- create_df_vp(kallisto_PBMCs_mouse_ed_filt,cellRanger_PBMCs_mouse_ed_filt,STARsolo_sce_filt_clus=cellRanger_PBMCs_mouse_ed_filt,alevin_sce_filt_clus=cellRanger_PBMCs_mouse_ed_filt, lncrna_names_mouse, protein_coding_names_mouse )
saveRDS(df_vp_PBMCs_mouse,"df_vp_PBMCs_mouse.rds")
df_vp_PBMCs_mouse <- readRDS("df_vp_PBMCs_mouse.rds")


##############################################################################################################################################################
#read all rds (9 in total, 2 in fig1, 7 in fig2)
df_vp_PBMC <- readRDS("df_vp_PBMC.rds")
df_vp_Mouse_Brain <- readRDS("df_vp_Mouse_Brain.rds")
df_vp_intestine_pool1 <- readRDS("df_vp_intestine_pool1.rds")
df_vp_intestine_pool2 <- readRDS("df_vp_intestine_pool2.rds")
df_vp_healthy_lung <- readRDS("df_vp_healthy_lung.rds")
df_vp_healthy_lung_GSM4037316 <- readRDS("df_vp_healthy_lung_GSM4037316.rds")
df_vp_pulmonary_fibrosis <- readRDS("df_vp_pulmonary_fibrosis.rds")
df_vp_PBMCs_5K <- readRDS("df_vp_PBMCs_5K.rds")
df_vp_PBMCs_mouse <- readRDS("df_vp_PBMCs_mouse.rds")

##############################################################################################################################################################
# 1. Expression 
all_datasets_250 <- rbind(main_expression_object(df_vp_PBMC, "hg_10k_PBMCs", 250),main_expression_object(df_vp_Mouse_Brain, "mm_1k_brain", 250),main_expression_object(df_vp_intestine_pool1, "hg_intestine_1", 250),main_expression_object(df_vp_intestine_pool2, "hg_intestine_2", 250),main_expression_object(df_vp_healthy_lung, "hg_lung_healthy_1", 250),main_expression_object(df_vp_healthy_lung_GSM4037316, "hg_lung_healthy_2", 250),main_expression_object(df_vp_pulmonary_fibrosis, "hg_lung_fibrosis", 250),main_expression_object(df_vp_PBMCs_5K, "hg_5k_PBMCs", 250),main_expression_object(df_vp_PBMCs_mouse, "mm_10k_PBMCs", 250))
pdf("violin_plot_expression_t250.pdf",width = 22, height = 8)
violin_plot_expression(all_datasets_250)
dev.off()

all_datasets_100 <- rbind(main_expression_object(df_vp_PBMC, "hg_10k_PBMCs", 100),main_expression_object(df_vp_Mouse_Brain, "mm_1k_brain", 100),main_expression_object(df_vp_intestine_pool1, "hg_intestine_1", 100),main_expression_object(df_vp_intestine_pool2, "hg_intestine_2", 100),main_expression_object(df_vp_healthy_lung, "hg_lung_healthy_1", 100),main_expression_object(df_vp_healthy_lung_GSM4037316, "hg_lung_healthy_2", 100),main_expression_object(df_vp_pulmonary_fibrosis, "hg_lung_fibrosis", 100),main_expression_object(df_vp_PBMCs_5K, "hg_5k_PBMCs", 100),main_expression_object(df_vp_PBMCs_mouse, "mm_10k_PBMCs", 100))
pdf("violin_plot_expression_t100.pdf",width = 22, height = 8)
violin_plot_expression(all_datasets_100)
dev.off()

all_datasets_50 <- rbind(main_expression_object(df_vp_PBMC, "hg_10k_PBMCs", 50),main_expression_object(df_vp_Mouse_Brain, "mm_1k_brain", 50),main_expression_object(df_vp_intestine_pool1, "hg_intestine_1", 50),main_expression_object(df_vp_intestine_pool2, "hg_intestine_2", 50),main_expression_object(df_vp_healthy_lung, "hg_lung_healthy_1", 50),main_expression_object(df_vp_healthy_lung_GSM4037316, "hg_lung_healthy_2", 50),main_expression_object(df_vp_pulmonary_fibrosis, "hg_lung_fibrosis", 50),main_expression_object(df_vp_PBMCs_5K, "hg_5k_PBMCs", 50),main_expression_object(df_vp_PBMCs_mouse, "mm_10k_PBMCs", 50))
pdf("violin_plot_expression_t50.pdf",width = 22, height = 8)
violin_plot_expression(all_datasets_50)
dev.off()

all_datasets_25 <- rbind(main_expression_object(df_vp_PBMC, "hg_10k_PBMCs", 25),main_expression_object(df_vp_Mouse_Brain, "mm_1k_brain", 25),main_expression_object(df_vp_intestine_pool1, "hg_intestine_1", 25),main_expression_object(df_vp_intestine_pool2, "hg_intestine_2", 25),main_expression_object(df_vp_healthy_lung, "hg_lung_healthy_1", 25),main_expression_object(df_vp_healthy_lung_GSM4037316, "hg_lung_healthy_2", 25),main_expression_object(df_vp_pulmonary_fibrosis, "hg_lung_fibrosis", 25),main_expression_object(df_vp_PBMCs_5K, "hg_5k_PBMCs", 25),main_expression_object(df_vp_PBMCs_mouse, "mm_10k_PBMCs", 25))
pdf("violin_plot_expression_t25.pdf",width = 22, height = 8)
violin_plot_expression(all_datasets_25)
dev.off()


##############################################################################################################################################################
# 2. Length
options(scipen=999)
all_datasets_250 <- rbind(main_object(length_PBMCs, "hg_10k_PBMCs", 250),main_object(length_Mouse_Brain, "mm_1k_brain", 250),main_object(length_intestine_pool1, "hg_intestine_1", 250),main_object(length_intestine_pool2, "hg_intestine_2", 250),main_object(length_healthy_lung, "hg_lung_healthy_1", 250),main_object(length_healthy_lung_GSM4037316, "hg_lung_healthy_2", 250),main_object(length_pulmonary_fibrosis, "hg_lung_fibrosis", 250),main_object(length_PBMCs_5K, "hg_5k_PBMCs", 250),main_object(length_PBMCs_mouse, "mm_10k_PBMCs", 250))
pdf("violin_plot_length_t250.pdf",width = 22, height = 8)
violin_plot_length(all_datasets_250)
dev.off()

all_datasets_100 <- rbind(main_object(length_PBMCs, "hg_10k_PBMCs", 100),main_object(length_Mouse_Brain, "mm_1k_brain", 100),main_object(length_intestine_pool1, "hg_intestine_1", 100),main_object(length_intestine_pool2, "hg_intestine_2", 100),main_object(length_healthy_lung, "hg_lung_healthy_1", 100),main_object(length_healthy_lung_GSM4037316, "hg_lung_healthy_2", 100),main_object(length_pulmonary_fibrosis, "hg_lung_fibrosis", 100),main_object(length_PBMCs_5K, "hg_5k_PBMCs", 100),main_object(length_PBMCs_mouse, "mm_10k_PBMCs", 100))
pdf("violin_plot_length_t100.pdf",width = 22, height = 8)
violin_plot_length(all_datasets_100)
dev.off()

all_datasets_50 <- rbind(main_object(length_PBMCs, "hg_10k_PBMCs", 50),main_object(length_Mouse_Brain, "mm_1k_brain", 50),main_object(length_intestine_pool1, "hg_intestine_1", 50),main_object(length_intestine_pool2, "hg_intestine_2", 50),main_object(length_healthy_lung, "hg_lung_healthy_1", 50),main_object(length_healthy_lung_GSM4037316, "hg_lung_healthy_2", 50),main_object(length_pulmonary_fibrosis, "hg_lung_fibrosis", 50),main_object(length_PBMCs_5K, "hg_5k_PBMCs", 50),main_object(length_PBMCs_mouse, "mm_10k_PBMCs", 50))
pdf("violin_plot_length_t50.pdf",width = 22, height = 8)
violin_plot_length(all_datasets_50)
dev.off()

all_datasets_25 <- rbind(main_object(length_PBMCs, "hg_10k_PBMCs", 25),main_object(length_Mouse_Brain, "mm_1k_brain", 25),main_object(length_intestine_pool1, "hg_intestine_1", 25),main_object(length_intestine_pool2, "hg_intestine_2", 25),main_object(length_healthy_lung, "hg_lung_healthy_1", 25),main_object(length_healthy_lung_GSM4037316, "hg_lung_healthy_2", 25),main_object(length_pulmonary_fibrosis, "hg_lung_fibrosis", 25),main_object(length_PBMCs_5K, "hg_5k_PBMCs", 25),main_object(length_PBMCs_mouse, "mm_10k_PBMCs", 25))
pdf("violin_plot_length_t25.pdf",width = 22, height = 8)
violin_plot_length(all_datasets_25)
dev.off()

##############################################################################################################################################################
# 3. Number of exons
all_datasets_250 <- rbind(main_object(number_exons_PBMCs, "hg_10k_PBMCs", 250),main_object(number_exons_Mouse_Brain, "mm_1k_brain", 250),main_object(number_exons_intestine_pool1, "hg_intestine_1", 250),main_object(number_exons_intestine_pool2, "hg_intestine_2", 250),main_object(number_exons_healthy_lung, "hg_lung_healthy_1", 250),main_object(number_exons_healthy_lung_GSM4037316, "hg_lung_healthy_2", 250),main_object(number_exons_pulmonary_fibrosis, "hg_lung_fibrosis", 250),main_object(number_exons_PBMCs_5K, "hg_5k_PBMCs", 250),main_object(number_exons_PBMCs_mouse, "mm_10k_PBMCs", 250))
pdf("violin_plot_number_exons_t250.pdf",width = 22, height = 8)
violin_plot_number_exons(all_datasets_250)
dev.off()

all_datasets_100 <- rbind(main_object(number_exons_PBMCs, "hg_10k_PBMCs", 100),main_object(number_exons_Mouse_Brain, "mm_1k_brain", 100),main_object(number_exons_intestine_pool1, "hg_intestine_1", 100),main_object(number_exons_intestine_pool2, "hg_intestine_2", 100),main_object(number_exons_healthy_lung, "hg_lung_healthy_1", 100),main_object(number_exons_healthy_lung_GSM4037316, "hg_lung_healthy_2", 100),main_object(number_exons_pulmonary_fibrosis, "hg_lung_fibrosis", 100),main_object(number_exons_PBMCs_5K, "hg_5k_PBMCs", 100),main_object(number_exons_PBMCs_mouse, "mm_10k_PBMCs", 100))
pdf("violin_plot_number_exons_t100.pdf",width = 22, height = 8)
violin_plot_number_exons(all_datasets_100)
dev.off()

all_datasets_50 <- rbind(main_object(number_exons_PBMCs, "hg_10k_PBMCs", 50),main_object(number_exons_Mouse_Brain, "mm_1k_brain", 50),main_object(number_exons_intestine_pool1, "hg_intestine_1", 50),main_object(number_exons_intestine_pool2, "hg_intestine_2", 50),main_object(number_exons_healthy_lung, "hg_lung_healthy_1", 50),main_object(number_exons_healthy_lung_GSM4037316, "hg_lung_healthy_2", 50),main_object(number_exons_pulmonary_fibrosis, "hg_lung_fibrosis", 50),main_object(number_exons_PBMCs_5K, "hg_5k_PBMCs", 50),main_object(number_exons_PBMCs_mouse, "mm_10k_PBMCs", 50))
pdf("violin_plot_number_exons_t50.pdf",width = 22, height = 8)
violin_plot_number_exons(all_datasets_50)
dev.off()

all_datasets_25 <- rbind(main_object(number_exons_PBMCs, "hg_10k_PBMCs", 25),main_object(number_exons_Mouse_Brain, "mm_1k_brain", 25),main_object(number_exons_intestine_pool1, "hg_intestine_1", 25),main_object(number_exons_intestine_pool2, "hg_intestine_2", 25),main_object(number_exons_healthy_lung, "hg_lung_healthy_1", 25),main_object(number_exons_healthy_lung_GSM4037316, "hg_lung_healthy_2", 25),main_object(number_exons_pulmonary_fibrosis, "hg_lung_fibrosis", 25),main_object(number_exons_PBMCs_5K, "hg_5k_PBMCs", 25),main_object(number_exons_PBMCs_mouse, "mm_10k_PBMCs", 25))
pdf("violin_plot_number_exons_t25.pdf",width = 22, height = 8)
violin_plot_number_exons(all_datasets_25)
dev.off()


##############################################################################################################################################################
# 4. Repeat content
repeat_ratio_df <- rbind(ratio_repeats_hg_10k_PBMCs, ratio_repeats_mm_1k_brain, ratio_repeats_intestine_pool1, ratio_repeats_intestine_pool2, ratio_repeats_healthy_lung, ratio_repeats_healthy_lung_GSM4037316, ratio_repeats_pulmonary_fibrosis, ratio_repeats_PBMCs_5K, ratio_repeats_PBMCs_mouse )
repeat_ratio_df$Threshold <- factor(repeat_ratio_df$Threshold, levels = c("25","50", "100","250"))
repeat_ratio_df$dataset <- factor(repeat_ratio_df$dataset, levels = c("hg_10k_PBMCs","mm_1k_brain", "hg_lung_healthy_1","hg_lung_healthy_2","hg_5k_PBMCs","hg_intestine_1","hg_intestine_2","mm_10k_PBMCs","hg_lung_fibrosis"))
#change dataset labels to match those of fig2
levels(repeat_ratio_df$dataset)[levels(repeat_ratio_df$dataset)=='hg_lung_healthy_1'] <- 'Healthy_lung_GSM5020383'
levels(repeat_ratio_df$dataset)[levels(repeat_ratio_df$dataset)=='hg_lung_healthy_2'] <- 'Healthy_lung_GSM4037316'
levels(repeat_ratio_df$dataset)[levels(repeat_ratio_df$dataset)=='hg_5k_PBMCs'] <- 'Human_5k_PBMCs'
levels(repeat_ratio_df$dataset)[levels(repeat_ratio_df$dataset)=='hg_10k_PBMCs'] <- 'Human_10k_PBMCs'
levels(repeat_ratio_df$dataset)[levels(repeat_ratio_df$dataset)=='hg_intestine_1'] <- 'Intestine_GSM4808339'
levels(repeat_ratio_df$dataset)[levels(repeat_ratio_df$dataset)=='hg_intestine_2'] <- 'Intestine_GSM4808348'
levels(repeat_ratio_df$dataset)[levels(repeat_ratio_df$dataset)=='mm_10k_PBMCs'] <- 'Mouse_10k_PBMCs'
levels(repeat_ratio_df$dataset)[levels(repeat_ratio_df$dataset)=='mm_1k_brain'] <- 'Mouse_1k_Brain'
levels(repeat_ratio_df$dataset)[levels(repeat_ratio_df$dataset)=='hg_lung_fibrosis'] <- 'Pulmonary_fibrosis_GSM4037320'
repeat_ratio_df$dataset <- factor(repeat_ratio_df$dataset, levels = c("Healthy_lung_GSM5020383","Healthy_lung_GSM4037316", "Human_5k_PBMCs","Human_10k_PBMCs","Intestine_GSM4808339","Intestine_GSM4808348","Mouse_1k_Brain","Mouse_10k_PBMCs","Pulmonary_fibrosis_GSM4037320"))

pdf("Repeat_content_all_lncRNAs.pdf")
p1 <- ggplot(repeat_ratio_df,aes(Threshold, y = LncRNAs, fill = dataset,color=dataset)) +geom_jitter(width=0.3, height=0,size=3) +theme_classic() +ggtitle("LncRNAs") +ylab("Ratio of repeat content exclusive vs common") + xlab("Defined threshold of expression") +theme(axis.text=element_text(size=14), axis.title=element_text(size=15), plot.title = element_text(size = 16)) +ylim(0,3.5) +coord_flip()
p1
dev.off()

pdf("Repeat_content_all_PCs.pdf")
p2 <- ggplot(repeat_ratio_df,aes(Threshold, y = Protein_coding, fill = dataset,color = dataset)) +geom_jitter(width=0.3, height=0,size=3) +theme_classic() +ggtitle("Protein-coding genes") +ylab("Ratio of repeat content exclusive vs common") + xlab("Defined threshold of expression") +theme(axis.text=element_text(size=14), axis.title=element_text(size=15), plot.title = element_text(size = 16))+ylim(0,3.5) +coord_flip()
p2 
dev.off()
t.test(repeat_ratio_df$Protein_coding, repeat_ratio_df$LncRNAs, alternetive = "greater")

##############################################################################################################################################################
# 5. K-mer content (SEEKR)
all_datasets_250 <- rbind(main_seekr_object(seekr_results_hg_10k_PBMCs, "hg_10k_PBMCs", 250),main_seekr_object(seekr_results_mm_1k_brain, "mm_1k_brain", 250),main_seekr_object(seekr_results_intestine_pool1, "hg_intestine_1", 250),main_seekr_object(seekr_results_intestine_pool2, "hg_intestine_2", 250),main_seekr_object(seekr_results_healthy_lung, "hg_lung_healthy_1", 250),main_seekr_object(seekr_results_healthy_lung_GSM4037316, "hg_lung_healthy_2", 250),main_seekr_object(seekr_results_pulmonary_fibrosis, "hg_lung_fibrosis", 250),main_seekr_object(seekr_results_PBMCs_5K, "hg_5k_PBMCs", 250),main_seekr_object(seekr_results_PBMCs_mouse, "mm_10k_PBMCs", 250))
pdf("seekr_barplot_t250.pdf",width = 22, height = 8)
seekr_barplot(all_datasets_250)
dev.off()

t=100
all_datasets_100 <- rbind(main_seekr_object(seekr_results_hg_10k_PBMCs, "hg_10k_PBMCs", t),main_seekr_object(seekr_results_mm_1k_brain, "mm_1k_brain", t),main_seekr_object(seekr_results_intestine_pool1, "hg_intestine_1", t),main_seekr_object(seekr_results_intestine_pool2, "hg_intestine_2", t),main_seekr_object(seekr_results_healthy_lung, "hg_lung_healthy_1", t),main_seekr_object(seekr_results_healthy_lung_GSM4037316, "hg_lung_healthy_2", t),main_seekr_object(seekr_results_pulmonary_fibrosis, "hg_lung_fibrosis", t),main_seekr_object(seekr_results_PBMCs_5K, "hg_5k_PBMCs", t),main_seekr_object(seekr_results_PBMCs_mouse, "mm_10k_PBMCs", t))
pdf("seekr_barplot_t100.pdf",width = 22, height = 8)
seekr_barplot(all_datasets_100)
dev.off()

t=50
all_datasets_50 <- rbind(main_seekr_object(seekr_results_hg_10k_PBMCs, "hg_10k_PBMCs", t),main_seekr_object(seekr_results_mm_1k_brain, "mm_1k_brain", t),main_seekr_object(seekr_results_intestine_pool1, "hg_intestine_1", t),main_seekr_object(seekr_results_intestine_pool2, "hg_intestine_2", t),main_seekr_object(seekr_results_healthy_lung, "hg_lung_healthy_1", t),main_seekr_object(seekr_results_healthy_lung_GSM4037316, "hg_lung_healthy_2", t),main_seekr_object(seekr_results_pulmonary_fibrosis, "hg_lung_fibrosis", t),main_seekr_object(seekr_results_PBMCs_5K, "hg_5k_PBMCs", t),main_seekr_object(seekr_results_PBMCs_mouse, "mm_10k_PBMCs", t))
pdf("seekr_barplot_t50.pdf",width = 22, height = 8)
seekr_barplot(all_datasets_50)
dev.off()

t=25
all_datasets_25 <- rbind(main_seekr_object(seekr_results_hg_10k_PBMCs, "hg_10k_PBMCs", t),main_seekr_object(seekr_results_mm_1k_brain, "mm_1k_brain", t),main_seekr_object(seekr_results_intestine_pool1, "hg_intestine_1", t),main_seekr_object(seekr_results_intestine_pool2, "hg_intestine_2", t),main_seekr_object(seekr_results_healthy_lung, "hg_lung_healthy_1", t),main_seekr_object(seekr_results_healthy_lung_GSM4037316, "hg_lung_healthy_2", t),main_seekr_object(seekr_results_pulmonary_fibrosis, "hg_lung_fibrosis", t),main_seekr_object(seekr_results_PBMCs_5K, "hg_5k_PBMCs", t),main_seekr_object(seekr_results_PBMCs_mouse, "mm_10k_PBMCs", t))
pdf("seekr_barplot_t25.pdf",width = 22, height = 8)
seekr_barplot(all_datasets_25)
dev.off()

##############################################################################################################################################################
# 6. SI of exclusive lncRNAs vs all protein-coding genes
t=250
all_datasets_250 <- rbind(main_SI_object(df_vp_PBMC, "hg_10k_PBMCs", t),main_SI_object(df_vp_Mouse_Brain, "mm_1k_brain", t),main_SI_object(df_vp_intestine_pool1, "hg_intestine_1", t),main_SI_object(df_vp_intestine_pool2, "hg_intestine_2", t),main_SI_object(df_vp_healthy_lung, "hg_lung_healthy_1", t),main_SI_object(df_vp_healthy_lung_GSM4037316, "hg_lung_healthy_2", t),main_SI_object(df_vp_pulmonary_fibrosis, "hg_lung_fibrosis", t),main_SI_object(df_vp_PBMCs_5K, "hg_5k_PBMCs", t),main_SI_object(df_vp_PBMCs_mouse, "mm_10k_PBMCs", t))
pdf("violin_plot_SI_ob2_ALLpcs_vs_EXCLUSIVElncRNAs_t250_less.pdf",width = 18, height = 10)
violin_plot_SI_ob2(all_datasets_250,hypothesis="less")
dev.off()

pdf("violin_plot_SI_ob2_ALLpcs_vs_EXCLUSIVElncRNAs_t250_greater.pdf",width = 18, height = 10)
violin_plot_SI_ob2(all_datasets_250,hypothesis="greater")
dev.off()

##############################################################################################################################################################
# 7. SI of exclusive lncRNAs vs all protein-coding genes
options(scipen=0)
all_datasets <- rbind(counts_crispr(all_crispr_data_intersection_PBMCs, "hg_10k_PBMCs"),counts_crispr(all_crispr_data_intersection_intestine_pool1, "hg_intestine_1"),counts_crispr(all_crispr_data_intersection_intestine_pool2, "hg_intestine_2"),counts_crispr(all_crispr_data_intersection_healthy_lung, "hg_lung_healthy_1"),counts_crispr(all_crispr_data_intersection_healthy_lung_GSM4037316, "hg_lung_healthy_2"),counts_crispr(all_crispr_data_intersection_pulmonary_fibrosis, "hg_lung_fibrosis"),counts_crispr(all_crispr_data_intersection_PBMCs_5K, "hg_5k_PBMCs"))
supp_table <- rbind(generate_supp_table(all_crispr_data_intersection_PBMCs, "hg_10k_PBMCs"),generate_supp_table(all_crispr_data_intersection_intestine_pool1, "hg_intestine_1"),generate_supp_table(all_crispr_data_intersection_intestine_pool2, "hg_intestine_2"),generate_supp_table(all_crispr_data_intersection_healthy_lung, "hg_lung_healthy_1"),generate_supp_table(all_crispr_data_intersection_healthy_lung_GSM4037316, "hg_lung_healthy_2"),generate_supp_table(all_crispr_data_intersection_pulmonary_fibrosis, "hg_lung_fibrosis"),generate_supp_table(all_crispr_data_intersection_PBMCs_5K, "hg_5k_PBMCs"))
write.table(supp_table, "intersection_lncRNAs_crispr_Liu_et_al.txt", row.names = F, quote = F, sep = "\t" )
write.table(crispr_data, "crispr_data_Liu_et_al.txt", row.names = F, quote = F, sep = "\t" ) # not sure if I need to add it as it is already published

pdf("barplot_intersection_crispr_Liu_et_al.pdf",width = 20, height = 7)
#my_df <- df_vp[df_vp$k_clus==25,]
colors <- c("#D4B996FF","#A07855FF")
a = all_datasets
a$Threshold <- as.factor(a$Threshold)
ggplot(a[a$features=="exclusive_kallisto",],aes(x = dataset, y = crispr_genes,fill = Threshold))  + geom_bar(stat="identity",width = 0.8,position=position_dodge())  + geom_text(aes(label=p_value_hypergeom), position = position_dodge(width = .8),vjust = -1)+ theme_classic() + theme(axis.text = element_text(size=13), axis.title.x = element_text(size=14))+ theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm')) + theme( strip.background = element_blank(),strip.placement = "outside",strip.text.x = element_text(angle = 45,size = 12) ) +scale_fill_manual(values=c("#2E5266FF","#6E8898FF","#9FB1BCFF","#D3D0CBFF")) + ggtitle("Exclusive lncRNAs found by Kallisto")
dev.off()
