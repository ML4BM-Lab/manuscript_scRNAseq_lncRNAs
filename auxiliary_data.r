######################################################################################################################################
#################################  SEEKR analysis: Does AL121895.1 cluster with the cis-repressors? ##################################
#################  Also: Is any of the SEEKR communities mainly formed by kallisto-exclusive lncRNAs  ################################
######################################################################################################################################
library(stringr)
library(factoextra)
gencode_path <- "/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/gencode.v37.annotation.gtf"
hg38_ensembl_gtf <- as.data.frame(rtracklayer::import(gencode_path))

kmers <- read.table("/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode/6mers.csv",header=T,row.names=1,sep=",")
# cluster and correlate the samples
new_rownames <- str_match(rownames(kmers), "ENST\\s*(.*?)\\s*..............")[,1]
new_rownames <- gsub("\\|","",new_rownames)
new_rownames <- gsub("_","",new_rownames)
kmers <- kmers[-which(duplicated(new_rownames)),]
new_rownames <- str_match(rownames(kmers), "ENST\\s*(.*?)\\s*..............")[,1]
new_rownames <- gsub("\\|","",new_rownames)
new_rownames <- gsub("_","",new_rownames)
gencode_path <- "/home/egonie/data/egonie/phd/reference_genomes_kike/GRCh38/gencode/gencode.v37.annotation.gtf"
hg38_ensembl_gtf <- as.data.frame(rtracklayer::import(gencode_path))
transcript_name <- hg38_ensembl_gtf$transcript_name[match(new_rownames,hg38_ensembl_gtf$transcript_id)]
rownames(kmers) = transcript_name

filtered_names <- c("AL121895.1-204","SCAANT1-201","GNAS-AS1-201","TSIX-201", "SNHG14-201","BDNF-AS-201","CDKN2B-AS1-201","KCNQ1OT1-201","XIST-201","PVT1-201","DBET-201","PCAT6-201" ,"HOTAIRM1-201","HOTTIP-201","LINC00570-201","AC005828.4-201")
kmers_reduced <- kmers[filtered_names,]
# k-means clustering
set.seed(100000000)

pdf("seekr_all_candidates_in_calabrese_supp1.pdf")
print(fviz_cluster(kmeans(z_scores_corrected,centers = 8, iter.max = 50, nstart = 1) , data = z_scores_corrected,ellipse.type = "convex",repel = T,ggtheme = theme_classic()))
dev.off()


######################################################################################################################################
##############################################  Parse NONCODE annotation  ############################################################
######################################################################################################################################
# Part of the parsing is indicated in auxiliary_data.r while other is in auxiliary_data.sh
# Download from http://v5.noncode.org/download.php
# I need to add the geneid to the gtf
noncode_gtf <- as.data.frame(rtracklayer::import("~/dato-activo/reference.genomes_kike/GRCh38/NONCODE/NONCODEv5_human_hg38_lncRNA.gtf"))
complete_gtf <- noncode_gtf[1,]
gene_id <- unique(noncode_gtf$gene_id)
for (i in 87855:length(gene_id))
{
  print(i)
  positions <- grep(gene_id[i], noncode_gtf$gene_id)
  subset_gtf <- noncode_gtf[positions,]
  max_value <- max(subset_gtf$end)
  min_value <- min(subset_gtf$start)
  new_gene <- subset_gtf[1,]
  new_gene$end <- max_value
  new_gene$start <- min_value
  new_gene$width <- max_value - min_value
  new_gene[11:13] <- NA
  new_gene$type <- "gene"

  complete_gtf <- rbind(complete_gtf,rbind(new_gene, subset_gtf))
  
}
all_gtf <- complete_gtf
all_gtf_unique <- all_gtf[!duplicated(all_gtf), ]
a <- all_gtf_unique[1,]
all_gtf_unique[1,] <- all_gtf_unique[2,]
all_gtf_unique[2,] <- a
rtracklayer::export(makeGRangesFromDataFrame(all_gtf_unique, keep.extra.columns =T), "~/dato-activo/reference.genomes_kike/GRCh38/NONCODE/NONCODEv5_human_hg38_lncRNA_complete.gtf")

# Remaining steps are performed in auxiliary_data.sh

