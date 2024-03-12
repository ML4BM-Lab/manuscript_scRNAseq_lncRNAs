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
