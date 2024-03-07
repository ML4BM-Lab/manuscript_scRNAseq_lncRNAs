############################################################################################################################
########################### Bibliographic investigation to identify functionally validated lncRNA ##########################
############################################################################################################################
#1. CRISPR public data from paper https://www.science.org/doi/10.1126/science.aah7111 (~500 lncRNAs proved to participate in cell growth)
# Download supplementary tables 1 (metadata information about hits) and supplementary tables 4 (hits information)
library("readxl")
crispr_screening_hits <- as.data.frame(read_excel("/home/egonie/kike/databases/aah7111-tables4.xlsx"))
crispr_screening_info <- as.data.frame(read_excel("/home/egonie/kike/databases/aah7111-tables1.xlsx"))
#lncRNAs hits
lncRNAs_hits <- c()
for (i in which(crispr_screening_hits[3,]=="lncRNA gene hit type"))
{
  if (i==101)
  {
    print("skipping")
  }
  else
  {
    print(i)
    lncRNAs_hits <- c(lncRNAs_hits,crispr_screening_hits[,1][which(crispr_screening_hits[,i]=="lncRNA hit")])
  }
}
unique_hits <- unique(lncRNAs_hits)
hits_info <- crispr_screening_info[crispr_screening_info[,1] %in% unique_hits,1:4]
write.table(hits_info, "hits_info_Liu_science_2015.txt",row.names = F, col.names = T)

# Get geneID from lncRNAs hits (not working so far) (authors of the Liu.et.al manuscript used gencode version 25)
crispr_data <- read.table("/home/egonie/kike/databases/hits_info_Liu_science_2015.txt", header=T)
hg38_gtf_v25 <- rtracklayer::import("/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode_release_25/gencode.v25.annotation.gtf.gz")
crispr_data$is_annotated_in_gencode <- crispr_data$gene.name %in% hg38_gtf_v25$gene_name
crispr_ens_ids <- c()
for (i in 1:nrow(crispr_data))
{
    if (crispr_data$is_annotated_in_gencode[i])
    {
        crispr_ens_ids <- c(crispr_ens_ids, unique(hg38_gtf_v25$gene_id[which(hg38_gtf_v25$gene_name==crispr_data$gene.name[i])])[1])
    }
    else
    {
        crispr_ens_ids <- c(crispr_ens_ids, NA)
    }
}
crispr_data$gene_id <- crispr_ens_ids
crispr_data$Transcript.ID[which(crispr_data$gene.name=="-")]
crispr_data$gene_id[which(crispr_data$gene.name=="-")]="-"
 # Some of them were annotated in genodev28 
hg38_gtf_v28 <- rtracklayer::import("/home/egonie/dato-activo/reference.genomes_kike/GRCh38/gencode_release_28/gencode.v28.annotation.gtf")
a <- which(crispr_data$gene.name[which(is.na(crispr_data$gene_id))] %in% hg38_gtf_v28$gene_name)
crispr_ens_ids_v28 <- c()
for (j in 1:length(crispr_data$gene.name[which(is.na(crispr_data$gene_id))]))
{
    if (j %in% a)
    {
        b <- unique(hg38_gtf_v28$gene_id[which(hg38_gtf_v28$gene_name==crispr_data$gene.name[which(is.na(crispr_data$gene_id))][j])])[1]
        crispr_ens_ids_v28 <- c(crispr_ens_ids_v28,b)
    }
}
crispr_data_complete <- crispr_data[which(is.na(crispr_data$gene_id)),]
crispr_data_complete$gene_id[a] <- crispr_ens_ids_v28
 # I have 34 genes that I've tried to manually annotate by bibliographic inspection (genome browser, genecards...) 
manually_ens_id_crispr <- c("-","ENSG00000226500.2","ENSG00000230623.4","ENSG00000234546.3", "-","ENSG00000234546.3", "ENSG00000235823.1","-","-","-","ENSG00000238045.9","ENSG00000175061.17","ENSG00000175061.17","ENSG00000267131.1","-","-","ENSG00000269107.1","-","-","-","-","-","-","-","-","ENSG00000198221.8","-","-","-","ENSG00000214783.9","-","-","-","ENSG00000269911.1")
crispr_data_complete$gene_id[setdiff(1:56,a)] <- manually_ens_id_crispr
crispr_data$gene_id[which(is.na(crispr_data$gene_id))] <- crispr_data_complete$gene_id
saveRDS(crispr_data,"/home/egonie/kike/databases/hits_info_Liu_science_2015_ensids.rds")

# In the future this list could be further expanded with more studies targeting lncRNAs functionality
