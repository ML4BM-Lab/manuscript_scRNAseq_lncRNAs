# Functions exclusive of the characterization
load_SEEKR_communities <- function(seekr_communities_path, gtf) 
{
  communities <- as.data.frame(read.csv(seekr_communities_path,header=T,row.names=1,sep=","))
  communities$row_names <- rownames(communities)
  a=rownames(communities)
  res <- str_match(a, "\\|\\s*(.*?)\\s*\\|")
  communities <- communities[!duplicated(res[,2],fromLast=T),]
  rownames(communities) <- res[,2][!duplicated(res[,2],fromLast=T)]
  #add gene names
  gene_name <- gtf$gene_name[match(rownames(communities),gtf$gene_id)]
  communities$gene_name <- gene_name

  return(communities)
}

# Characterization functions
# Get exonic length for each transcript. 
exonic_length_transcript <- function(gtf)
{
    exons = gtf[gtf$type=="exon",]
    current_gene <- exons$gene_id[1]
    current_transcript <- exons$transcript_id[1]
    length_transcript = 0
    t_lengths_all <- c()
    transcripts_all <- c()
    genes_all <- c()
    genes_names_all <- c()

    for (i in 1:nrow(exons))
    {
        print(i)
        gene_id <- exons$gene_id[i]
        transcript_id <- exons$transcript_id[i]

        if(current_transcript == transcript_id)
        {
            length_transcript <- exons$width[i] + length_transcript
        }
        else
        {
            t_lengths_all <- c(t_lengths_all, length_transcript)
            current_transcript <- transcript_id

            length_transcript <- exons$width[i]
            transcripts_all <- c(transcripts_all, exons$transcript_id[i-1])
            genes_all <- c(genes_all, exons$gene_id[i-1])
            genes_names_all <- c(genes_names_all, exons$gene_name[i-1])
        }
    }
    all_info <- as.data.frame(cbind(genes_all,genes_names_all,transcripts_all,t_lengths_all))
    #add last isoform manually
    last_isoform <- exons[nrow(exons),]
    all_info <- rbind(all_info,c(last_isoform$gene_id,last_isoform$gene_name,last_isoform$transcript_id,last_isoform$width))
    return(all_info)
}

# Get longest transcript for every gene and their exon coordinates
get_longest_transcript <- function(all_info_path)
{
    all_info <- readRDS(all_info_path)
    longest_transcripts <- all_info[1,]
    for (i in unique(all_info$genes_all))
    {
    print(nrow(longest_transcripts))
    transcripts <- all_info[all_info$genes_all == i,]
    longest_transcripts <- rbind(longest_transcripts,transcripts[which.max(transcripts$t_lengths_all),])
    }
    longest_transcripts_final <- longest_transcripts[-1,]
    return(longest_transcripts_final)
}

# get for every transcript the list of coordinates of its exons
transcripts_exons_coordinates <- function(gtf_path, longest_transcripts)
{
    gtf = rtracklayer::import(gtf_path)
    exons = gtf[gtf$type=="exon"]
    exons_longest_transcripts <- exons[which(exons$transcript_id %in% longest_transcripts$transcripts_all),]
    aux <- as.data.frame(exons_longest_transcripts)
    total_exon_length <- abs(aux$start-1-aux$end)
    exons_longest_transcripts$total_exon_length <- total_exon_length
    return(exons_longest_transcripts)
}

# Get the number of exons for the longest isoform of each gene
number_of_exons_longest_isoform <- function(gtf, longest_transcripts)
{
    exons = gtf[gtf$type=="exon",]
    number_exons <- c()
    gene_id <- c()
    transcript_id <- c()
    g_n <- exons[,12]
    g_i <- exons[,10]
    gene_name <- c()
    for (i in unique(exons$transcript_id))
    { 
        if (i %in% longest_transcripts$transcripts_all)
        {
            number_exons <- c(number_exons,nrow(exons[exons$transcript_id==i,]))
            transcript_id <- c(transcript_id,i)
            gene_id <- c(gene_id,unique(g_i[exons$transcript_id==i]))
            gene_name <- c(gene_name,unique(g_n[exons$transcript_id==i]))
            print(length(number_exons))
        }
    }
    n_exons_gene <- as.data.frame(cbind(gene_id,gene_name,transcript_id,number_exons))
    colnames(n_exons_gene) <- c("gene_id","gene_exons","transcript_id","number_exons")

    return(n_exons_gene)
}

number_proximal_exons <- function(gtf, longest_transcripts, distance_3UTR)
{
    exons = gtf[gtf$type=="exon",]
    number_exons <- c()
    gene_id <- c()
    transcript_id <- c()
    g_n <- exons[,12]
    g_i <- exons[,10]
    gene_name <- c()
    for (i in unique(exons$transcript_id))
    { 
        if (i %in% longest_transcripts$transcripts_all)
        {
            subset_exons <- exons[exons$transcript_id==i,]
            if(subset_exons$strand[1] == "+")
            {
                position_3UTR <- max(subset_exons$end)
                number_proximal_exons <- sum(subset_exons$end > (position_3UTR -  distance_3UTR))
            }
            else
            {
                position_3UTR <- min(subset_exons$start)
                number_proximal_exons <- sum(subset_exons$start < (position_3UTR+ distance_3UTR))
            }
            number_exons <- c(number_exons,number_proximal_exons)
            transcript_id <- c(transcript_id,i)
            gene_id <- c(gene_id,unique(g_i[exons$transcript_id==i]))
            gene_name <- c(gene_name,unique(g_n[exons$transcript_id==i]))
            print(length(number_exons))
        }
    }
    n_exons_gene <- as.data.frame(cbind(gene_id,gene_name,transcript_id,number_exons))
    colnames(n_exons_gene) <- c("gene_id","gene_exons","transcript_id","number_proximal_exons")

    return(n_exons_gene)
}


length_distributions <- function(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_sce_filt_clus, cellRanger_sce_filt_clus, STARsolo_sce_filt_clus, alevin_sce_filt_clus, lncrna_names, protein_coding_names,longest_transcripts,gene_name="gene_name")
{
  final_df <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(final_df) <- c("features","t_lengths_all","biotype","Threshold")
  for (j in 1:length(threshold_minumun_gene_counts_v))
  {
    threshold_minumun_gene_counts <- threshold_minumun_gene_counts_v[j]
    threshold_cells_detected <- threshold_cells_detected_v[j]
    print(paste("Threshold of minimun expression:",threshold_minumun_gene_counts))

    kallisto_top_genes <- top_genes(kallisto_sce_filt_clus,threshold_minumun_gene_counts,threshold_cells_detected )
    cellRanger_top_genes <- top_genes(cellRanger_sce_filt_clus,threshold_minumun_gene_counts,threshold_cells_detected )
    STARsolo_top_genes <- top_genes(STARsolo_sce_filt_clus,threshold_minumun_gene_counts,threshold_cells_detected )
    alevin_top_genes <- top_genes(alevin_sce_filt_clus,threshold_minumun_gene_counts,threshold_cells_detected )

    # Uniquely vs common
    input_list <- list(CellRanger = rownames(cellRanger_top_genes), STARsolo = rownames(STARsolo_top_genes), Kallisto = rownames(kallisto_top_genes), Salmon = rownames(alevin_top_genes))
    candidates_kallisto <- setdiff(input_list[[3]], c(input_list[[1]],input_list[[2]],input_list[[4]]))
    common_genes <- unique(intersect(input_list[[3]],intersect(input_list[[1]],intersect(input_list[[2]],input_list[[4]]))))
    if(gene_name=="gene_name")
    {
      longest_transcripts_lncRNAs <- longest_transcripts[longest_transcripts$genes_names_all %in% lncrna_names,]
      longest_transcripts_ex_lncRNAs <- longest_transcripts_lncRNAs[longest_transcripts_lncRNAs$genes_names_all %in% candidates_kallisto,]
      longest_transcripts_common_lncRNAs <- longest_transcripts_lncRNAs[longest_transcripts_lncRNAs$genes_names_all %in% common_genes,]

      longest_transcripts_PCs <- longest_transcripts[longest_transcripts$genes_names_all %in% protein_coding_names,]
      longest_transcripts_ex_PCs <- longest_transcripts_PCs[longest_transcripts_PCs$genes_names_all %in% candidates_kallisto,]
      longest_transcripts_common_PCs <- longest_transcripts_PCs[longest_transcripts_PCs$genes_names_all %in% common_genes,]
      
    }
    else
    {
      longest_transcripts_lncRNAs <- longest_transcripts[longest_transcripts$genes_all %in% lncrna_names,]
      longest_transcripts_ex_lncRNAs <- longest_transcripts_lncRNAs[longest_transcripts_lncRNAs$genes_all %in% candidates_kallisto,]
      longest_transcripts_common_lncRNAs <- longest_transcripts_lncRNAs[longest_transcripts_lncRNAs$genes_all %in% common_genes,]

      longest_transcripts_PCs <- longest_transcripts[longest_transcripts$genes_all %in% protein_coding_names,]
      longest_transcripts_ex_PCs <- longest_transcripts_PCs[longest_transcripts_PCs$genes_all %in% candidates_kallisto,]
      longest_transcripts_common_PCs <- longest_transcripts_PCs[longest_transcripts_PCs$genes_all %in% common_genes,]
    }
    all_lncRNAs <- data.frame(features = c(rep("exclusive_kallisto",nrow(longest_transcripts_ex_lncRNAs)),rep("common",nrow(longest_transcripts_common_lncRNAs))), t_lengths_all=c(as.numeric(longest_transcripts_ex_lncRNAs$t_lengths_all),as.numeric(longest_transcripts_common_lncRNAs$t_lengths_all)))
    all_lncRNAs$biotype = "LncRNA"
    all_PCs <- data.frame(features = c(rep("exclusive_kallisto",nrow(longest_transcripts_ex_PCs)),rep("common",nrow(longest_transcripts_common_PCs))), t_lengths_all=c(as.numeric(longest_transcripts_ex_PCs$t_lengths_all),as.numeric(longest_transcripts_common_PCs$t_lengths_all)))
    all_PCs$biotype = "Protein-coding"
    my_df <- rbind(all_lncRNAs,all_PCs)
    my_df$Threshold = threshold_minumun_gene_counts_v[j]
    final_df <- rbind(final_df,my_df)
  }
  final_df$Threshold <- as.factor(final_df$Threshold)
  return(final_df)
}

number_exons_distributions <- function(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_sce_filt_clus, cellRanger_sce_filt_clus, STARsolo_sce_filt_clus, alevin_sce_filt_clus, lncrna_names, protein_coding_names,n_exons_all,gene_name="gene_name")
{
  final_df <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(final_df) <- c("features","number_exons","biotype","Threshold")
  for (j in 1:length(threshold_minumun_gene_counts_v))
  {
    threshold_minumun_gene_counts <- threshold_minumun_gene_counts_v[j]
    threshold_cells_detected <- threshold_cells_detected_v[j]
    print(paste("Threshold of minimun expression:",threshold_minumun_gene_counts))

    kallisto_top_genes <- top_genes(kallisto_sce_filt_clus,threshold_minumun_gene_counts,threshold_cells_detected )
    cellRanger_top_genes <- top_genes(cellRanger_sce_filt_clus,threshold_minumun_gene_counts,threshold_cells_detected )
    STARsolo_top_genes <- top_genes(STARsolo_sce_filt_clus,threshold_minumun_gene_counts,threshold_cells_detected )
    alevin_top_genes <- top_genes(alevin_sce_filt_clus,threshold_minumun_gene_counts,threshold_cells_detected )

    # Uniquely vs common
    input_list <- list(CellRanger = rownames(cellRanger_top_genes), STARsolo = rownames(STARsolo_top_genes), Kallisto = rownames(kallisto_top_genes), Salmon = rownames(alevin_top_genes))
    candidates_kallisto <- setdiff(input_list[[3]], c(input_list[[1]],input_list[[2]],input_list[[4]]))
    common_genes <- unique(intersect(input_list[[3]],intersect(input_list[[1]],intersect(input_list[[2]],input_list[[4]]))))

    if (gene_name == "gene_name")
    {
      exons_lncRNAs <- n_exons_all[n_exons_all$gene_exons %in% lncrna_names,]
      exons_lncRNAs_ex <- exons_lncRNAs[exons_lncRNAs$gene_exons %in% candidates_kallisto,]
      exons_lncRNAs_common <- exons_lncRNAs[exons_lncRNAs$gene_exons %in% common_genes,]
      
      exons_PCs <- n_exons_all[n_exons_all$gene_exons %in% protein_coding_names,]
      exons_PC_ex <- exons_PCs[exons_PCs$gene_exons %in% candidates_kallisto,]
      exons_PCs_common <- exons_PCs[exons_PCs$gene_exons %in% common_genes,]
    }
    else
    {
      exons_lncRNAs <- n_exons_all[n_exons_all$gene_id %in% lncrna_names,]
      exons_lncRNAs_ex <- exons_lncRNAs[exons_lncRNAs$gene_id %in% candidates_kallisto,]
      exons_lncRNAs_common <- exons_lncRNAs[exons_lncRNAs$gene_id %in% common_genes,]
      
      exons_PCs <- n_exons_all[n_exons_all$gene_id %in% protein_coding_names,]
      exons_PC_ex <- exons_PCs[exons_PCs$gene_id %in% candidates_kallisto,]
      exons_PCs_common <- exons_PCs[exons_PCs$gene_id %in% common_genes,]
    }

    all_lncRNAs <- data.frame(features = c(rep("exclusive_kallisto",nrow(exons_lncRNAs_ex)),rep("common",nrow(exons_lncRNAs_common))), number_exons=c(as.numeric(exons_lncRNAs_ex$number_exons),as.numeric(exons_lncRNAs_common$number_exons)))
    all_lncRNAs$biotype = "LncRNA"
    all_PCs <- data.frame(features = c(rep("exclusive_kallisto",nrow(exons_PC_ex)),rep("common",nrow(exons_PCs_common))), number_exons=c(as.numeric(exons_PC_ex$number_exons),as.numeric(exons_PCs_common$number_exons)))
    all_PCs$biotype = "Protein-coding"
    my_df <- rbind(all_lncRNAs,all_PCs)
    my_df$Threshold = threshold_minumun_gene_counts_v[j]
    final_df <- rbind(final_df,my_df)
  }
  final_df$Threshold <- as.factor(final_df$Threshold)
  return(final_df)
}

percentage_of_repeats <- function(repeats_GR, GR_object)
{
    repeats_common_genes <- GenomicRanges::intersect(repeats_GR,GR_object,ignore.strand = F)
    a=as.data.frame(repeats_common_genes)
    b=as.data.frame(GR_object)
    return(100*sum(a$width)/sum(b$width))
}

repeats_results <- function(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_sce_filt_clus, cellRanger_sce_filt_clus, STARsolo_sce_filt_clus, alevin_sce_filt_clus, lncrna_names, protein_coding_names,gene_name="gene_name",hg38_repeats,exons_longest_transcript)
{
  final_repeats <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(final_repeats) <- c("ALL_genes","LncRNAs","Protein-coding")
  for (j in 1:length(threshold_minumun_gene_counts_v))
  {
    threshold_minumun_gene_counts <- threshold_minumun_gene_counts_v[j]
    threshold_cells_detected <- threshold_cells_detected_v[j]
    print(paste("Threshold of minimun expression:",threshold_minumun_gene_counts))

    kallisto_top_genes <- top_genes(kallisto_sce_filt_clus,threshold_minumun_gene_counts,threshold_cells_detected )
    cellRanger_top_genes <- top_genes(cellRanger_sce_filt_clus,threshold_minumun_gene_counts,threshold_cells_detected )
    STARsolo_top_genes <- top_genes(STARsolo_sce_filt_clus,threshold_minumun_gene_counts,threshold_cells_detected )
    alevin_top_genes <- top_genes(alevin_sce_filt_clus,threshold_minumun_gene_counts,threshold_cells_detected )

    # Uniquely vs common
    input_list <- list(CellRanger = rownames(cellRanger_top_genes), STARsolo = rownames(STARsolo_top_genes), Kallisto = rownames(kallisto_top_genes), Salmon = rownames(alevin_top_genes))
    candidates_kallisto <- setdiff(input_list[[3]], c(input_list[[1]],input_list[[2]],input_list[[4]]))
    common_genes <- unique(intersect(input_list[[3]],intersect(input_list[[1]],intersect(input_list[[2]],input_list[[4]]))))

    if(gene_name=="gene_name")
    {
      exons_common_genes <- exons_longest_transcript[which(exons_longest_transcript$gene_name %in%common_genes)]
      exons_common_genes_lncRNAs <- exons_common_genes[which(exons_common_genes$gene_name %in% lncrna_names)]
      exons_common_genes_pc <- exons_common_genes[which(exons_common_genes$gene_name %in% protein_coding_names)]

      exons_candidates_kallisto <- exons_longest_transcript[which(exons_longest_transcript$gene_name %in%candidates_kallisto)]
      exons_candidates_kallisto_lncRNAs <- exons_candidates_kallisto[which(exons_candidates_kallisto$gene_name %in% lncrna_names)]
      exons_candidates_kallisto_pc <- exons_candidates_kallisto[which(exons_candidates_kallisto$gene_name %in% protein_coding_names)]
    }
    else
    {
      exons_longest_transcript$gene_id <- gsub("\\..*","",exons_longest_transcript$gene_id)
      protein_coding_names <- gsub("\\..*","",protein_coding_names) 
      lncrna_names <- gsub("\\..*","",lncrna_names) 
      exons_common_genes <- exons_longest_transcript[which(exons_longest_transcript$gene_id %in% gsub("\\..*","",common_genes))]
      exons_common_genes_lncRNAs <- exons_common_genes[which(exons_common_genes$gene_id %in% gsub("\\..*","",lncrna_names))]
      exons_common_genes_pc <- exons_common_genes[which(exons_common_genes$gene_id %in% gsub("\\..*","",protein_coding_names))]

      exons_candidates_kallisto <- exons_longest_transcript[which(exons_longest_transcript$gene_id %in% gsub("\\..*","",candidates_kallisto))]
      exons_candidates_kallisto_lncRNAs <- exons_candidates_kallisto[which(exons_candidates_kallisto$gene_id %in% gsub("\\..*","",lncrna_names))]
      exons_candidates_kallisto_pc <- exons_candidates_kallisto[which(exons_candidates_kallisto$gene_id %in% gsub("\\..*","",protein_coding_names))]
    }

    repeat_content_common_genes <- c(percentage_of_repeats(hg38_repeats,exons_common_genes),percentage_of_repeats(hg38_repeats,exons_common_genes_lncRNAs),percentage_of_repeats(hg38_repeats,exons_common_genes_pc))

    repeat_content_candidates_kallisto <- c(percentage_of_repeats(hg38_repeats,exons_candidates_kallisto),percentage_of_repeats(hg38_repeats,exons_candidates_kallisto_lncRNAs),percentage_of_repeats(hg38_repeats,exons_candidates_kallisto_pc))

    repeats_out <- as.data.frame(rbind(repeat_content_common_genes,repeat_content_candidates_kallisto))
    rownames(repeats_out) <- c(paste(threshold_minumun_gene_counts,"Common",sep="_"),paste(threshold_minumun_gene_counts,"Exclusive",sep="_"))
    repeats_out$threshold = threshold_minumun_gene_counts
    final_repeats <- rbind(final_repeats,repeats_out)
  }
  colnames(final_repeats) <- c("ALL_genes","LncRNAs","Protein_coding","Threshold")
  final_repeats$type <- gsub(".*_","",rownames(final_repeats))
  final_repeats$Threshold <- as.factor(final_repeats$Threshold)

  return(final_repeats)
}

ratios_repeats <- function(final_repeats_percentage, dataset_name)
{
  final_repeats <- data.frame(matrix(ncol = 3, nrow = 0))
  for (i in levels(final_repeats_percentage$Threshold))
  {
    ratio_exclusive_common_lncRNAs <- final_repeats_percentage$LncRNAs[(final_repeats_percentage$Threshold==i) & (final_repeats_percentage$type=="Exclusive")]/final_repeats_percentage$LncRNAs[(final_repeats_percentage$Threshold==i) & (final_repeats_percentage$type=="Common")]
    print(ratio_exclusive_common_lncRNAs)

    ratio_exclusive_common_PCs <- final_repeats_percentage$Protein_coding[(final_repeats_percentage$Threshold==i) & (final_repeats_percentage$type=="Exclusive")]/final_repeats_percentage$Protein_coding[(final_repeats_percentage$Threshold==i) & (final_repeats_percentage$type=="Common")]
    print(ratio_exclusive_common_PCs)

    final_repeats <- rbind(final_repeats,c(ratio_exclusive_common_PCs,ratio_exclusive_common_lncRNAs,i))
  }
  colnames(final_repeats) <- c("Protein_coding","LncRNAs","Threshold")
  final_repeats$Protein_coding <- as.numeric(final_repeats$Protein_coding)
  final_repeats$LncRNAs <- as.numeric(final_repeats$LncRNAs)
  final_repeats$dataset <- dataset_name

  return(final_repeats)
}


SEEKR_communities <- function(communities, candidates_kallisto, common_genes,ens_id=F,title,lncrnas_ids)
{
  p_values_hypergeometric <- c()
  if(ens_id==F)
  {
    communities <- communities[which(communities$gene_name %in% intersect(c(candidates_kallisto,common_genes), lncrnas_ids)),]
    communities$is_kallisto_candidates <- communities$gene_name %in% candidates_kallisto
  }
  else
  {
    communities <- communities[which(rownames(communities) %in% intersect(c(candidates_kallisto,common_genes), lncrnas_ids)),]
    communities$is_kallisto_candidates <- rownames(communities) %in% candidates_kallisto
  }
  df1 <- as.data.frame(table(communities$Group, communities$is_kallisto_candidate))
  colnames(df1) <- c("Community", "Kallisto_candidates","Counts")
  df1$Percentage <- 100*df1$Counts/(df1$Counts[1:length(table(df1$Community))]+df1$Counts[(length(table(df1$Community))+1):nrow(df1)])
  print(df1)
  a1=c("ALL_comms",FALSE,sum(df1$Counts[df1$Kallisto_candidates==F]),(100*sum(df1$Counts[df1$Kallisto_candidates==F]))/(sum(df1$Counts[df1$Kallisto_candidates==F])+sum(df1$Counts[df1$Kallisto_candidates==T])))
  a2=c("ALL_comms",TRUE,sum(df1$Counts[df1$Kallisto_candidates==T]),(100*sum(df1$Counts[df1$Kallisto_candidates==T]))/(sum(df1$Counts[df1$Kallisto_candidates==F])+sum(df1$Counts[df1$Kallisto_candidates==T])))
  a3=as.data.frame(rbind(a1,a2))
  colnames(a3) = colnames(df1)
  df1 = rbind(df1,a3)
  df1$Counts = as.numeric(df1$Counts)
  df1$Percentage = as.numeric(df1$Percentage)
  #hypergeometric analysis. What is the probably that if I have in total 1923 exclusive lncRNAs out of 2550 lncRNAs, I have more than X lncRNAs in each category randomly
  for (i in 1:(length(table(df1$Community))-1))
  {
    p <- phyper(df1$Counts[i+length(table(df1$Community))-1]-1,sum(communities$is_kallisto_candidate),nrow(communities)-sum(communities$is_kallisto_candidate),df1$Counts[i]+df1$Counts[i+length(table(df1$Community))-1],lower.tail=FALSE)
    p_values_hypergeometric = c(p_values_hypergeometric,p)
  }

  return(list(df_all=df1,p_vals=p_values_hypergeometric))
}

SEEKR_results <- function(threshold_minumun_gene_counts, threshold_cells_detected_v, kallisto_sce_filt_clus, cellRanger_sce_filt_clus, STARsolo_sce_filt_clus, alevin_sce_filt_clus, seekr_communities,lncrnas_ids,ens_id=F)
{
  p_values_all  = c()
  all_df =  data.frame(matrix(nrow = 0, ncol = 5))
  colnames(all_df) <- c("Community","Kallisto_candidates","Counts","Percentage")
  for (j in 1:length(threshold_minumun_gene_counts_v))
  {
  threshold_minumun_gene_counts <- threshold_minumun_gene_counts_v[j]
  threshold_cells_detected <- threshold_cells_detected_v[j]
  print(paste("Threshold of minimun expression:",threshold_minumun_gene_counts))
  kallisto_top_genes <- top_genes(kallisto_sce_filt_clus,threshold_minumun_gene_counts,threshold_cells_detected )
  cellRanger_top_genes <- top_genes(cellRanger_sce_filt_clus,threshold_minumun_gene_counts,threshold_cells_detected )
  STARsolo_top_genes <- top_genes(STARsolo_sce_filt_clus,threshold_minumun_gene_counts,threshold_cells_detected )
  alevin_top_genes <- top_genes(alevin_sce_filt_clus,threshold_minumun_gene_counts,threshold_cells_detected )

  # Uniquely vs common
  input_list <- list(CellRanger = rownames(cellRanger_top_genes), STARsolo = rownames(STARsolo_top_genes), Kallisto = rownames(kallisto_top_genes), Salmon = rownames(alevin_top_genes))
  candidates_kallisto <- setdiff(input_list[[3]], c(input_list[[1]],input_list[[2]],input_list[[4]]))
  common_genes <- unique(intersect(input_list[[3]],intersect(input_list[[1]],intersect(input_list[[2]],input_list[[4]]))))

  SEEKR_results <- SEEKR_communities(seekr_communities,candidates_kallisto,common_genes,ens_id, title=paste("Min expression Threshold:",threshold_minumun_gene_counts,"UMI counts"),lncrnas_ids)
  a1 <- SEEKR_results[[2]]
  names(a1) <- c(paste(paste("Comm_",names(table(SEEKR_results[[1]][,1]))[-length(names(table(SEEKR_results[[1]][,1])))],sep=""),threshold_minumun_gene_counts,sep=":threshold "))
  p_values_all <- c(p_values_all,a1)
  SEEKR_results_df <- SEEKR_results[[1]]
  SEEKR_results_df$threshold = threshold_minumun_gene_counts
  all_df <- rbind(all_df,SEEKR_results_df)
  }
  return(list(all_df=all_df,p_vals=p_values_all))
}

evolution_plots_over_clusters <- function(kallisto_sce_filt_clus,k_clus,threshold_minumun_gene_counts,threshold_cells_detected,cellRanger_sce_filt_clus,STARsolo_sce_filt_clus,alevin_sce_filt_clus, lncrna_names, protein_coding_names,n_clusters)
{
  SI <- cell_type_specific_score_function(kallisto_sce_filt_clus,group_by="louvain_clusters", average_by="mean")
  cell_type_specific_score <- SI[["cell_type_specificity_score"]]
  counts_cell_specificity_index <- SI[["counts_cell_specificity_index"]]

  print(paste(k_clus," clusters: ",threshold_minumun_gene_counts," in ",threshold_cells_detected," cells",sep=""))   

  kallisto_top_genes <- top_genes(kallisto_sce_filt_clus,threshold_minumun_gene_counts,threshold_cells_detected )
  cellRanger_top_genes <- top_genes(cellRanger_sce_filt_clus,threshold_minumun_gene_counts,threshold_cells_detected )
  STARsolo_top_genes <- top_genes(STARsolo_sce_filt_clus,threshold_minumun_gene_counts,threshold_cells_detected )
  alevin_top_genes <- top_genes(alevin_sce_filt_clus,threshold_minumun_gene_counts,threshold_cells_detected )

  # Uniquely vs common
  input_list <- list(CellRanger = rownames(cellRanger_top_genes), STARsolo = rownames(STARsolo_top_genes), Kallisto = rownames(kallisto_top_genes), Salmon = rownames(alevin_top_genes))
  candidates_kallisto <- setdiff(input_list[[3]], c(input_list[[1]],input_list[[2]],input_list[[4]]))
  common_genes <- unique(intersect(input_list[[3]],intersect(input_list[[1]],intersect(input_list[[2]],input_list[[4]]))))

  specificity_index <- c(mean(cell_type_specific_score[common_genes]),mean(cell_type_specific_score[candidates_kallisto]),mean(cell_type_specific_score[intersect(common_genes,lncrna_names)]),mean(cell_type_specific_score[intersect(candidates_kallisto,lncrna_names)]),mean(cell_type_specific_score[intersect(common_genes,protein_coding_names)]),mean(cell_type_specific_score[intersect(candidates_kallisto,protein_coding_names)]))

  feature <- c(c(rep("All_genes",2)),c(rep("LncRNAs",2)),c(rep("Protein_coding",2)))
  type <- rep(c("common","kallisto_exclusive"),3)
  final_df <- as.data.frame(cbind(specificity_index,type,feature,threshold_minumun_gene_counts,k_clus ))
  final_df$specificity_index = as.numeric(final_df$specificity_index)

  #save all the values for creating a evolutionary violin plot
  specificity_index_vp<- c(cell_type_specific_score[common_genes],cell_type_specific_score[candidates_kallisto],cell_type_specific_score[intersect(common_genes,lncrna_names)],cell_type_specific_score[intersect(candidates_kallisto,lncrna_names)],cell_type_specific_score[intersect(common_genes,protein_coding_names)],cell_type_specific_score[intersect(candidates_kallisto,protein_coding_names)])

  feature_vp <- c(rep("All_genes",length(common_genes)+length(candidates_kallisto)),rep("LncRNAs",length(intersect(common_genes,lncrna_names))+length(intersect(candidates_kallisto,lncrna_names))),rep("Protein_coding",length(intersect(common_genes,protein_coding_names))+length(intersect(candidates_kallisto,protein_coding_names))))

  type_vp <- c(rep("common_genes",length(common_genes)),rep("kallisto_exclusive",length(candidates_kallisto)), rep("common_genes",length(intersect(common_genes,lncrna_names))), rep("kallisto_exclusive",length(intersect(candidates_kallisto,lncrna_names))),  rep("common_genes",length(intersect(common_genes,protein_coding_names))), rep("kallisto_exclusive",length(intersect(candidates_kallisto,protein_coding_names)))   )

  final_df_vp <- as.data.frame(cbind(specificity_index_vp,type_vp,feature_vp,threshold_minumun_gene_counts,k_clus ))
  final_df_vp$specificity_index = as.numeric(final_df_vp$specificity_index)
  final_df_vp$n_clusters <- n_clusters
  final_df_vp$total_counts <- c(rowSums(logcounts(kallisto_sce_filt_clus))[common_genes],rowSums(logcounts(kallisto_sce_filt_clus))[candidates_kallisto],rowSums(logcounts(kallisto_sce_filt_clus))[intersect(common_genes,lncrna_names)],rowSums(logcounts(kallisto_sce_filt_clus))[intersect(candidates_kallisto,lncrna_names)],rowSums(logcounts(kallisto_sce_filt_clus))[intersect(common_genes,protein_coding_names)],rowSums(logcounts(kallisto_sce_filt_clus))[intersect(candidates_kallisto,protein_coding_names)])

  final_df_list <- list("final_df"=final_df,"final_df_vp"=final_df_vp)
  return(final_df_list)
}


create_df_vp <- function(kallisto_sce_filt_clus,cellRanger_sce_filt_clus,STARsolo_sce_filt_clus,alevin_sce_filt_clus,lncrna_names,protein_coding_names)
{
  k_clusters <- c(25,10,5,3,2)
  df <- data.frame(matrix(ncol = 5, nrow = 0))
  df_vp <- data.frame(matrix(ncol = 5, nrow = 0))

  threshold_minumun_gene_counts_v <- c(250,100,50,25)
  threshold_cells_detected_v <- c(25,10,5,3)

  for (i in k_clusters)
  {
    kallisto_clustered <- clustering(kallisto_sce_filt_clus,i, change_rownames = F)
    n_clusters <- length(table(kallisto_clustered$louvain_clusters))

    for (j in 1:length(threshold_minumun_gene_counts_v))
    {   
        threshold_minumun_gene_counts <- threshold_minumun_gene_counts_v[j]
        threshold_cells_detected <- threshold_cells_detected_v[j]
        
        t=evolution_plots_over_clusters(kallisto_clustered,k_clus=i,threshold_minumun_gene_counts=threshold_minumun_gene_counts,threshold_cells_detected=threshold_cells_detected,cellRanger_sce_filt_clus,STARsolo_sce_filt_clus,alevin_sce_filt_clus, lncrna_names, protein_coding_names,n_clusters)
        t1 <- t[["final_df"]]
        df <- rbind(df,t1)
        t2 <- t[["final_df_vp"]]
        df_vp <- rbind(df_vp,t2)
        print(paste("number of neighbors:",i,"and number of clusters:",n_clusters))
    }
  }
  df_vp$k_clus <- factor(df_vp$k_clus, levels = rev(c("2","3","5","10","25")))
  df_vp$threshold_minumun_gene_counts <- factor(df_vp$threshold_minumun_gene_counts,levels = c("250","100","50","25"))

  return(df_vp)
}

main_expression_object <- function(dataset, dataset_name, threshold)
{
  a = dataset[(dataset$k_clus==25) & (dataset$threshold_minumun_gene_counts==threshold),]
  a$dataset <- dataset_name
  return(a)
}

violin_plot_expression <- function(my_df)
{
  colors <- c("#D4B996FF","#A07855FF")
  p2 <- ggplot(my_df[my_df$feature_vp=="LncRNAs",],aes(x = type_vp, y = total_counts,fill = type_vp)) + scale_y_continuous(trans='log10',breaks=c(10,100,1000, 10000),limits = range(my_df$total_counts))+   facet_wrap(~dataset, nrow = 1,strip.position="bottom") + geom_violin() +xlab("") + ylab("") + theme_classic() + theme(axis.text = element_text(size=13), axis.title.x = element_text(size=14))+ theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) + theme(strip.background = element_blank(),strip.placement = "outside",strip.text.x = element_text(angle = 45,size = 12)) + stat_compare_means(method = "wilcox.test",symnum.args = list(cutpoints = c(0, 0.0005, 0.005, 0.05, 0.1, 1), symbols = c("****", "***", "**", "*", "ns")),hjust = -1,aes(label = ..p.signif..)) +scale_fill_manual(values=colors) 

  p3 <- ggplot(my_df[my_df$feature_vp=="Protein_coding",],aes(x = type_vp, y = total_counts,fill = type_vp)) + scale_y_continuous(trans='log10',breaks=c(10,100,1000, 10000),limits = range(my_df$total_counts)) +   facet_wrap(~dataset, nrow = 1,strip.position="bottom") + geom_violin() +xlab("") + ylab("") + theme_classic() + theme(axis.text = element_text(size=13), axis.title.x = element_text(size=14)) + theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) +theme( strip.background = element_blank(),strip.placement = "outside",strip.text.x = element_text(angle = 45,size = 12) ) + stat_compare_means(method = "wilcox.test",symnum.args = list(cutpoints = c(0, 0.0005, 0.005, 0.05, 0.1, 1), symbols = c("****", "***", "**", "*", "ns")),hjust = -1,aes(label = ..p.signif..))  +scale_fill_manual(values=colors) 

  figure <- ggarrange(p2,p3,labels = c("Expression of exclusive lncRNAs vs common lncRNAs","Expression of exclusive PCs vs common PCs"),ncol = 2,nrow=1,legend = "bottom",common.legend = TRUE, hjust=c(-0.15,-0.25), vjust= c(0.15,0.15)) + theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) + theme( strip.background = element_blank(),strip.text = element_text(size = 15),strip.placement = "outside" )
  print(annotate_figure(figure,top = text_grob(""),left = text_grob("Expression (logcounts)", rot = 90, vjust = 2,hjust = -0.01,size = 13)))

  p2 <- ggplot(my_df[my_df$feature_vp=="LncRNAs",],aes(x = type_vp, y = total_counts,fill = type_vp)) + scale_y_continuous(trans='log10',breaks=c(10,100,1000, 10000),limits = range(my_df$total_counts))+   facet_wrap(~dataset, nrow = 1,strip.position="bottom") + geom_violin() +xlab("") + ylab("") + theme_classic() + theme(axis.text = element_text(size=13), axis.title.x = element_text(size=14))+ theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) +theme( strip.background = element_blank(),strip.placement = "outside",strip.text.x = element_text(angle = 45,size = 12) ) + stat_compare_means(aes(label = ifelse( p < 2.e-16, p.format, ifelse( p < 1.e-2, sprintf("p = %2.1e", as.numeric(..p.format..)), sprintf("p = %5.4f", as.numeric(..p.format..))))))  +scale_fill_manual(values=colors) 
  
  p3 <- ggplot(my_df[my_df$feature_vp=="Protein_coding",],aes(x = type_vp, y = total_counts,fill = type_vp)) + scale_y_continuous(trans='log10',breaks=c(10,100,1000, 10000),limits = range(my_df$total_counts)) + facet_wrap(~dataset, nrow = 1,strip.position="bottom") + geom_violin() +xlab("") + ylab("") + theme_classic() + theme(axis.text = element_text(size=13), axis.title.x = element_text(size=14)) + theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) + theme(strip.background = element_blank(),strip.placement = "outside",strip.text.x = element_text(angle = 45,size = 12)) + stat_compare_means(aes(label = ifelse( p < 2.e-16, p.format, ifelse( p < 1.e-2, sprintf("p = %2.1e", as.numeric(..p.format..)), sprintf("p = %5.4f", as.numeric(..p.format..))))))  +scale_fill_manual(values=colors) 
  
  figure <- ggarrange(p2,p3,labels = c("Expression of exclusive lncRNAs vs common lncRNAs","Expression of exclusive PCs vs common PCs"),ncol = 2,nrow=1,align = c("h"),legend = "bottom",common.legend = TRUE, hjust=c(-0.15,-0.25), vjust= c(0.15,0.15)) + theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) + theme( strip.background = element_blank(),strip.text = element_text(size = 15),strip.placement = "outside" )
  print(annotate_figure(figure,top = text_grob(""),left = text_grob("Expression (logcounts)", rot = 90, vjust = 2,hjust = -0.01,size = 13)))
}

main_object <- function(dataset, dataset_name, Threshold)
{
  a = dataset[(dataset$Threshold==Threshold),]
  a$dataset <- dataset_name
  return(a)
}

violin_plot_length <- function(my_df)
{
  colors <- c("#D4B996FF","#A07855FF")
  p2 <- ggplot(my_df[my_df$biotype=="LncRNA",],aes(x = features, y = t_lengths_all,fill = features)) + scale_y_continuous(trans='log10',breaks=c(10,100,1000, 10000,100000),limits = range(my_df$t_lengths_all))+   facet_wrap(~dataset, nrow = 1,strip.position="bottom") + geom_violin() + geom_boxplot(width=0.2,position=position_dodge(width = 0.9))+xlab("") + ylab("") + theme_classic() + theme(axis.text = element_text(size=13), axis.title.x = element_text(size=14))+ theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) + theme( strip.background = element_blank(),strip.placement = "outside",strip.text.x = element_text(angle = 45,size = 12) ) + stat_compare_means(method = "wilcox.test",symnum.args = list(cutpoints = c(0, 0.0005, 0.005, 0.05, 0.1, 1), symbols = c("****", "***", "**", "*", "ns")),hjust = -1,aes(label = ..p.signif..)) +scale_fill_manual(values=colors) 

  p3 <- ggplot(my_df[my_df$biotype=="Protein-coding",],aes(x = features, y = t_lengths_all,fill = features)) + scale_y_continuous(trans='log10',breaks=c(10,100,1000, 10000,100000),limits = range(my_df$t_lengths_all)) +   facet_wrap(~dataset, nrow = 1,strip.position="bottom") + geom_violin() + geom_boxplot(width=0.2,position=position_dodge(width = 0.9)) +xlab("") + ylab("") + theme_classic() + theme(axis.text = element_text(size=13), axis.title.x = element_text(size=14)) + theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) +theme( strip.background = element_blank(),strip.placement = "outside",strip.text.x = element_text(angle = 45,size = 12) ) + stat_compare_means(method = "wilcox.test",symnum.args = list(cutpoints = c(0, 0.0005, 0.005, 0.05, 0.1, 1), symbols = c("****", "***", "**", "*", "ns")),hjust = -1,aes(label = ..p.signif..))  +scale_fill_manual(values=colors) 

  figure <- ggarrange(p2,p3,labels = c("Length of lncRNAs","Length of Protein-coding genes"),ncol = 2,nrow=1,legend = "bottom",common.legend = TRUE, hjust=c(-0.15,-0.25), vjust= c(0.15,0.15)) + theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) + theme( strip.background = element_blank(),strip.text = element_text(size = 15),strip.placement = "outside" )
  print(annotate_figure(figure,top = text_grob(""),left = text_grob("Length", rot = 90, vjust = 2,hjust = -0.01,size = 13)))

  p2 <- ggplot(my_df[my_df$biotype=="LncRNA",],aes(x = features, y = t_lengths_all,fill = features)) + scale_y_continuous(trans='log10',breaks=c(10,100,1000, 10000,100000),limits = range(my_df$t_lengths_all))+   facet_wrap(~dataset, nrow = 1,strip.position="bottom") + geom_violin() + geom_boxplot(width=0.2,position=position_dodge(width = 0.9))+xlab("") + ylab("") + theme_classic() + theme(axis.text = element_text(size=13), axis.title.x = element_text(size=14))+ theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) + theme( strip.background = element_blank(),strip.placement = "outside",strip.text.x = element_text(angle = 45,size = 12) ) + stat_compare_means(aes(label = ifelse( p < 2.e-16, p.format, ifelse( p < 1.e-2, sprintf("p = %2.1e", as.numeric(..p.format..)), sprintf("p = %5.4f", as.numeric(..p.format..)))))) +scale_fill_manual(values=colors) 

  p3 <- ggplot(my_df[my_df$biotype=="Protein-coding",],aes(x = features, y = t_lengths_all,fill = features)) + scale_y_continuous(trans='log10',breaks=c(10,100,1000, 10000,100000),limits = range(my_df$t_lengths_all)) +   facet_wrap(~dataset, nrow = 1,strip.position="bottom") + geom_violin() + geom_boxplot(width=0.2,position=position_dodge(width = 0.9)) +xlab("") + ylab("") + theme_classic() + theme(axis.text = element_text(size=13), axis.title.x = element_text(size=14)) + theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) +theme( strip.background = element_blank(),strip.placement = "outside",strip.text.x = element_text(angle = 45,size = 12) ) + stat_compare_means(aes(label = ifelse( p < 2.e-16, p.format, ifelse( p < 1.e-2, sprintf("p = %2.1e", as.numeric(..p.format..)), sprintf("p = %5.4f", as.numeric(..p.format..))))))  +scale_fill_manual(values=colors) 

  figure <- ggarrange(p2,p3,labels = c("Length of lncRNAs","Length of Protein-coding genes"),ncol = 2,nrow=1,legend = "bottom",common.legend = TRUE, hjust=c(-0.15,-0.25), vjust= c(0.15,0.15)) + theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) + theme( strip.background = element_blank(),strip.text = element_text(size = 15),strip.placement = "outside" )
  print(annotate_figure(figure,top = text_grob(""),left = text_grob("Length", rot = 90, vjust = 2,hjust = -0.01,size = 13)))

}


violin_plot_number_exons <- function(my_df)
{
  colors <- c("#D4B996FF","#A07855FF")
  p2 <- ggplot(my_df[my_df$biotype=="LncRNA",],aes(x = features, y = number_exons,fill = features)) + scale_y_continuous(trans='log10',breaks=c(10,100,1000),limits = range(my_df$number_exons))+   facet_wrap(~dataset, nrow = 1,strip.position="bottom") + geom_violin() + geom_boxplot(width=0.2,position=position_dodge(width = 0.9))+xlab("") + ylab("") + theme_classic() + theme(axis.text = element_text(size=13), axis.title.x = element_text(size=14))+ theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) + theme( strip.background = element_blank(),strip.placement = "outside",strip.text.x = element_text(angle = 45,size = 12) ) + stat_compare_means(method = "wilcox.test",symnum.args = list(cutpoints = c(0, 0.0005, 0.005, 0.05, 0.1, 1), symbols = c("****", "***", "**", "*", "ns")),hjust = -1,aes(label = ..p.signif..)) +scale_fill_manual(values=colors) 

  p3 <- ggplot(my_df[my_df$biotype=="Protein-coding",],aes(x = features, y = number_exons,fill = features)) + scale_y_continuous(trans='log10',breaks=c(10,100,1000),limits = range(my_df$number_exons)) +   facet_wrap(~dataset, nrow = 1,strip.position="bottom") + geom_violin() + geom_boxplot(width=0.2,position=position_dodge(width = 0.9)) +xlab("") + ylab("") + theme_classic() + theme(axis.text = element_text(size=13), axis.title.x = element_text(size=14)) + theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) +theme( strip.background = element_blank(),strip.placement = "outside",strip.text.x = element_text(angle = 45,size = 12) ) + stat_compare_means(method = "wilcox.test",symnum.args = list(cutpoints = c(0, 0.0005, 0.005, 0.05, 0.1, 1), symbols = c("****", "***", "**", "*", "ns")),hjust = -1,aes(label = ..p.signif..))  +scale_fill_manual(values=colors) 

  figure <- ggarrange(p2,p3,labels = c("Number of exons of lncRNAs","Number of exons of Protein-coding genes"),ncol = 2,nrow=1,legend = "bottom",common.legend = TRUE, hjust=c(-0.15,-0.25), vjust= c(0.15,0.15)) + theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) + theme( strip.background = element_blank(),strip.text = element_text(size = 15),strip.placement = "outside" )
  print(annotate_figure(figure,top = text_grob(""),left = text_grob("Number of exons", rot = 90, vjust = 2,hjust = -0.01,size = 13)))

  p2 <- ggplot(my_df[my_df$biotype=="LncRNA",],aes(x = features, y = number_exons,fill = features)) + scale_y_continuous(trans='log10',breaks=c(10,100,1000),limits = range(my_df$number_exons))+   facet_wrap(~dataset, nrow = 1,strip.position="bottom") + geom_violin() + geom_boxplot(width=0.2,position=position_dodge(width = 0.9))+xlab("") + ylab("") + theme_classic() + theme(axis.text = element_text(size=13), axis.title.x = element_text(size=14))+ theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) + theme( strip.background = element_blank(),strip.placement = "outside",strip.text.x = element_text(angle = 45,size = 12) ) + stat_compare_means(aes(label = ifelse( p < 2.e-16, p.format, ifelse( p < 1.e-2, sprintf("p = %2.1e", as.numeric(..p.format..)), sprintf("p = %5.4f", as.numeric(..p.format..)))))) +scale_fill_manual(values=colors) 

  p3 <- ggplot(my_df[my_df$biotype=="Protein-coding",],aes(x = features, y = number_exons,fill = features)) + scale_y_continuous(trans='log10',breaks=c(10,100,1000),limits = range(my_df$number_exons)) +   facet_wrap(~dataset, nrow = 1,strip.position="bottom") + geom_violin() + geom_boxplot(width=0.2,position=position_dodge(width = 0.9)) +xlab("") + ylab("") + theme_classic() + theme(axis.text = element_text(size=13), axis.title.x = element_text(size=14)) + theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) +theme( strip.background = element_blank(),strip.placement = "outside",strip.text.x = element_text(angle = 45,size = 12) ) + stat_compare_means(aes(label = ifelse( p < 2.e-16, p.format, ifelse( p < 1.e-2, sprintf("p = %2.1e", as.numeric(..p.format..)), sprintf("p = %5.4f", as.numeric(..p.format..))))))  +scale_fill_manual(values=colors) 

  figure <- ggarrange(p2,p3,labels = c("Number of exons of lncRNAs","Number of exons of Protein-coding genes"),ncol = 2,nrow=1,legend = "bottom",common.legend = TRUE, hjust=c(-0.15,-0.25), vjust= c(0.15,0.15)) + theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) + theme( strip.background = element_blank(),strip.text = element_text(size = 15),strip.placement = "outside" )
  print(annotate_figure(figure,top = text_grob(""),left = text_grob("Number of exons", rot = 90, vjust = 2,hjust = -0.01,size = 13)))
}

violin_plot_number_proximal_exons <- function(my_df)
{
    # Hypothesis 1: Exclusive is higher than common
    colors <- c("#D4B996FF","#A07855FF")
    p2 <- ggplot(my_df[my_df$biotype=="LncRNA",],aes(x = features, y = number_exons,fill = features)) + scale_y_continuous(trans='log10',limits = range(my_df$number_exons))+   facet_wrap(~dataset, nrow = 1,strip.position="bottom") + geom_violin() + geom_boxplot(width=0.2,position=position_dodge(width = 0.9))+xlab("") + ylab("") + theme_classic() + theme(axis.text = element_text(size=13), axis.title.x = element_text(size=14))+ theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) + theme( strip.background = element_blank(),strip.placement = "outside",strip.text.x = element_text(angle = 45,size = 12) ) + stat_compare_means(method = "wilcox.test",method.args = list(alternative = "less"),symnum.args = list(cutpoints = c(0, 0.0005, 0.005, 0.05, 0.1, 1), symbols = c("****", "***", "**", "*", "ns")),hjust = -1,aes(label = ..p.signif..)) +scale_fill_manual(values=colors) 

    p3 <- ggplot(my_df[my_df$biotype=="Protein-coding",],aes(x = features, y = number_exons,fill = features)) + scale_y_continuous(trans='log10',limits = range(my_df$number_exons)) +   facet_wrap(~dataset, nrow = 1,strip.position="bottom") + geom_violin() + geom_boxplot(width=0.2,position=position_dodge(width = 0.9)) +xlab("") + ylab("") + theme_classic() + theme(axis.text = element_text(size=13), axis.title.x = element_text(size=14)) + theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) +theme( strip.background = element_blank(),strip.placement = "outside",strip.text.x = element_text(angle = 45,size = 12) ) + stat_compare_means(method = "wilcox.test",method.args = list(alternative = "less"),symnum.args = list(cutpoints = c(0, 0.0005, 0.005, 0.05, 0.1, 1), symbols = c("****", "***", "**", "*", "ns")),hjust = -1,aes(label = ..p.signif..))  +scale_fill_manual(values=colors) 

    figure <- ggarrange(p2,p3,labels = c("Number of proximal exons of lncRNAs: Hyp. Exclusive have more exons than common","Number of exons of proximal protein-coding genes: Hyp. Exclusive have more exons than common"),ncol = 2,nrow=1,legend = "bottom",common.legend = TRUE, hjust=c(-0.15,-0.25), vjust= c(0.15,0.15)) + theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) + theme( strip.background = element_blank(),strip.text = element_text(size = 15),strip.placement = "outside" )
    print(annotate_figure(figure,top = text_grob(""),left = text_grob("Number of exons", rot = 90, vjust = 2,hjust = -0.01,size = 15)))

    p2 <- ggplot(my_df[my_df$biotype=="LncRNA",],aes(x = features, y = number_exons,fill = features)) + scale_y_continuous(trans='log10',limits = range(my_df$number_exons))+   facet_wrap(~dataset, nrow = 1,strip.position="bottom") + geom_violin() + geom_boxplot(width=0.2,position=position_dodge(width = 0.9))+xlab("") + ylab("") + theme_classic() + theme(axis.text = element_text(size=13), axis.title.x = element_text(size=14))+ theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) + theme( strip.background = element_blank(),strip.placement = "outside",strip.text.x = element_text(angle = 45,size = 12) ) + stat_compare_means(method = "wilcox.test",method.args = list(alternative = "less"),aes(label = ifelse( p < 2.e-16, p.format, ifelse( p < 1.e-2, sprintf("p = %2.1e", as.numeric(..p.format..)), sprintf("p = %5.4f", as.numeric(..p.format..)))))) +scale_fill_manual(values=colors) 

    p3 <- ggplot(my_df[my_df$biotype=="Protein-coding",],aes(x = features, y = number_exons,fill = features)) + scale_y_continuous(trans='log10',limits = range(my_df$number_exons)) +   facet_wrap(~dataset, nrow = 1,strip.position="bottom") + geom_violin() + geom_boxplot(width=0.2,position=position_dodge(width = 0.9)) +xlab("") + ylab("") + theme_classic() + theme(axis.text = element_text(size=13), axis.title.x = element_text(size=14)) + theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) +theme( strip.background = element_blank(),strip.placement = "outside",strip.text.x = element_text(angle = 45,size = 12) ) + stat_compare_means(method = "wilcox.test",method.args = list(alternative = "less"),aes(label = ifelse( p < 2.e-16, p.format, ifelse( p < 1.e-2, sprintf("p = %2.1e", as.numeric(..p.format..)), sprintf("p = %5.4f", as.numeric(..p.format..))))))  +scale_fill_manual(values=colors) 

    figure <- ggarrange(p2,p3,labels = c("Number of proximal exons of lncRNAs: Hyp. Exclusive have more exons than common","Number of exons of proximal protein-coding genes: Hyp. Exclusive have more exons than common"),ncol = 2,nrow=1,legend = "bottom",common.legend = TRUE, hjust=c(-0.15,-0.25), vjust= c(0.15,0.15)) + theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) + theme( strip.background = element_blank(),strip.text = element_text(size = 15),strip.placement = "outside" )
    print(annotate_figure(figure,top = text_grob(""),left = text_grob("Number of exons", rot = 90, vjust = 2,hjust = -0.01,size = 15)))

    # Hypothesis 2: Exclusive is smaller than common
    colors <- c("#D4B996FF","#A07855FF")
    p2 <- ggplot(my_df[my_df$biotype=="LncRNA",],aes(x = features, y = number_exons,fill = features)) + scale_y_continuous(trans='log10',limits = range(my_df$number_exons))+   facet_wrap(~dataset, nrow = 1,strip.position="bottom") + geom_violin() + geom_boxplot(width=0.2,position=position_dodge(width = 0.9))+xlab("") + ylab("") + theme_classic() + theme(axis.text = element_text(size=13), axis.title.x = element_text(size=14))+ theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) + theme( strip.background = element_blank(),strip.placement = "outside",strip.text.x = element_text(angle = 45,size = 12) ) + stat_compare_means(method = "wilcox.test",method.args = list(alternative = "greater"),symnum.args = list(cutpoints = c(0, 0.0005, 0.005, 0.05, 0.1, 1), symbols = c("****", "***", "**", "*", "ns")),hjust = -1,aes(label = ..p.signif..)) +scale_fill_manual(values=colors) 

    p3 <- ggplot(my_df[my_df$biotype=="Protein-coding",],aes(x = features, y = number_exons,fill = features)) + scale_y_continuous(trans='log10',limits = range(my_df$number_exons)) +   facet_wrap(~dataset, nrow = 1,strip.position="bottom") + geom_violin() + geom_boxplot(width=0.2,position=position_dodge(width = 0.9)) +xlab("") + ylab("") + theme_classic() + theme(axis.text = element_text(size=13), axis.title.x = element_text(size=14)) + theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) +theme( strip.background = element_blank(),strip.placement = "outside",strip.text.x = element_text(angle = 45,size = 12) ) + stat_compare_means(method = "wilcox.test",method.args = list(alternative = "greater"),symnum.args = list(cutpoints = c(0, 0.0005, 0.005, 0.05, 0.1, 1), symbols = c("****", "***", "**", "*", "ns")),hjust = -1,aes(label = ..p.signif..))  +scale_fill_manual(values=colors) 

    figure <- ggarrange(p2,p3,labels = c("Number of proximal exons of lncRNAs: Hyp. Exclusive have less exons than common","Number of exons of proximal protein-coding genes: Hyp. Exclusive have less exons than common"),ncol = 2,nrow=1,legend = "bottom",common.legend = TRUE, hjust=c(-0.15,-0.25), vjust= c(0.15,0.15)) + theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) + theme( strip.background = element_blank(),strip.text = element_text(size = 15),strip.placement = "outside" )
    print(annotate_figure(figure,top = text_grob(""),left = text_grob("Number of exons", rot = 90, vjust = 2,hjust = -0.01,size = 15)))

    p2 <- ggplot(my_df[my_df$biotype=="LncRNA",],aes(x = features, y = number_exons,fill = features)) + scale_y_continuous(trans='log10',limits = range(my_df$number_exons))+   facet_wrap(~dataset, nrow = 1,strip.position="bottom") + geom_violin() + geom_boxplot(width=0.2,position=position_dodge(width = 0.9))+xlab("") + ylab("") + theme_classic() + theme(axis.text = element_text(size=13), axis.title.x = element_text(size=14))+ theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) + theme( strip.background = element_blank(),strip.placement = "outside",strip.text.x = element_text(angle = 45,size = 12) ) + stat_compare_means(method = "wilcox.test",method.args = list(alternative = "greater"),aes(label = ifelse( p < 2.e-16, p.format, ifelse( p < 1.e-2, sprintf("p = %2.1e", as.numeric(..p.format..)), sprintf("p = %5.4f", as.numeric(..p.format..)))))) +scale_fill_manual(values=colors) 

    p3 <- ggplot(my_df[my_df$biotype=="Protein-coding",],aes(x = features, y = number_exons,fill = features)) + scale_y_continuous(trans='log10',limits = range(my_df$number_exons)) +   facet_wrap(~dataset, nrow = 1,strip.position="bottom") + geom_violin() + geom_boxplot(width=0.2,position=position_dodge(width = 0.9)) +xlab("") + ylab("") + theme_classic() + theme(axis.text = element_text(size=13), axis.title.x = element_text(size=14)) + theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) +theme( strip.background = element_blank(),strip.placement = "outside",strip.text.x = element_text(angle = 45,size = 12) ) + stat_compare_means(method = "wilcox.test",method.args = list(alternative = "greater"),aes(label = ifelse( p < 2.e-16, p.format, ifelse( p < 1.e-2, sprintf("p = %2.1e", as.numeric(..p.format..)), sprintf("p = %5.4f", as.numeric(..p.format..))))))  +scale_fill_manual(values=colors) 

    figure <- ggarrange(p2,p3,labels = c("Number of proximal exons of lncRNAs: Hyp. Exclusive have less exons than common","Number of exons of proximal protein-coding genes: Hyp. Exclusive have less exons than common"),ncol = 2,nrow=1,legend = "bottom",common.legend = TRUE, hjust=c(-0.15,-0.25), vjust= c(0.15,0.15)) + theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm'),axis.ticks.x=element_blank(),axis.text.x=element_blank()) + theme( strip.background = element_blank(),strip.text = element_text(size = 15),strip.placement = "outside" )
    print(annotate_figure(figure,top = text_grob(""),left = text_grob("Number of exons", rot = 90, vjust = 2,hjust = -0.01,size = 15)))
}


main_seekr_object <- function(dataset, dataset_name, threshold)
{
    dataset <- dataset[[1]]
    a = dataset[(dataset$threshold==threshold),]
    a$dataset <- dataset_name
    return(a)
}

seekr_barplot <- function(my_df)
{
  levels(my_df$Kallisto_candidates)[levels(my_df$Kallisto_candidates)=="TRUE"] = "exclusive_kallisto"
  levels(my_df$Kallisto_candidates)[levels(my_df$Kallisto_candidates)=="FALSE"] = "common"
  levels(my_df$Community)[levels(my_df$Community)=="ALL_comms"] = "ALL"

  colors <- c("#D4B996FF","#A07855FF")
  p1 <- ggplot(my_df,aes(x=Community, y=Percentage, fill=Kallisto_candidates)) +geom_bar(stat="identity") +   facet_wrap(~dataset, nrow = 3,scales="free") +xlab("SEEKR communities") + ylab("Percentage") + theme_classic() + theme(axis.text = element_text(size=13), axis.title.x = element_text(size=14))+ theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm')) +scale_fill_discrete(labels = c("Common_genes", "Kallisto_exclusive")) +  theme(strip.background = element_blank(), strip.placement = "outside",strip.text.x = element_text(size = 14)) +scale_fill_manual(values=colors) +ylim(0,101)
  print(p1)
}


main_SI_object <- function(dataset, dataset_name, threshold)
{
    a = dataset[(dataset$threshold_minumun_gene_counts==threshold),]
    a$dataset <- dataset_name

    a$cluster_categories <- as.factor(a$n_clusters)
    levels(a$cluster_categories)[((levels(a$cluster_categories)=="4")) | (levels(a$cluster_categories)=="5") | ((levels(a$cluster_categories)=="6")) |  ((levels(a$cluster_categories)=="7") | (levels(a$cluster_categories)=="8") | (levels(a$cluster_categories)=="9"))] = "[5-9]"
    levels(a$cluster_categories)[(levels(a$cluster_categories)=="10") | ((levels(a$cluster_categories)=="11") | (levels(a$cluster_categories)=="12") | (levels(a$cluster_categories)=="13")) | ((levels(a$cluster_categories)=="14")) | ((levels(a$cluster_categories)=="15"))] = "[10-15]"
    levels(a$cluster_categories)[(levels(a$cluster_categories)=="17") | (levels(a$cluster_categories)=="18") | ((levels(a$cluster_categories)=="19") | (levels(a$cluster_categories)=="20") | (levels(a$cluster_categories)=="21")) | ((levels(a$cluster_categories)=="22"))] = "[16-22]"
    levels(a$cluster_categories)[(levels(a$cluster_categories)=="23") | ((levels(a$cluster_categories)=="24")) | ((levels(a$cluster_categories)=="25")) | (levels(a$cluster_categories)=="26") | ((levels(a$cluster_categories)=="27")) | ((levels(a$cluster_categories)=="28"))| ((levels(a$cluster_categories)=="30")) | ((levels(a$cluster_categories)=="31")) | ((levels(a$cluster_categories)=="33")) | ((levels(a$cluster_categories)=="34"))] = "[>23]"
    return(a)
}

violin_plot_SI_ob2 <- function(my_df, hypothesis = "two.sided")
{
  my_df$k_clus = factor(as.character(my_df$k_clus),levels = c("25", "10", "5","3","2"))

  colors <- c("#E3CD81FF","#B1B3B3FF")
  if (hypothesis=="less")
  {
    aux = "more"
  }
  else
  {
    aux = "less"
  }
  my_df <- my_df[(my_df$feature_vp == "Protein_coding") | ((my_df$feature_vp == "LncRNAs") & (my_df$type_vp == "kallisto_exclusive")) ,]
  p1 <- ggplot(my_df,aes(x = cluster_categories, y = specificity_index,fill = feature_vp)) +   facet_wrap(~dataset, nrow = 3,scales="free") + geom_violin() + xlab("Number of clusters") + ylab("SI") + theme_classic() + theme(axis.text = element_text(size=13), axis.title.x = element_text(size=14))+ theme(legend.text=element_text(size=12),legend.title=element_blank(),legend.key.size = unit(1.5, 'cm')) + theme( strip.background = element_blank(),strip.placement = "outside",strip.text.x = element_text(size = 14) ) + stat_compare_means(method = "wilcox.test",method.args = list(alternative = hypothesis),symnum.args = list(cutpoints = c(0, 0.0005, 0.005, 0.05, 0.1, 1), symbols = c("****", "***", "**", "*", "ns")),hjust = -0,aes(label = ..p.signif..)) +scale_fill_manual(values=colors) + theme(legend.position="bottom") +ggtitle(paste("Specificity index: Exclusive LncRNAs have", aux, "SI than ALL protein-coding genes")) + theme(plot.title = element_text(size=16,hjust = 0.5))

  print(p1)
} 

counts_crispr <- function(dataset_all, dataset_name)
{
  dataset <- dataset_all[[1]]
  a = table(dataset$Threshold, dataset$features)
  reshaped_a <- melt(a)
  colnames(reshaped_a) <- c("Threshold","features","crispr_genes")
  reshaped_a$dataset <- dataset_name

  round(dataset_all[[2]],4)
  reshaped_a$p_value_hypergeom <- as.numeric(c(signif(dataset_all[[2]][,3],3),signif(dataset_all[[2]][,4],3)))
  return(reshaped_a)
}
generate_supp_table <- function(dataset_all, dataset_name)
{
    dataset <- dataset_all[[1]]
    dataset <- dataset[dataset$feature=="exclusive_kallisto",]
    dataset <- dataset[,1:4]
    dataset$dataset <- dataset_name
    return(dataset)
}
