# Manuscript functions
qc_metrics <- function(pipeline, mitochondrial_ens_ids )
{
  sce <- pipeline  
  stats <- perCellQCMetrics(sce, subsets=list(Mito=which(rownames(sce)%in% mitochondrial_ens_ids)))
  sum_log=data.frame(log10(stats$sum))
  detected_log=data.frame(log10(stats$detected))
  colData(sce) <- c(colData(sce),stats,sum_log,detected_log)
  
  return(sce)
}

emptydrops_filt <- function(sce, lower_ED, EmptyDrops_FDR_thres)
{
  e.out <- DropletUtils::emptyDrops(counts(sce), lower = lower_ED, niters = 10000)
  is.cell <- e.out$FDR <= EmptyDrops_FDR_thres
  print(paste("cells_detected: ",sum(is.cell, na.rm=TRUE)))
  non_empty_drops <- which(e.out$FDR < EmptyDrops_FDR_thres)
  sce_filt <- sce[,non_empty_drops]

  return(sce_filt)
}

doublet_analysis <- function(data)
{
  data2 <- scDblFinder(data, clusters=TRUE)
  print(table(data2$scDblFinder.class))
  data$isDoublet <- data2$scDblFinder.class == "doublet"
  return(data)
}

Filtering <- function(sce, cells_mito_threshold, cells_max_threshold, cells_min_genes_detected_threshold)
{
    total_max_mito_content <- sce$subsets_Mito_percent > cells_mito_threshold
    total_max_expression <- sce$sum > cells_max_threshold
    total_min_genes_detected <- sce$detected < cells_min_genes_detected_threshold

    discard <-  total_max_expression | total_max_mito_content | total_min_genes_detected
    print(table(discard))
    print(DataFrame(total_max_expression=sum(total_max_expression),total_max_mito_content=sum(total_max_mito_content),total_min_genes_detected=sum(total_min_genes_detected), Total=sum(discard)))
    sce$discard <- discard
    #filter and normalize data
    sce_filt <- sce[,discard==F]

    sce_filt <- logNormCounts(sce_filt)

    return(sce_filt)
}

top_genes <- function(data,threshold_minumun_gene_counts,threshold_cells_detected)
{
    out_data <- data[((rowSums(counts(data)) > threshold_minumun_gene_counts) & (rowSums(counts(data) != 0) > threshold_cells_detected)),]
    return(out_data)
}

red_dim <- function(data)
{
  set.seed(100100100)
  data <- runPCA(data)
  data <- runTSNE(data , ntop = 1000)
  data <- runUMAP(data, ntop = 1000)
  return(data)
}

clustering <- function(data, k, hg38_ensembl_gtf)
{
  g <- buildSNNGraph(data, use.dimred = "PCA", k = k ) # k is the number of nearest neighbors used to construct the graph. Higher k
  clust <- igraph::cluster_louvain(g)$membership
  print(table(clust))
  data$louvain_clusters <- factor(clust)

  # uniquifyFeatures
  gene_name <- hg38_ensembl_gtf$gene_name[match(rownames(data),hg38_ensembl_gtf$gene_id)]
  rowData(data) <- cbind(ens_id = rownames(data),gene_name)
  rownames(data) <- uniquifyFeatureNames(rownames(data), rowData(data)$gene_name)

  return(data)
}

red_dim_plots <- function(data)
{
  b=as.Seurat(data)
  b[['ident']] = b$louvain_clusters
  Idents(b) =  b$louvain_clusters
  print(DimPlot(b, reduction = "UMAP", label = TRUE, pt.size = 0.5) + NoLegend())
  print(plotUMAP(data, colour_by = "louvain_clusters", point_size = 0.5)  + guides(colour = guide_legend(override.aes = list(size=2))))
  print(DimPlot(b, reduction = "TSNE", label = TRUE, pt.size = 0.5) + NoLegend())
  print(plotTSNE(data, colour_by = "louvain_clusters", point_size = 0.5)  + guides(colour = guide_legend(override.aes = list(size=2))))
}

plot_markers <- function(sce, markers, title, group)
{
    print(DotPlot(as.Seurat(sce), assay = "RNA", features = intersect(markers,rownames(sce)), scale.by = "size",col.min =  0, cols = c("lightgrey", "darkred"), group.by = c(group)) + ggtitle(title)  +geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+ theme(legend.position="bottom", legend.direction="horizontal")) 
}

marker_umaps <- function(sce,clusters, annotation, subclustering = F, cells_subcluster = NA, new_label = NA, title, reduction = "UMAP",final_labels_order )
{
  #normally clusters=sce$louvain_clusters
  map <- data.frame(find=names(table(clusters)),replace=annotation)
  a=as.Seurat(sce)
  Idents(a) =  as.character(map[match(clusters, map$find), "replace"])
  sce$cell_types = Idents(a)
  sce$cell_types <- factor(sce$cell_types, levels = final_labels_order)
  Idents(a) = sce$cell_types
  if(subclustering == T)
  {
    sce$cell_types <- as.character(sce$cell_types)
    positions <- which(colnames(sce) %in% cells_subcluster)
    sce$cell_types[positions] = rep(new_label,length(cells_subcluster))
    sce$cell_types <- as.factor(sce$cell_types)
    Idents(a) <- sce$cell_types
  }
  #data_percentage <- (100*table(sce$cell_types)) / sum(table(sce$cell_types))
  #f <- melt(data_percentage)
  #print(ggplot(f,aes(Var.1,value))+geom_bar(stat="identity", position = "stack") +ggtitle(title))
  print(DimPlot(a, reduction = reduction, label = TRUE, pt.size = 0.3, repel = T,label.size = 5) +ggtitle(title))
  return(Idents(a))
}

kallisto_processing <- function(kallisto_dir,name, mito_gene)
{
    kallisto_data <- BUSpaRse::read_count_output(kallisto_dir, name = name)
    kallisto <- CreateSeuratObject(kallisto_data, project = "kallisto")
    kallisto <- as.SingleCellExperiment(kallisto)
    kallisto_sce <- qc_metrics(pipeline = kallisto, mitochondrial_ens_ids = mito_gene, ribosomal_ens_ids = ribosomal_ens_ids)

    return(kallisto_sce)
}

test_match_order <- function(x,y) {

if (all(x==y)) print('Perfect match in same order')

if (!all(x==y) && all(sort(x)==sort(y))) print('Perfect match in wrong order')

if (!all(x==y) && !all(sort(x)==sort(y))) print('No match')
}

cellRanger_processing <- function(cell_ranger_dir,mito_gene, multiplexing = F)
{
    cellRanger_data <- Read10X(data.dir = cell_ranger_dir, gene.column = 1)
    if (multiplexing == T)
    {
        cellRanger <- CreateSeuratObject(cellRanger_data[[1]], project = "cellRanger")
    }
    else {
       cellRanger <- CreateSeuratObject(cellRanger_data, project = "cellRanger")
    }
    cellRanger <- as.SingleCellExperiment(cellRanger)
    cellRanger_sce <- qc_metrics(pipeline = cellRanger,mitochondrial_ens_ids = mito_gene, ribosomal_ens_ids = ribosomal_ens_ids)

    return(cellRanger_sce)
}

# Function to generate the UpSet plots in the extended benchmarking
comp_pct <- function(count, PANEL, cut) {
  data.frame(count, PANEL, cut) %>% 
    group_by(PANEL, cut) %>% 
    mutate(pct = 100*count / sum(count)) %>% 
    pull(pct)
}


wnn_clustering <- function(seurat_object)
{
    seurat_object <- RunTFIDF(seurat_object)
    seurat_object <- FindTopFeatures(seurat_object, min.cutoff = 'q0')
    seurat_object <- RunSVD(seurat_object)
    seurat_object <- RunUMAP(seurat_object, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
    seurat_object <- FindMultiModalNeighbors(seurat_object, reduction.list = list("PCA", "lsi"), dims.list = list(1:50, 2:50))
    seurat_object <- RunUMAP(seurat_object, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
    seurat_object <- FindClusters(seurat_object, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
    DefaultAssay(seurat_object) <- "RNA"
    seurat_object_sce <- as.SingleCellExperiment(seurat_object)
    return(list(seurat_object,seurat_object_sce))
}

boxplot_nuclei_scMULTIOME <- function(input_RDS)
{
    common_expressed_genes_threshold <- readRDS(input_RDS)
    cellRanger_ATAC_expressed_genes_threshold <- common_expressed_genes_threshold[["cellRanger_ATAC"]]
    kallisto_ATAC_expressed_genes_threshold <- common_expressed_genes_threshold[["kallisto_nuclei_ATAC"]]
    cellRanger_ATAC_expressed_genes_threshold_lncRNAs <- common_expressed_genes_threshold[["cellRanger_nuclei_ATAC_lncRNAs"]]
    kallisto_ATAC_expressed_genes_threshold_lncRNAs <- common_expressed_genes_threshold[["kallisto_nuclei_ATAC_lncRNAs"]]
    cellRanger_ATAC_expressed_genes_threshold_PCs_genes <- common_expressed_genes_threshold[["cellRanger_nuclei_ATAC_PCs_genes"]]
    kallisto_ATAC_expressed_genes_threshold_PCs_genes <- common_expressed_genes_threshold[["kallisto_nuclei_ATAC_PCs_genes"]]

    df_threshold <- data.frame(pipeline = c(rep("Kallisto",length(kallisto_ATAC_expressed_genes_threshold)),rep("Cell Ranger",length(cellRanger_ATAC_expressed_genes_threshold)) ),common_expression = c(kallisto_ATAC_expressed_genes_threshold,cellRanger_ATAC_expressed_genes_threshold),common_expression_lncRNAs = c(kallisto_ATAC_expressed_genes_threshold_lncRNAs,cellRanger_ATAC_expressed_genes_threshold_lncRNAs), common_expression_PC_genes = c(kallisto_ATAC_expressed_genes_threshold_PCs_genes, cellRanger_ATAC_expressed_genes_threshold_PCs_genes))
    aux <- gsub("validations_paper/common_expressed_genes_threshold_","t",input_RDS)
    df_threshold$threshold <- gsub(".rds","",aux)

    p <- return(df_threshold)
}

odds_ratio_genes <- function(input_RDS, gene.activities)
{
  common_expressed_genes_threshold <- readRDS(input_RDS)
  cellRanger_ATAC_expressed_genes_threshold <- common_expressed_genes_threshold[["cellRanger_ATAC"]]
  kallisto_ATAC_expressed_genes_threshold <- common_expressed_genes_threshold[["kallisto_nuclei_ATAC"]]
  cellRanger_ATAC_expressed_genes_threshold_lncRNAs <- common_expressed_genes_threshold[["cellRanger_nuclei_ATAC_lncRNAs"]]
  kallisto_ATAC_expressed_genes_threshold_lncRNAs <- common_expressed_genes_threshold[["kallisto_nuclei_ATAC_lncRNAs"]]
  cellRanger_ATAC_expressed_genes_threshold_PCs_genes <- common_expressed_genes_threshold[["cellRanger_nuclei_ATAC_PCs_genes"]]
  kallisto_ATAC_expressed_genes_threshold_PCs_genes <- common_expressed_genes_threshold[["kallisto_nuclei_ATAC_PCs_genes"]]

  OR_all <- c()
  OR_lncRNAs <- c()
  OR_PCs_genes <- c()
   
  for(i in 1:length(cellRanger_ATAC_expressed_genes_threshold))
    {
      print(i)

      # with the if statements I am ignoring huge True Negatives
      if(cellRanger_ATAC_expressed_genes_threshold[i] != 0 | kallisto_ATAC_expressed_genes_threshold[i] != 0 )
      {
        m <- matrix(c(kallisto_ATAC_expressed_genes_threshold[i] + 0.5,(nrow(gene.activities) +0.5) - kallisto_ATAC_expressed_genes_threshold[i],cellRanger_ATAC_expressed_genes_threshold[i] + 0.5,(nrow(gene.activities) +0.5) - cellRanger_ATAC_expressed_genes_threshold[i] ), ncol = 2)
        OR_all <- c(OR_all,log2((m[1,1]/m[2,1])/(m[1,2]/m[2,2])))
      }
      if(cellRanger_ATAC_expressed_genes_threshold_lncRNAs[i] != 0 | kallisto_ATAC_expressed_genes_threshold_lncRNAs[i] != 0 )
      {
        m <- matrix(c(kallisto_ATAC_expressed_genes_threshold_lncRNAs[i] + 0.5,(nrow(gene.activities) +0.5) - kallisto_ATAC_expressed_genes_threshold_lncRNAs[i],cellRanger_ATAC_expressed_genes_threshold_lncRNAs[i] + 0.5,(nrow(gene.activities) +0.5) - cellRanger_ATAC_expressed_genes_threshold_lncRNAs[i] ), ncol = 2)
        OR_lncRNAs <- c(OR_lncRNAs,log2((m[1,1]/m[2,1])/(m[1,2]/m[2,2])))
      }
      if(cellRanger_ATAC_expressed_genes_threshold_PCs_genes[i] != 0 | kallisto_ATAC_expressed_genes_threshold_PCs_genes[i] != 0 )
      {
        m <- matrix(c(kallisto_ATAC_expressed_genes_threshold_PCs_genes[i] + 0.5,(nrow(gene.activities) +0.5) - kallisto_ATAC_expressed_genes_threshold_PCs_genes[i],cellRanger_ATAC_expressed_genes_threshold_PCs_genes[i] + 0.5,(nrow(gene.activities) +0.5) - cellRanger_ATAC_expressed_genes_threshold_PCs_genes[i] ), ncol = 2)
        OR_PCs_genes <- c(OR_PCs_genes,log2((m[1,1]/m[2,1])/(m[1,2]/m[2,2])))
      }
    }
    odds_results <- list(OR_all, OR_lncRNAs , OR_PCs_genes)
    return(odds_results)
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


