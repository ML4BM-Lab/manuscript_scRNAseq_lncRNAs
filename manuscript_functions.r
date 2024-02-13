# Manuscript functions
qc_metrics <- function(pipeline, mitochondrial_ens_ids, ribosomal_ens_ids )
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
  print(DimPlot(a, reduction = reduction, label = TRUE, pt.size = 0.3, repel = T,label.size = 5, cols = c("#c8723e","#9071c7","#84a142","#cb547e","#4daf95")) +ggtitle(title))
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




