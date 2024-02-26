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



