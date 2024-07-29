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

filtered_gene_expression <- function(sce, max_gene_expression, lncrna_ens_ids,  protein_coding_ens_ids)
{
  lncRNAs_expression <- (rowSums(logcounts(sce))[lncrna_ens_ids[lncrna_ens_ids %in% rownames(sce)==T]])
  PCGs_expression <- (rowSums(logcounts(sce))[protein_coding_ens_ids[protein_coding_ens_ids %in% rownames(sce)==T]])
  filtered_lncRNAs <- names(lncRNAs_expression)[lncRNAs_expression < max_gene_expression]
  filtered_PCGs <- names(PCGs_expression)[PCGs_expression < max_gene_expression]

  return(sce[c(filtered_lncRNAs,filtered_PCGs ),])
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

Filtering_TNBC <- function(sce, cells_mito_threshold, cells_max_threshold, cells_min_threshold, cells_min_genes_detected_threshold)
{
    total_min_expression <- sce$sum < cells_min_threshold
    total_max_mito_content <- sce$subsets_Mito_percent > cells_mito_threshold
    total_max_expression <- sce$sum > cells_max_threshold
    total_min_genes_detected <- sce$detected < cells_min_genes_detected_threshold

    discard <-  total_min_expression | total_max_expression | total_min_genes_detected | total_max_mito_content
    print(table(discard))
    print(DataFrame(total_max_expression=sum(total_max_expression),total_max_mito_content=sum(total_max_mito_content),total_min_genes_detected=sum(total_min_genes_detected), Total=sum(discard)))
    sce$discard <- discard
    #filter and normalize data
    sce_filt <- sce[,discard==F]

    sce_filt <- logNormCounts(sce_filt)

    return(sce_filt)
}

get_sequences <- function(features, transcripts)
{
  positions_features <- c()
  for (i in 1:length(features))
  {
    print(i)
    positions_features <- c(positions_features,grep(features[i],names(transcripts))) 
  }
  sequences <- transcripts[positions_features]
  return(sequences)
}


combine_4sce <- function(sce1, sce2, sce3, sce4, type_eval)
{
  combined_sce <- as.data.frame(cbind(c(sce1, sce2, sce3, sce4),c(rep("CellRanger",length(sce1)),rep("STARsolo",length(sce2)),rep("Kallisto",length(sce3)), rep("Salmon",length(sce4)))))
  colnames(combined_sce) <- c(type_eval,"ident")
  combined_sce$ident <- as.character(combined_sce$ident)
  combined_sce$ident <- factor(combined_sce$ident,  levels = c("CellRanger","STARsolo","Kallisto","Salmon"))
  combined_sce[,1] <- as.numeric(combined_sce[,1])
  return(combined_sce)
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

clustering <- function(data, k, hg38_ensembl_gtf, change_rownames = T)
{
  g <- buildSNNGraph(data, use.dimred = "PCA", k = k ) # k is the number of nearest neighbors used to construct the graph. Higher k
  clust <- igraph::cluster_louvain(g)$membership
  print(table(clust))
  data$louvain_clusters <- factor(clust)

  if (change_rownames == T)
  {
    # uniquifyFeatures
    gene_name <- hg38_ensembl_gtf$gene_name[match(rownames(data),hg38_ensembl_gtf$gene_id)]
    rowData(data) <- cbind(ens_id = rownames(data),gene_name)
    rownames(data) <- uniquifyFeatureNames(rownames(data), rowData(data)$gene_name)
  }

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
    kallisto_sce <- qc_metrics(pipeline = kallisto, mitochondrial_ens_ids = mito_gene)

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


crispr_data_intersection <- function(threshold_minumun_gene_counts_v,threshold_cells_detected_v, kallisto_sce_filt_clus, cellRanger_sce_filt_clus, STARsolo_sce_filt_clus, alevin_sce_filt_clus, lncrna_names,crispr_data,gene_name="gene_name",hg38_ensembl_gtf)
{
  final_df <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(final_df) <- c("features","gene_id","gene_name","Threshold")
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

    candidates_kallisto_lncRNAs <- intersect(candidates_kallisto, lncrna_names)
    common_lncRNAs <- intersect(common_genes, lncrna_names)

    if(gene_name=="gene_name")
    {
      candidates_kallisto_lncRNAs_ens_ids <- unique(hg38_ensembl_gtf$gene_id[hg38_ensembl_gtf$gene_name %in% candidates_kallisto_lncRNAs])
      candidates_kallisto_crispr <- intersect(gsub("\\..*","",candidates_kallisto_lncRNAs_ens_ids), gsub("\\..*","",crispr_data$gene_id))
      common_lncRNAs_ens_ids <- unique(hg38_ensembl_gtf$gene_id[hg38_ensembl_gtf$gene_name %in% common_lncRNAs])
      common_lncRNAs_crispr <- intersect(gsub("\\..*","",common_lncRNAs_ens_ids), gsub("\\..*","",crispr_data$gene_id))
    }

    else
    {
      candidates_kallisto_crispr <- intersect(gsub("\\..*","",candidates_kallisto_lncRNAs), gsub("\\..*","",crispr_data$gene_id))
      common_lncRNAs_crispr <- intersect(gsub("\\..*","",common_lncRNAs), gsub("\\..*","",crispr_data$gene_id))
    }
    candidates_kallisto_crispr_gene_name <- unique(hg38_ensembl_gtf$gene_name[gsub("\\..*","",hg38_ensembl_gtf$gene_id) %in% candidates_kallisto_crispr])
    common_lncRNAs_crispr_gene_name <- unique(hg38_ensembl_gtf$gene_name[gsub("\\..*","",hg38_ensembl_gtf$gene_id) %in% common_lncRNAs_crispr])

    all_lncRNAs <- data.frame(features = c(rep("exclusive_kallisto",length(candidates_kallisto_crispr)),rep("common",length(common_lncRNAs_crispr))), gene_id=c(candidates_kallisto_crispr,common_lncRNAs_crispr), gene_name = c(candidates_kallisto_crispr_gene_name,common_lncRNAs_crispr_gene_name))
    all_lncRNAs$Threshold = threshold_minumun_gene_counts_v[j]
    all_lncRNAs$candidates_kallisto_all_lncRNAs <- length(candidates_kallisto_lncRNAs)
    all_lncRNAs$common_all_lncRNAs <- length(common_lncRNAs)
    final_df <- rbind(final_df,all_lncRNAs)
  }
  final_df$Threshold <- as.factor(final_df$Threshold)

  # #hypergeometric_test
  hypergeometric_test_table <- table(final_df$Threshold,final_df$features)
  common_hyper_p_value_greater <- c()
  candidates_hyper_p_value_greater <- c()
  for (i in 1:nrow(hypergeometric_test_table))
  {
    common_hyper_p_value_greater <- c(common_hyper_p_value_greater, phyper(hypergeometric_test_table[i,1], 287, 16871-287, as.numeric(names(rev(table(final_df$common_all_lncRNAs)))[i]),lower.tail=F))
    candidates_hyper_p_value_greater <- c(candidates_hyper_p_value_greater, phyper(hypergeometric_test_table[i,2], 287, 16871-287, as.numeric(names(rev(table(final_df$candidates_kallisto_all_lncRNAs)))[i]),lower.tail=F))
  }
  hypergeometric_test_table <- cbind(hypergeometric_test_table,common_hyper_p_value_greater)
  hypergeometric_test_table <- cbind(hypergeometric_test_table,candidates_hyper_p_value_greater)

  to_return_list <- list(final_df, hypergeometric_test_table)
  return(to_return_list)
}

gene_names_sce <- function(sce, gtf)
{
  gene_name <- gtf$gene_name[match(gsub("\\..*","",rownames(sce)),gsub("\\..*","",gtf$gene_id))]
  rowData(sce) <- cbind(ens_id = rownames(sce),gene_name)
  rownames(sce) <- uniquifyFeatureNames(rownames(sce), rowData(sce)$gene_name)
  return(sce)
}

cell_type_specific_score_function <- function(sce,group_by, average_by)
{
  out <- sumCountsAcrossCells(logcounts(sce), sce[[group_by]], average=average_by)
  counts_cell = as.data.frame(assay(out))
  counts_cell <- counts_cell[rowSums(counts_cell[])>0,]
  counts_cell_specificity_index=t(apply(counts_cell, 1, function(x) x/(sum(x))))
  colnames(counts_cell_specificity_index) = out$ids

  cell_type_specificity_score <- c()
  for (i in 1:nrow(counts_cell_specificity_index))
  {
    s=0
    for (j in 1:ncol(counts_cell_specificity_index))
    {
      if(counts_cell_specificity_index[i,j]!=0)
      {
        s=s+(log(counts_cell_specificity_index[i,j],ncol(counts_cell_specificity_index))*counts_cell_specificity_index[i,j]) # Same results checked as the Shannon Entropy Specificity (HS): https://apcamargo.github.io/tspex/metrics/#fnref:9
      }
    }
    cell_type_specificity_score <- c(cell_type_specificity_score,1+s)
  }
  names(cell_type_specificity_score) <- rownames(counts_cell_specificity_index)

  return(list("cell_type_specificity_score" = cell_type_specificity_score,"counts_cell_specificity_index" = counts_cell_specificity_index))
}

# TNBC specific functions
filter_top_genes <- function(sce, features)
{
    sce <- sce[which(rownames(sce) %in% features)]
    return(sce)
}

uncorrected_integration_TNBC <- function(AA017, AA024, AA025, AA038, AA051)
{
    universe <- intersect(rownames(AA017), intersect(rownames(AA024), intersect(rownames(AA025), intersect(rownames(AA038),rownames(AA051) ))))
    length(universe)
    AA017 <- AA017[universe,]
    AA017 <- red_dim(AA017)
    AA024 <- AA024[universe,]
    AA024 <- red_dim(AA024)
    AA025 <- AA025[universe,]
    AA025 <- red_dim(AA025)
    AA038 <- AA038[universe,]
    AA051 <- AA051[universe,]
    rescaled <- multiBatchNorm(AA017, AA024, AA025, AA038,AA051)  # recomputes log-normalized expression values after adjusting the size factors for systematic differences in coverage between. This improves the quality of the correction by removing one aspect of the technical differences between batches
    AA017 <- rescaled[[1]]
    AA024 <- rescaled[[2]]
    AA025 <- rescaled[[3]]
    AA038 <- rescaled[[4]]
    AA051 <- rescaled[[5]]

    # No correction
    AA017$batch <- "AA017"
    AA024$batch <- "AA024"
    AA025$batch <- "AA025"
    AA038$batch <- "AA038"
    AA051$batch <- "AA051"

    universe_coldata <- intersect(colnames(colData(AA025)),intersect(colnames(colData(AA038)),colnames(colData(AA051))))
    colData(AA017) <- colData(AA017)[,universe_coldata]
    colData(AA024) <- colData(AA024)[,universe_coldata]
    colData(AA025) <- colData(AA025)[,universe_coldata]
    colData(AA038) <- colData(AA038)[,universe_coldata]
    colData(AA051) <- colData(AA051)[,universe_coldata]
    uncorrected <- cbind(AA017,AA024,AA025,cbind(AA038,AA051))
    uncorrected <- runPCA(uncorrected)
    uncorrected <- runTSNE(uncorrected, dimred = "PCA")
    uncorrected <- runUMAP(uncorrected, dimred = "PCA")
    snn.gr <- buildSNNGraph(uncorrected, use.dimred="PCA",k=30)
    clusters <- igraph::cluster_walktrap(snn.gr)$membership
    tab <- table(Cluster=clusters, Batch=uncorrected$batch)
    print(tab)
    uncorrected$combined_clustering <- as.factor(clusters)

    return(uncorrected)
}

mnn_integration_TNBC <- function(AA017, AA024, AA025, AA038, AA051)
{
    #MNN correction (OSCA)
    universe <- intersect(rownames(AA017), intersect(rownames(AA024), intersect(rownames(AA025), intersect(rownames(AA038),rownames(AA051) ))))
    length(universe)
    AA017 <- AA017[universe,]
    AA017 <- red_dim(AA017)
    AA024 <- AA024[universe,]
    AA024 <- red_dim(AA024)
    AA025 <- AA025[universe,]
    AA025 <- red_dim(AA025)
    AA038 <- AA038[universe,]
    AA051 <- AA051[universe,]
    rescaled <- multiBatchNorm(AA017, AA024, AA025, AA038,AA051)  # recomputes log-normalized expression values after adjusting the size factors for systematic differences in coverage between. This improves the quality of the correction by removing one aspect of the technical differences between batches
    AA017 <- rescaled[[1]]
    AA024 <- rescaled[[2]]
    AA025 <- rescaled[[3]]
    AA038 <- rescaled[[4]]
    AA051 <- rescaled[[5]]

    # No correction
    AA017$batch <- "AA017"
    AA024$batch <- "AA024"
    AA025$batch <- "AA025"
    AA038$batch <- "AA038"
    AA051$batch <- "AA051"

    universe_coldata <- intersect(colnames(colData(AA025)),intersect(colnames(colData(AA038)),colnames(colData(AA051))))
    colData(AA017) <- colData(AA017)[,universe_coldata]
    colData(AA024) <- colData(AA024)[,universe_coldata]
    colData(AA025) <- colData(AA025)[,universe_coldata]
    colData(AA038) <- colData(AA038)[,universe_coldata]
    colData(AA051) <- colData(AA051)[,universe_coldata]

    var_AA017 <- modelGeneVar(AA017)
    var_AA024 <- modelGeneVar(AA024)
    var_AA025 <- modelGeneVar(AA025)
    var_AA038 <- modelGeneVar(AA038)
    var_AA051 <- modelGeneVar(AA051)
    combined.dec <- combineVar(var_AA017, var_AA024, var_AA025, var_AA038, var_AA051)
    chosen.hvgs <- combined.dec$bio > 0
    sum(chosen.hvgs)
    #fastMNN to do MNN correction
    set.seed(100100100)
    mnn.out <- fastMNN(AA017, AA024,AA025,AA038,AA051, subset.row = chosen.hvgs, d = 50, k = 20) #don't subset genes
    mnn.out
    metadata(mnn.out)$merge.info$lost.var # percentage of variance lost after correction (have less than 10% [REF OSCA: http://bioconductor.org/books/3.16/OSCA.multisample/correction-diagnostics.html#mnn-specific-diagnostics])
    snn.gr <- buildSNNGraph(mnn.out, use.dimred="corrected", k = 10) #k=10
    clusters.mnn <- igraph::cluster_louvain(snn.gr)$membership
    tab.mnn <- table(Cluster=clusters.mnn, Batch=mnn.out$batch)
    tab.mnn
    mnn.out$louvain_clusters <- factor(clusters.mnn)
    mnn.out$batch <- factor(mnn.out$batch)
    mnn.out$identity <- (mnn.out$batch)
    levels(mnn.out$identity) <- c("AA017", "AA024","AA025","AA038","AA051")

    #dim reduction
    set.seed(0010101010)
    mnn.out <- runTSNE(mnn.out, dimred = "corrected")
    mnn.out <- runUMAP(mnn.out, dimred = "corrected")
    mnn.out_backup <- mnn.out
    colnames(mnn.out)[which(mnn.out$identity=="AA017")] <- paste(colnames(mnn.out)[which(mnn.out$identity=="AA017")],"AA017",sep="_")
    colnames(mnn.out)[which(mnn.out$identity=="AA024")] <- paste(colnames(mnn.out)[which(mnn.out$identity=="AA024")],"AA024",sep="_")
    colnames(mnn.out)[which(mnn.out$identity=="AA025")] <- paste(colnames(mnn.out)[which(mnn.out$identity=="AA025")],"AA025",sep="_")
    colnames(mnn.out)[which(mnn.out$identity=="AA038")] <- paste(colnames(mnn.out)[which(mnn.out$identity=="AA038")],"AA038",sep="_")
    colnames(mnn.out)[which(mnn.out$identity=="AA051")] <- paste(colnames(mnn.out)[which(mnn.out$identity=="AA051")],"AA051",sep="_")

    return(mnn.out)
}

plot_markers_integrated <- function(sce,mnn_object, markers, title, group)
{
  a=sce
  colnames(a)=colnames(mnn_object)
  print(DotPlot(as.Seurat(a), assay = "RNA", features = intersect(markers,rownames(sce)), scale.by = "size",col.min =  0, cols = c("lightgrey", "darkred"), group.by = c(group)) + ggtitle(title)  +geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+ theme(legend.position="bottom", legend.direction="horizontal")) 
}

to_seurat <- function(real_expression_object, mnn_object)
{
    universe=intersect(rownames(real_expression_object),rownames(mnn_object))
    b=real_expression_object[universe,]
    # We do that in order to color the UNCORRECTED expression!!! Very important! Never color CORRECTED values
    assay(mnn_object, withDimnames=FALSE,"counts")=assay(b, "counts")
    assay(mnn_object, withDimnames=FALSE,"logcounts")=assay(b, "logcounts")
    t=as.Seurat(mnn_object)
    return(t)
}

JIND_simplificate_cell_types <- function(JIND_raw_predictions)
{
  JIND_raw_predictions = gsub("B cells.*","B cells",JIND_raw_predictions)
  JIND_raw_predictions = gsub("CAFs.*","CAFs",JIND_raw_predictions)
  JIND_raw_predictions = gsub("Cancer.*","Cancer Epithelial",JIND_raw_predictions)
  JIND_raw_predictions = gsub("Endothelial.*","Endothelial",JIND_raw_predictions)
  JIND_raw_predictions = gsub("Macrophage","Myeloid",JIND_raw_predictions)
  JIND_raw_predictions = gsub("Monocyte","Myeloid",JIND_raw_predictions)
  JIND_raw_predictions = gsub("Mature Luminal","Epithelial",JIND_raw_predictions)
  JIND_raw_predictions = gsub("Luminal Progenitors","Epithelial",JIND_raw_predictions)
  JIND_raw_predictions = gsub("Myoepithelial","Epithelial",JIND_raw_predictions)
  JIND_raw_predictions = gsub("Cycling PVL","Mesenchymal",JIND_raw_predictions)
  JIND_raw_predictions = gsub("PVL.*","Mesenchymal",JIND_raw_predictions)

  return(JIND_raw_predictions)
}

COL2 = function(diverging = c('RdBu', 'BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdYlBu'),
                n = 200) {

  diverging = match.arg(diverging)

  colors = switch(
    diverging,
    RdBu = c('#67001F', '#B2182B', '#D6604D', '#F4A582', '#FDDBC7', '#FFFFFF',
             '#D1E5F0', '#92C5DE', '#4393C3', '#2166AC', '#053061'),
    BrBG = c('#543005', '#8C510A', '#BF812D', '#DFC27D', '#F6E8C3', '#FFFFFF',
             '#C7EAE5', '#80CDC1', '#35978F', '#01665E', '#003C30'),
    PiYG = c('#8E0152', '#C51B7D', '#DE77AE', '#F1B6DA', '#FDE0EF', '#FFFFFF',
             '#E6F5D0', '#B8E186', '#7FBC41', '#4D9221', '#276419'),
    PRGn = c('#40004B', '#762A83', '#9970AB', '#C2A5CF', '#E7D4E8', '#FFFFFF',
             '#D9F0D3', '#A6DBA0', '#5AAE61', '#1B7837', '#00441B'),
    PuOr = c('#7F3B08', '#B35806', '#E08214', '#FDB863', '#FEE0B6', '#FFFFFF',
             '#D8DAEB', '#B2ABD2', '#8073AC', '#542788', '#2D004B'),
    RdYlBu = c('#A50026', '#D73027', '#F46D43', '#FDAE61', '#FEE090', '#FFFFFF',
               '#E0F3F8', '#ABD9E9', '#74ADD1', '#4575B4', '#313695')
  )

  return(colorRampPalette(colors)(n))

}

uniquifyFeatures <- function(data, hg38_ensembl_gtf)
{
    gene_name <- hg38_ensembl_gtf$gene_name[match(rownames(data),hg38_ensembl_gtf$gene_id)]
    rowData(data) <- cbind(ens_id = rownames(data),gene_name)
    rownames(data) <- uniquifyFeatureNames(rownames(data), rowData(data)$gene_name)
    return(data)
}

uncorrected_integration_2_samples <- function(s1, s2, label1, label2)
{
    rescaled <- multiBatchNorm(s1, s2)
    s1 <- rescaled[[1]]
    s2 <- rescaled[[2]]
    # No correction
    s1$batch <- label1
    s2$batch <- label2
    universe_coldata <- intersect(colnames(colData(s1)),colnames(colData(s2)))
    colData(s1) <- colData(s1)[,universe_coldata]
    colData(s2) <- colData(s2)[,universe_coldata]

    uncorrected <- cbind(s1,s2)
    uncorrected <- red_dim(uncorrected)
    uncorrected <- clustering(uncorrected, k=20, "PCA", change_rownames = F)

    return(uncorrected)
}

mnn_integration_2_samples <- function(s1,s2, label1, label2)
{
    rescaled <- multiBatchNorm(s1, s2)
    s1 <- rescaled[[1]]
    s2 <- rescaled[[2]]
    # No correction
    s1$batch <- label1
    s2$batch <- label2
    universe_coldata <- intersect(colnames(colData(s1)),colnames(colData(s2)))
    colData(s1) <- colData(s1)[,universe_coldata]
    colData(s2) <- colData(s2)[,universe_coldata]

    #MNN correction (OSCA)
    var_s1 <- modelGeneVar(s1)
    var_s2 <- modelGeneVar(s2)
    combined.dec <- combineVar(var_s1, var_s2)
    chosen.hvgs <- combined.dec$bio > 0
    sum(chosen.hvgs)
    #fastMNN to do MNN correction
    set.seed(100100100)
    mnn.out <- fastMNN(s1, s2, subset.row = chosen.hvgs, d = 50, k = 20) 
    metadata(mnn.out)$merge.info$lost.var # percentage of variance lost after correction (have less than 10% [REF OSCA: http://bioconductor.org/books/3.16/OSCA.multisample/correction-diagnostics.html#mnn-specific-diagnostics])
    
    snn.gr <- buildSNNGraph(mnn.out, use.dimred="corrected", k = 10) #k=10
    clusters.mnn <- igraph::cluster_louvain(snn.gr)$membership
    tab.mnn <- table(Cluster=clusters.mnn, Batch=mnn.out$batch)
    tab.mnn
    
    mnn.out$louvain_clusters <- factor(clusters.mnn)
    mnn.out$batch <- factor(mnn.out$batch)
    mnn.out$identity <- (mnn.out$batch)
    levels(mnn.out$identity) <- c(label1, label2)

    #dim reduction
    set.seed(0010101010)
    mnn.out <- runTSNE(mnn.out, dimred = "corrected")
    mnn.out <- runUMAP(mnn.out, dimred = "corrected")
    colnames(mnn.out)[which(mnn.out$identity==label1)] <- paste(colnames(mnn.out)[which(mnn.out$identity==label1)],label1,sep="_")
    colnames(mnn.out)[which(mnn.out$identity==label2)] <- paste(colnames(mnn.out)[which(mnn.out$identity==label2)],label2,sep="_")

    return(mnn.out)
}

get_index_df <- function(a)
{
  new_ct <- c()
  for (i in 1:ncol(a))
  {
    new_ct <- c(new_ct,names(which.max(a[,i])))
  }
  return(new_ct)
}

simulate_reads <- function(lncrna_ens_ids, modified_gtf, protein_coding_ens_ids, file_name, CB_10X, simulation_iteration, genome_fasta)
{
    lncrna_subset <- sample(lncrna_ens_ids, 1500)
    lncRNA_gtf <- modified_gtf[modified_gtf$gene_id %in% lncrna_subset]
    protein_coding_subset <- sample(protein_coding_ens_ids, 1500)
    protein_coding_gtf <- modified_gtf[modified_gtf$gene_id %in% protein_coding_subset]
    subsetted_ids <- c(lncrna_subset, protein_coding_subset)

    coordinates_to_get_sequences <- GRanges(seqnames = "chr1", strand = c("+"),ranges = IRanges(start = c(1), width = 3))
    coordinates_to_get_sequences$gene_id <- "ENSG1"

    for (j in 1:3000)
    {
        print(j)
        exons_subset <- exons[exons$gene_id == subsetted_ids[j]]
        ranges(exons_subset)
        number_of_reads_exon = round(2300/length(exons_subset))
        for (k in 1:length(exons_subset))
        {
            possible_positions <- seq(start(exons_subset[k]), end(exons_subset[k])-92)
            new_start = sample(possible_positions, number_of_reads_exon, replace = T)
            chromosome <- rep(as.character(unique(seqnames(exons_subset))),length(new_start))
            strand = rep(unique(as.character(strand(exons_subset))),length(new_start))
            gr <- GRanges(seqnames = chromosome, strand = strand,ranges = IRanges(start = new_start, width = 91))
            gr$gene_id <- subsetted_ids[j]
            coordinates_to_get_sequences <- c(coordinates_to_get_sequences, gr)
        }
    }
    paste("/home/egonie/data/egonie/phd/test_data/Minju_ha_nat_comm/", file_name, sep="")
    saveRDS(list(coordinates_to_get_sequences, subsetted_ids),paste("/home/egonie/data/egonie/phd/test_data/Minju_ha_nat_comm/", file_name, sep=""))

    # remove first row
    genes_subsetted <- subsetted_ids
    coordinates_to_get_sequences <- coordinates_to_get_sequences[-1]
    # ADD CB
    CB <- sample(CB_10X$V1,100)
    CB_sequences <- sample(CB, length(coordinates_to_get_sequences), replace = T)
    coordinates_to_get_sequences$CB <- CB_sequences

    # Create REAL Matrix of expression that will be compared against Kallisto and Cell Ranger
    gene_ids <- unique(hg38_ensembl_gtf$gene_id)
    real_df <- matrix(0, nrow = length(gene_ids), ncol = length(CB))
    rownames(real_df) <- gene_ids
    counts <- table(coordinates_to_get_sequences$gene_id, coordinates_to_get_sequences$CB)
    colnames(real_df) <- colnames(counts)
    real_df[rownames(counts),] <- counts
    saveRDS(real_df, paste("simulation",simulation_iteration,"_real_count_matrix.rds",sep=""))

    #generate sequences
    sequences = BSgenome::getSeq(genome_fasta, coordinates_to_get_sequences)

    seqnames=seq(1452, 1452+length(sequences))
    seqnames_final <- paste("A00519:593:H7LNLDSXY:1:1101:",seqnames,":1000 2:N:0:TCCGTTGGAT+NCGTTCTCGC",sep="")
    names(sequences) <- seqnames_final[1:length(sequences)]

    
    writeXStringSet(sequences, paste("simulation",simulation_iteration,"_S1_L001_R2_001.fastq.gz",sep=""),append=FALSE,compress=TRUE, format="fastq")

    UMI_B <- replicate(length(coordinates_to_get_sequences),randDNA(12))
    sequence_R1 <- paste(CB_sequences, UMI_B, sep ="")
    sequence_R1_DNA <- DNAStringSet(sequence_R1)
    names(sequence_R1_DNA) <- seqnames_final[1:length(sequences)]
    writeXStringSet(sequence_R1_DNA, paste("simulation",simulation_iteration,"_S1_L001_R1_001.fastq.gz",sep=""),append=FALSE,compress=TRUE, format="fastq")
}
