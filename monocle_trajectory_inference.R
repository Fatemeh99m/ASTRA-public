obj <- all_samples.combined

gene_annotation <- as.data.frame(rownames(obj@reductions[["pca"]]@feature.loadings), row.names = rownames(obj@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"
head(gene_annotation)
#
cell_metadata <- as.data.frame(obj[["RNA"]]$counts@Dimnames[[2]], row.names = obj[["RNA"]]$counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"
head(cell_metadata)
#
New_matrix <- obj[["RNA"]]$data
New_matrix <- New_matrix[rownames(obj@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix
New_matrix[1:10,1:10]
#
cds <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)
#
### Construct and assign the made up partition
recreate.partition <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition
#
list_clusters <-Idents(obj)
cds@colData$cca_clusters <- list_clusters[rownames(cds@colData)]
names(list_clusters) <- obj[["RNA"]]$data@Dimnames[[2]]
cds@clusters@listData[["UMAP"]][["cca_clusters"]] <- list_clusters
#
cds@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
cds@int_colData@listData$reducedDims$UMAP <-obj@reductions[["umap.cca"]]@cell.embeddings
#

cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds, verbose = FALSE, learn_graph_control = list(geodesic_distance_ratio=1))

plot_cells(cds,
           color_cells_by = "cca_clusters",
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE, , group_label_size = 3)


cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)