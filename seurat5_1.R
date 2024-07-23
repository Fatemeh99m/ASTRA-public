#Loading libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(gdata)
library(reshape)
library(ggplot2)
library(data.table)
library(monocle3)
#reading data
#creating seurat object
#creating seurat object for each dataset :aging 1 to alz 8
aging1 <- CreateSeuratObject(aging1, min.cells = 3, min.features = 400, project = "aging 73_88")
#subsetting seurat objects
aging1 <- subset(aging1, subset = nFeature_RNA > 400 & nFeature_RNA < 5000)
###

###merging all seurat objects
all_samples.list <- c(aging1, aging2, aging3, aging4, aging5, aging6, aging7, aging8,
                  adult1, adult2,adult3, adult4, adult5,
                  adolescent1, adolescent2, adolescent3, adolescent4,
                  child1, child2, child3, child4,
                  infant1, infant2, infant3, infant4,
                  alz1, alz2, alz3, alz4, alz5, alz6, alz7, alz8)
names(all_samples.list) <- c("aging1", "aging2", "aging3", "aging4", "aging5", "aging6", "aging7", "aging8",
                         "adult1", "adult2", "adult3", "adult4", "adult5",
                         "adolescent1", "adolescent2", "adolescent3", "adolescent4",
                         "child1", "child2", "child3", "child4",
                         "infant1", "infant2", "infant3", "infant4",
                         "alz1", "alz2", "alz3", "alz4", "alz5", "alz6", "alz7", "alz8")
all_samples.combined <- merge(aging1, y = c(aging2, aging3, aging4, aging5, aging6, aging7, aging8,
                                adult1, adult2,adult3, adult4, adult5,
                                adolescent1, adolescent2, adolescent3, adolescent4,
                                child1, child2, child3, child4,
                                infant1, infant2, infant3, infant4,
                                alz1, alz2, alz3, alz4, alz5, alz6, alz7, alz8),
                  add.cell.ids = c("aging1", "aging2", "aging3", "aging4", "aging5", "aging6", "aging7", "aging8",
                                   "adult1", "adult2", "adult3", "adult4", "adult5",
                                   "adolescent1", "adolescent2", "adolescent3", "adolescent4",
                                   "child1", "child2", "child3", "child4",
                                   "infant1", "infant2", "infant3", "infant4",
                                   "alz1", "alz2", "alz3", "alz4", "alz5", "alz6", "alz7", "alz8"))



######data normalization
all_samples.combined <- NormalizeData(all_samples.combined)

######finding variable genes
all_samples.combined <- FindVariableFeatures(all_samples.combined, selection.method = "vst", nfeatures = 2000)

######scaling the data

all_samples.combined <- ScaleData(all_samples.combined, verbose = FALSE)
#####running PCA

all_samples.combined <- RunPCA(all_samples.combined, features = VariableFeatures(object = all_samples.combined), npcs = 100)

##elbow plot for finding number of pc
png("ElbPCA.png")
ElbowPlot(all_samples.combined,ndims = 100)
dev.off()

###JackStraw analysis for finding optimal number of pcs

all_samples.combined <- JackStraw(all_samples.combined,
                               dims=100,
                               num.replicate = 100)
all_samples.combined <- ScoreJackStraw(alz_heal.combined, dims = 1:100)
# determine pcs to use for downstream

pdf("jackplot1.pdf")
JackStrawPlot(all_samples.combined, dims = 1:50)
dev.off()
pdf("jackplot2.pdf")
jackplot2 <- JackStrawPlot(all_samples.combined, dims = 50:100)
dev.off()
#####################
#findneighbors
all_samples.combined <- FindNeighbors(all_samples.combined, dims = 1:100, reduction = "pca")
#find best resolution for clustering

resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
# create a list of Seurat objects for each resolution
seurat_list <- lapply(resolutions, function(resolution) {
  seurat <- FindClusters(all_samples.combined, resolution = resolution)
  return(seurat)
})
## do clustering
all_samples.combined <- FindClusters(all_samples.combined, resolution = 0.1, cluster.name = "unintegrated_clusters")

## runiing UMAP
all_samples.combined <- RunUMAP(all_samples.combined, dims = 1:100, reduction = "pca", reduction.name = "umap.unintegrated")
pdf("UMAP_ALL_res01.pdf")
DimPlot(all_samples.combined, reduction = "umap.unintegrated"و label=T)
dev.off()

#differential expression genes analysis for each cluster

#Using differentially expressed genes to help identify cell subsets
MIN_LOGFOLD_CHANGE = 0.25 # set to minimum required average log fold change in gene expression.
MIN_PCT_CELLS_EXPR_GENE = 0.1  # minimum percent of cells that must express gene in either clstr
all.markers = FindAllMarkers(all_samples.combined,
                             min.pct = MIN_PCT_CELLS_EXPR_GENE,
                             logfc.threshold = MIN_LOGFOLD_CHANGE,
                             slot = "data",
                             test.use = "wilcox",
                             return.thresh=0.05,
                             verbose = T,
                             only.pos = TRUE)
# sort all the markers by p-value
all.markers.sortedByPval = all.markers[order(all.markers$p_val),]

# take a look at the top most significant markers
head(all.markers.sortedByPval)
#a heatmap showing the top 10 marker genes for each cluster
pdf("topmarkerhitmap_allsamples.pdf")
top10 <- all.markers.sortedByPval %>%  group_by(cluster)  %>% do(head(., n=10))
DoHeatmap(object = all_samples.combined, features = top10$gene)+ NoLegend() +
    theme(axis.text.y = element_text(color = "black", size = 5))
dev.off()

#write an exel file of top 100 markers
library(writexl)
top100 <- all.markers.sortedByPval %>%  group_by("cca_clusters")  %>% do(head(., n=100))
write_xlsx(top100,"/home/top100marker_allsamples.xlsx")
#Get genes uniquely differentially expressed in each cluster
genes_uniquely_DE = all.markers.sortedByPval %>% dplyr::filter(avg_log2FC >= MIN_LOGFOLD_CHANGE) %>% group_by(gene) %>%  summarize(n=n()) %>%  dplyr::filter(n==1)
# Sort DE genes by p value.
genes_uniquely_DE.markers.sortedByPval =
  all.markers.sortedByPval[all.markers.sortedByPval$gene %in% genes_uniquely_DE$gene &
                             all.markers.sortedByPval$avg_log2FC >= MIN_LOGFOLD_CHANGE,]
top_marker_each = genes_uniquely_DE.markers.sortedByPval %>% dplyr::group_by("cca_clusters") %>% do(head(., n=100))
#write an exel file of top markers sorted
library(writexl)
write_xlsx(top_marker_each,"/home/top_marker_each_cluster_allsamples.xlsx")
###
saveRDS(all_samples.combined , file = "all_samples.rds")
