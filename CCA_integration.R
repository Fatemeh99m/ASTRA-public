#subset Granule cell, immature granule cells and astrocytes using wang 2022 and zhou 2022 cell annotation markers
zhou_markers <- c("SLC17A7", "EPHA5", "PROX1", "CALB1", "GAD1", "MOBP", "AQP4","TRHDE",
                  "SATB2", "PDGFRA", "APBB1IP", "TLE4", "RELN", "CLIC6", "FLT1", "FOXJ1")
wang_markers <- c("MPPED1", "SATB2", "PCP4", "RGS4", "CPNE4", "NETO1", "PROX1", "BCL11B",
                  "DCX", "SEMA3C", "GAD1", "GAD2", "AQP4", "MOBP", "MOG", "PDGFRA", "OLIG2",
                  "CX3CR1", "PTPPRC", "FLI1", "FN1", "DNAH9", "GFAP54", "DCN", "COL1A2")

# cluster 3 belongs to astrocytes
all_samples_sub <- subset(x = all_samples.combined, idents = c("0", "3", "4", "9"))

#integration and reclustering to find subclusters
#CCA
all_samples_sub <- IntegrateLayers(
  object = all_samples_sub, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

all_samples_sub <- FindNeighbors(all_samples_sub, reduction = "integrated.cca", dims = 1:100)

# create a vector of different resolutions of all_samples_sub object
resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
# create a list of Seurat objects for each resolution
seurat_list <- lapply(resolutions, function(resolution) {
  seurat <- FindClusters(all_samples_sub, resolution = resolution)
  return(seurat)
})
# res 0.3 for clustering
all_samples_sub <- FindClusters(all_samples_sub, resolution = 0.3, cluster.name = "cca_clusters")


#Runing UMAP
all_samples_sub <- RunUMAP(all_samples_sub, reduction = "integrated.cca", dims = 1:50, reduction.name = "umap.cca")
pdf("integrated_all_res03.pdf")
DimPlot(
  alz,
  reduction = "umap.cca",
  group.by = "cca_clusters",
  combine = FALSE, label.size = 3, label=T
)
dev.off()

all_samples_sub <- JoinLayers(all_samples_sub)
