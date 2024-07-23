library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)
library(KEGGREST)
library(readxl)
library(enrichplot)

# Load marker gene data
file_path <- "path/top_marker_each_cluster_allsamples.xlsx"
marker_genes <- read_excel(file_path)

# Filter data for cluster 9
cluster_9_genes <- marker_genes$gene[marker_genes$cluster == 9]

# GO Enrichment Analysis
go_results <- enrichGO(gene = cluster_9_genes,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)

# KEGG Pathway Enrichment Analysis
kegg_results <- enrichKEGG(gene = cluster_9_genes,
                           organism = 'hsa',
                           keyType = 'kegg',
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)

# Visualize and save results
pdf("GO_Enrichment_Results_Cluster_9.pdf")
dotplot(go_results, showCategory = 20) + ggtitle("GO Enrichment for Cluster 9")
dev.off()

pdf("KEGG_Enrichment_Results_Cluster_9.pdf")
dotplot(kegg_results, showCategory = 20) + ggtitle("KEGG Enrichment for Cluster 9")
dev.off()

# Save GO Results
write.csv(as.data.frame(go_results), file = "GO_Cluster_9.csv")

# Save KEGG Results
write.csv(as.data.frame(kegg_results), file = "KEGG_Cluster_9.csv")
#barchart
barplot(go_results, showCategory = 20) + ggtitle("GO Enrichment Bar Chart for Cluster 9")

# Network Plot for GO Enrichment Results
pdf("GO_Enrichment_Network_Cluster_9.pdf")
cnetplot(go_results, showCategory = 20) + ggtitle("GO Enrichment Network for Cluster 9")
dev.off()