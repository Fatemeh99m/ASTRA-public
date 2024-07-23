#find the genes that are differentially expressed on the different paths through the trajectory
#using graph_test
#neighbor_graph="principal_graph",
#tells to test whether cells at similar positions on the trajectory have correlated expression
library(monocle)

cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=1)
pr_deg_ids <- row.names(subset(cds_pr_test_res, q_value < 0.05))

#finding trajectory-variable genes into modules
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=c(10^seq(-6,-1)))

#plot the aggregate module scores within each group of cell types or cell clusters
#we plot based on cca clusters

cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$cca_cluater)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")

#plot gene_modules related to clusters that are the subjects of study

plot_cells(cds,
           genes=gene_module_df %>% filter(module %in% c(2, 3, 9, 11)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

#filter gene modules to just show a specific module genes
library(dplyr)

# Assuming your dataframe is named 'df'
filtered_df_2 <- gene_module_df %>% 
  filter(module == "2") # Replace "2" with the desired module in your dataframe

#save filtered gene modules
write_xlsx(filtered_df_2,"/home/filtered_df_2.xlsx")

#save gene_module_df
write_xlsx(gene_module_df,"/home/gene_module_df.xlsx")