library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(patchwork)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(patchwork)
library(gridExtra)
library(scCustomize)
set.seed(1234)

##################################################################################################
get_earliest_principal_node <- function(cds,  cluster){
  print(paste0("In helper function for samples ",sample))
  
  cell_ids <- which( colData(cds)[, "clusters"] == cluster )
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  print(root_pr_nodes)
  
  root_pr_nodes
}
##################################################################################################
seurat.file="Notochord_20PCA_22_Nov_2025.rds"

seurat.obj <-readRDS(seurat.file)
DimPlot_scCustom(seurat.obj,pt.size = 1)
Idents(seurat.obj)<-seurat.obj$genotype
seurat.obj<-subset(seurat.obj,idents="control")
Idents(seurat.obj)<-seurat.obj$biol.names
DimPlot_scCustom(seurat.obj,pt.size = 1)
####################################################################################################
#prepare monocle
#https://rpubs.com/mahima_bose/Seurat_and_Monocle3_p
#cds <- as.cell_data_set(seurat.obj,assay="Xenium",group.by = 'seurat_clusters')
seurat.obj<-JoinLayers(seurat.obj)
cds <- as.cell_data_set(seurat.obj,assay="RNA",group.by = 'seurat_clusters')
colData(cds)$cluster <- cds@colData$seurat_clusters
fData(cds)$gene_short_name <- rownames(fData(cds))
##### assign partition
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
##########################
list.cluster <- seurat.obj@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
###########################################################
# UMAP
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- seurat.obj@reductions$umap@cell.embeddings
cluster.before.traj <-plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F, 
                                 group_label_size = 5) + theme(legend.position = "right")


#Learn Trajectory
cds <- learn_graph(cds, use_partition = T)
plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5)

cds <- order_cells(cds, reduction_method = "UMAP",  root_cells = colnames(cds[, clusters(cds) == "0"]))
p1<-plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = F, label_leaves = F,
           cell_size=1,graph_label_size = 3)+ggtitle("Pseudotime Trajectory, branch points are labeled")
p2<-plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
               label_branch_points = F, label_roots = F, label_leaves = T,
               cell_size=1,graph_label_size = 3)+ggtitle("Pseudotime Trajectory, final states are labeled")
p3<-DimPlot_scCustom(seurat.obj,group.by='RNA_snn_res.1')+ggtitle("clusters, res=1")
p4<-DimPlot_scCustom(seurat.obj,group.by='RNA_snn_res.1.5',label=F)+ggtitle("clusters, res=1.5")

grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
##########################################
# save pseudotime to seurat object
############################################
seurat.obj$pseudotime <- pseudotime(cds)
p1<-FeaturePlot_scCustom(seurat.obj,features = "pseudotime",pt.size=0.5)
p2<-DimPlot_scCustom(seurat.obj,pt.size=0.5,group.by = 'stage')
p1|p2
###################################################
#saveRDS(seurat.obj,file=paste(save.dir,"/notochord_subclusters_25Nov25.rds",sep=""))
# old










cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

osteoclust.combined.cds <- cluster_cells(cds = osteoclust.combined.cds, reduction_method = "UMAP")
osteoclust.combined.cds <- learn_graph(osteoclust.combined.cds, use_partition = TRUE)

## Calculate size factors using built-in function in monocle3
osteoclust.combined.cds<- estimate_size_factors(osteoclust.combined.cds)
## Add gene names into CDS
osteoclust.combined.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(osteoclust.c[["SCT"]])
# plot
plot_cells(osteoclust.combined.cds,
           color_cells_by = "orig.ident",
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE,
           cell_size = 2,
           group_label_size =4,
           graph_label_size=4,
           label_cell_groups = TRUE,
           label_branch_points=TRUE)

# order cells
osteoclust.combined_T.cds <- order_cells(osteoclust.combined.cds, root_pr_nodes=get_earliest_principal_node(osteoclust.combined.cds,"mono","1"))

plot_cells(osteoclust.combined_T.cds ,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_roots=TRUE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=5,
           cell_size = 0.5,
           show_trajectory_graph = TRUE)


# save trajectory
osteoclust.c<- AddMetaData(
  object = osteoclust.c,
  metadata = osteoclust.combined_T.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "traj_20PCA"
)
# plot trajectory
FeaturePlot(osteoclust.c, c("traj_20PCA"), pt.size = 0.1)& scale_color_viridis_c()
#
write.csv(osteoclust.c$traj_20PCA,file="pseudotime_20PCA.csv")

saveRDS(osteoclust.c, file = "osteoclust_30Nov21.rds")
