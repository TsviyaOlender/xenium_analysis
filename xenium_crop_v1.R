library(Seurat)
library(sf)
library(ggplot2)
library(dplyr)
library(patchwork)
library(scCustomize)
library(jsonlite)
library(gridExtra)
library(grid)
###################################################################################
source("/home/labs/olenderlab/lvzvia/MyRScripts/scRNA_pipline/xenium_functions_V1.R")
###################################################################################
qc_params=c("nCount_RNA","nFeature_RNA","percent_apop","percent_dna_repair","percent_oxphos", "percent_hemo")
# RUN params file!!!
######################################################################################
# create subfolders for the plots
solution.folder=create_folders(save.dir,pca.num,res.num)
folder.1=paste(save.dir,"/1_QC_plots",sep="")
folder.2=paste(solution.folder,"/2_markers",sep="")
folder.3=paste(solution.folder,"/3_automatic_annotation",sep="")
folder.4=paste(solution.folder,"/4_comparison_across_conditions",sep="")
##################################################################
xenium.obj.1<-create_expression_object("/home/labs/olenderlab/lvzvia/zelzer/Mor_grinshtein/data/output-XETG00188_0037810_Region_1_20250108_143713",
                                       folder.1,"region1")
xenium.obj.3<-create_expression_object("/home/labs/olenderlab/lvzvia/zelzer/Mor_grinshtein/data/output-XETG00188_0037810_Region_3_20250108_143713",
                                       folder.1,"region3")
xenium.obj.4<-create_expression_object("/home/labs/olenderlab/lvzvia/zelzer/Mor_grinshtein/data/output-XETG00188_0037810_Region_4_20250108_143713",
                                       folder.1,"region4")
xenium.obj.5<-create_expression_object("/home/labs/olenderlab/lvzvia/zelzer/Mor_grinshtein/data/output-XETG00188_0037810_Region_5_20250108_143713",
                                       folder.1,"region5")
xenium.obj.6<-create_expression_object("/home/labs/olenderlab/lvzvia/zelzer/Mor_grinshtein/data/embryo_4/output-XETG00188_0029402_Region_1_20250612_101514",
                                       folder.1,"region6")
########################################################################################
#read metadata
for(i in 1:nrow(metadata_all)){
  metadata_all$obj_name[i]=paste(metadata_all$sub_structure[i],"_",metadata_all$Number[i],
                                 "_",metadata_all$side[i],sep="")
}

########################################################################################
####################################################################
seurat_list <- list(xenium.obj.1,xenium.obj.3,xenium.obj.4,xenium.obj.5,xenium.obj.6)
region_i <- c(1, 3, 4, 5, 6)
#seurat_list<-list(xenium.obj.6)
#region_i <- c(6)
for (j in seq_along(region_i)) {
  # filter metadata for this region
  metadata <- subset(metadata_all, region == region_i[j])
  message("Region: ", region_i[j], " | entries: ", nrow(metadata))
  if(nrow(metadata) > 0){
    # pull the correct Seurat object from the list
    xenium.obj <- seurat_list[[j]]
    
    # loop over rows in metadata
    for (i in seq_len(nrow(metadata))) {
      file <- file.path(csv_files_path, metadata$file[i])
      print(file)
      cells_to_extract <- read.table(file, header = TRUE, sep = ",", skip = 2)
      subset_obj <- subset(xenium.obj, cells = cells_to_extract$Cell.ID)
      # add metadata columns (recycling automatically expands to all cells)
      subset_obj$orig.ident <- metadata$sub_structure[i]
      subset_obj$embryo     <- metadata$embryo[i]
      subset_obj$region     <- metadata$region[i]
      subset_obj$side       <- metadata$side[i]
      subset_obj$stage      <- metadata$stage[i]
      subset_obj$name       <- metadata$obj_name[i]
      
      message(" -> Created subset: ", metadata$obj_name[i], " (", ncol(subset_obj), " cells)")
      assign(metadata$obj_name[i], subset_obj, envir = .GlobalEnv)
    }
  }
}
#somites
#experiments=c('EM','PSM','ES','DMScl','DS','DTS')
#notochord.1
experiments=c("Early","Intermediate","Late_brachial","Late_trunk","Late")   
#mesonephros
#experiments<-c("EM","ENS","IM","LNBL","LNBR","LNS","LNTL","LNTR")
for(i in 1:length(experiments)){
  print(i)
  metadata=subset(metadata_all,metadata_all$sub_structure==experiments[i])
  seurat_objs <- lapply(metadata$obj_name, get)
  # the below is a function of ScCustomize. Must be installed
  seurat.obj <-Merge_Seurat_List(seurat_objs,add.cell.ids = metadata$obj_name,merge.data = TRUE)
  assign(experiments[i],seurat.obj)
}
#experiments=c("limb_114_R","limb_115_L","limb_116_R","limb_117_R")
seurat_objs <- lapply(experiments, get)
seurat.obj <-Merge_Seurat_List(seurat_objs,merge.data = TRUE)
seurat.obj$name_short <- sub("_[^_]+$", "", seurat.obj$name)
seurat.obj <- JoinLayers(seurat.obj)
VlnPlot(seurat.obj  , features = c("nFeature_RNA", "nCount_RNA"), ncol = 3,pt.size=0.05)
##################################################################

##################################################################
#normalize, re-scale
seurat.obj<-NormalizeData(seurat.obj,scale.factor=median(xenium.obj$nCount_RNA))
VariableFeatures(seurat.obj) <- rownames(seurat.obj) 
seurat.obj <- ScaleData(seurat.obj,features = row.names(seurat.obj))
seurat.obj <- RunPCA(seurat.obj)
elbow_plot<-ElbowPlot(seurat.obj,ndims=50)
file.name <- paste(folder.1,"/elbow_plot.jpg",sep="")
ggsave(file.name, plot = elbow_plot, width = 15, height = 10, dpi = 150)
#####################################################
#Score cc
#cell cycle
calc.cc=1
if(calc.cc==1){
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  seurat.obj<-CellCycleScoring(seurat.obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  #seurat.obj <- ScaleData(seurat.obj,vars.to.regress = c("S.Score", "G2M.Score"), features = row.names(seurat.obj))
  DimPlot_scCustom(seurat.obj,group.by = 'Phase')
  seurat.obj <- RunPCA(seurat.obj, features = c(s.genes, g2m.genes))
  DimPlot(seurat.obj,reduction = 'pca')
  elbow_plot<-ElbowPlot(seurat.obj,ndims=50)
  file.name<-paste(folder.1,"/elbow_plot_after_cc.jpg",sep="")
  ggsave(file.name, plot = elbow_plot, width = 15, height = 10, dpi = 150)
}
# If need to regress cc
calc.cc=1
if(calc.cc==1){
  seurat.obj <- ScaleData(seurat.obj,vars.to.regress = c("S.Score", "G2M.Score"), features = row.names(seurat.obj))
  seurat.obj <- RunPCA(seurat.obj,features = row.names(seurat.obj))
}
################################################################################
## get this list as a vector and save to disk for QC
# Get the PCA loadings
top_genes=get_top_PC_genes(seurat.obj)#in the seurat_functions file
top_genes=unique(top_genes)
plot.loading=VizDimLoadings(seurat.obj, dims = 1:2, reduction = "pca")
ggsave(paste(folder.1,"/VizDimLoadings.png",sep=""), plot = plot.loading, width = 15, height = 10, dpi = 150)
DimHeatmap(seurat.obj, dims = 25:35, cells = 1000, balanced = TRUE)
# save
top_genes_f=paste(folder.1,"/top_genes_loadings_",pca.num,"_PCA.csv",sep="")
write.csv(top_genes,file=top_genes_f)
#plot
top_BP_Plot<-singleList_enrichR(top_genes)
#file.name<-paste(folder.1,"/GO_BP_top_PCA_genes.jpg",sep="")
#ggsave(file.name, plot = top_BP_Plot, width = 15, height = 10, dpi = 150)
##########################################################################
# UMAP reduction
set.seed(1234)
seurat.obj <- RunUMAP(seurat.obj, dims = 1:pca.num, verbose = FALSE)
################################################################################
p1<-DimPlot_scCustom(seurat.obj,label=T,  group.by = 'Phase',label.size = 6,repel = T)
p2<-DimPlot_scCustom(seurat.obj,label=T,  group.by = 'stage',label.size = 6,repel=T)
p3<-DimPlot_scCustom(seurat.obj,label=T,  group.by = 'orig.ident',label.size = 6,repel=T)
p4<-DimPlot_scCustom(seurat.obj,label=T,  group.by = 'side',label.size = 6,repel=T)
multi_umap=grid.arrange(p1,p2,p3,p4, ncol = 2)
ggsave(paste(folder.1,"/multi_umap.png",sep=""), plot = multi_umap, width = 15, height = 13, dpi = 150)
###################################################################################
res_grid <- c(1,1.5)

for (res.num in res_grid) {
  solution.folder <- create_folders(save.dir, pca.num, res.num)
  folder.2 <- file.path(solution.folder, "2_markers")
  folder.3 <- file.path(solution.folder, "3_automatic_annotation")
  folder.4 <- file.path(solution.folder, "4_comparison_across_conditions")
  
  #8. cluster cells
  to.save=1# if 0 the UMAP plot will not be saved to disk. 
  #if 1= the UMAP plot will be saved.
  # this is good at the beginning of the project when we try several criteria
  seurat.obj<-cluster_cells(seurat.obj,pca.num,res.num,folder.2,to.save,'pca')
  markers.file = paste(folder.2,"/","cell_markers_",pca.num,"PCA_res_",res.num,".csv",sep="")
  find_markers(seurat.obj,folder.2,normalization.method,pca.num,res.num)
  ## plot clustering results in UMAP
  
  ###############################################################################
  qc_per_cluster=VlnPlot_scCustom(seurat.obj,features = c("nFeature_RNA","percent_oxphos","percent_apop",
                                           "percent_dna_repair","percent_ieg","percent_hemo"),pt.size = 0.1)
  cell_count_proportions=Proportion_Plot(seurat_object = seurat.obj, plot_type = "bar", split.by = "Phase", plot_scale = "count")
  cell.count <- Cluster_Stats_All_Samples(seurat_object = seurat.obj, group_by_var = "Phase")
  cell.count.file = paste(folder.2,"/","cell.count_stage_byPhase_",pca.num,"PCA_res_",res.num,".csv",sep="")
  write.csv(cell.count,file=cell.count.file)
  ggsave(paste(folder.2,"/cell_count_proportions.png",sep=""), plot = cell_count_proportions, width = 15, height = 8, dpi = 150)
  ggsave(paste(folder.2,"/QC_per_cluster.png",sep=""), plot = qc_per_cluster, width = 20, height = 15, dpi = 150)
  ##############################################################################
  # export clustering for Xenium explorer
  # generate a table with two columns - one for cell id and the other for group (clustering)
  clustering_res <- FetchData(seurat.obj, vars = c("seurat_clusters"))
  clustering_res$cell_id <- rownames(clustering_res)
  clustering_res <- clustering_res[,c("cell_id","seurat_clusters")]
  clustering_res$cell_id <- sub(".*_", "", clustering_res$cell_id)
  row.names(clustering_res) <- clustering_res$cell_id 
  head(clustering_res, n = 3)
  clustering_cell_barcodes_file = paste(folder.2,"/","cell.barcodes_clusters_",pca.num,"PCA_res_",res.num,".csv",sep="")
  write.csv(clustering_res,file=clustering_cell_barcodes_file)
}
#################################################################################
# QC final
p1<-VlnPlot(seurat.obj,features = "nFeature_RNA",group.by = "RNA_snn_res.1.5",pt.size=0.1)
p2<-DimPlot_scCustom(seurat.obj,group.by = "RNA_snn_res.1.5",label=T)
p1|p2
#################################################################################
#save
formatted_date <- format(Sys.Date(), "%d_%b_%Y")
seuratFileName = paste(save.dir,"/",experiment_name,"_",pca.num,"PCA_",formatted_date,".rds",sep="")
saveRDS(seurat.obj,file=seuratFileName)
session_file_name=paste(save.dir,"/session_info.txt",sep="")
sessionInfo() %>% capture.output(file=session_file_name)
#################################################################################
#additional
#plot UMAP by cut number
library(randomcoloR)

n_groups <- length(levels(factor(seurat.obj$name_short)))
cols <- distinctColorPalette(n_groups)
p1<-DimPlot_scCustom(seurat.obj,colors_use = cols,group.by = "name_short")
ggsave(paste(folder.1,"/UMAP_by_slice.png",sep=""), plot = p1, width = 20, height = 17, dpi = 150)

#################################################################################
#################################################################################
#Highlight cells from previous cluster10 in current analysis
cells_by_cluster <- lapply(levels(analysis1$seurat_clusters), function(cl) {
  cells <- WhichCells(analysis1, idents = cl)
  intersect(cells, Cells(seurat.obj))
})
cells_by_cluster<-function(cluster_n){
  cells_c10 <- WhichCells(analysis1, idents = cluster_n)
  cells_c10 <- intersect(cells_c10, Cells(seurat.obj))
  cells0 <- list(cells_c10 = cells_c10)
 return(cells0)
}

p1<-DimPlot_scCustom(seurat.obj)+ggtitle("current, Res 1.0")
p2<-Cell_Highlight_Plot(seurat_object = seurat.obj, cells_highlight = cells_by_cluster("9"),pt.size = 0.2)+ggtitle("Presomitic mesoderm")
p3<-Cell_Highlight_Plot(seurat_object = seurat.obj, cells_highlight = cells_by_cluster("10"),pt.size = 0.2)+ggtitle("Lateral plate mesoderm")
p4<-Cell_Highlight_Plot(seurat_object = seurat.obj, cells_highlight = cells_by_cluster("12"),pt.size = 0.2)+ggtitle("Ectorderm")
grid.arrange(p1,p2,p3,p4, ncol = 2)



p4<-FeaturePlot_scCustom(seurat.obj,features="CDH3")+ggtitle("CDH3, current")



FeaturePlot_scCustom(analysis1,features="CDH3")
