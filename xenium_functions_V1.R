library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggplotify)
library(gridExtra)
library(scCustomize)
find_markers<-function(the.object,save.dir,normalization.method,pca.num,res.num){
  if(normalization.method == "SCT"){ 
    the.object <- PrepSCTFindMarkers(the.object)
  }else if(normalization.method == "RNA"){
    the.object <- JoinLayers(the.object)
  }
  the.object.markers <- FindAllMarkers(the.object, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2) %>%
    Add_Pct_Diff()#,recorrect_umi=FALSE
  
  the.object.markers <-subset(the.object.markers,the.object.markers$p_val_adj<0.05)
  #sort
  the.object.markers <- the.object.markers %>%
    arrange(cluster, desc(avg_log2FC))
  #TOP10
  top10 <- the.object.markers %>%
    filter(avg_log2FC > 0) %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 10)
  #top5
  top5 <- the.object.markers %>%
    filter(avg_log2FC > 0) %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 5)
  
  HP2<-DoHeatmap(the.object, features = c(top5$gene)) + NoLegend()
  # save. file names
  
  top5.file = paste(save.dir,"/","top5.markers_",pca.num,"PCA_res_",res.num,".csv",sep="")
  top10.file = paste(save.dir,"/","top10.markers_",pca.num,"PCA_res_",res.num,".csv",sep="")
  top5.HM.file = paste(save.dir,"/","top5.markers.heatmap_",pca.num,"PCA_res_",res.num,".png",sep="")
  
  # write
  markers.file=paste(save.dir,"/","cell_markers_",pca.num,"PCA_res_",res.num,".csv",sep="")
  write.csv(the.object.markers,file=markers.file)
  write.csv(top5,file=top5.file)
  write.csv(top10,file=top10.file)
  ggplot2::ggsave(filename = top5.HM.file, plot = HP2,units ="cm",width=30,height=20)
  
}
cluster_cells<-function(the.object,pcas,res,save.dir,to.save,object.reduction){
  set.seed(1234)
  the.object <- FindNeighbors(the.object, dims = 1:pcas,reduction =object.reduction)
  the.object <- FindClusters(the.object, resolution = res )
  p1<-DimPlot_scCustom(the.object,label.size = 6, color_seed = 2)
  p2<-DimPlot_scCustom(the.object,group.by = 'orig.ident')
  print (p1|p2)
  if(to.save == 1){
    umap.file=paste(save.dir,"/umap_",pcas,"PCAs_",res,"res.png",sep='')
    ggplot2::ggsave(filename = umap.file, plot = p1|p2,width = 20, height = 15)
  }
  return(the.object)
}
get_top_PC_genes<-function(the.object){
  pca_results <- the.object[["pca"]]
  loadings <- pca_results@feature.loadings
  top_genes <- loadings[, 1:5] %>% as.data.frame()  # Convert to data frame
  top_genes<-as.data.frame(top_genes)
  top_genes$names=row.names(top_genes)
  # Sort the data frame by the second column in descending order and extract names
  top_n <- 10 # Specify how many top values you want to extract
  top_names<-vector()
  
  for(i in 1:5){
    pc.num=paste("PC_",i,sep="")
    # Sort the data frame by the 'Age' column in descending order
    top_genes <- top_genes %>%
      arrange(desc(!!sym(pc.num)))  # Use desc() for descending order
    top_names<-c(top_names,top_genes$names[1:top_n])
    # Sort the data frame by the 'Age' column in ascending order
    top_genes <- top_genes %>%
      arrange(!!sym(pc.num))  # Use sym and !! to refer to the column
    top_names<-c(top_names,top_genes$names[1:top_n])
  }
  
  return(top_names)
  
}
create_folders <-function(save.dir,pca.num,res.num){
  #create main folder
  ifelse(!dir.exists(save.dir), dir.create(save.dir, showWarnings = FALSE), FALSE)
  # cerate folders for the results
  folder.1=paste(save.dir,"/1_QC_plots",sep="")
  ifelse(!dir.exists(folder.1), dir.create(folder.1, showWarnings = FALSE), FALSE)
  #
  solution.folder=paste(save.dir,"/solution_pca",pca.num,"_res",res.num,sep="")
  ifelse(!dir.exists(solution.folder), dir.create(solution.folder, showWarnings = FALSE), FALSE)
  #
  folder.2=paste(solution.folder,"/2_markers",sep="")
  ifelse(!dir.exists(folder.2), dir.create(folder.2, showWarnings = TRUE), FALSE)
  #
  folder.3=paste(solution.folder,"/3_automatic_annotation",sep="")
  ifelse(!dir.exists(folder.3), dir.create(folder.3, showWarnings = FALSE), FALSE)
  #
  folder.4=paste(solution.folder,"/4_comparison_across_conditions",sep="")
  ifelse(!dir.exists(folder.4), dir.create(folder.4, showWarnings = FALSE), FALSE)
  return(solution.folder)
}
Do_QC_plots<-function(the.object,folder.1,plot_name){
  qc_params=c("nCount_RNA","nFeature_RNA","percent_apop","percent_dna_repair","percent_oxphos", "percent_hemo")
  vln_1=VlnPlot_scCustom(the.object,features=qc_params,pt.size=0)
  num_cells_report<-as.data.frame(table(the.object$orig.ident))
  text <- tableGrob(num_cells_report)
  tbl_plot <- as.ggplot(text)
  grid_plot.1=vln_1 / tbl_plot + plot_layout(widths = c(3, 1))
  file1.name<-paste(folder.1,"/",plot_name,sep="")
  ggsave(file1.name, plot = grid_plot.1, width = 20, height = 30, dpi = 150)
}
create_expression_object<-function(data.dir,folder.1,obj_name){
  expression.matrix=paste(data.dir,"/","cell_feature_matrix.h5",sep="")
  data<-Read10X_h5(expression.matrix)
  expression.obj<-CreateSeuratObject(counts = data[["Gene Expression"]])
  expression.obj<-Add_Cell_QC_Metrics(expression.obj, species = "human")
  thres <- quantile(expression.obj$nCount_RNA, c(0.98))
  expression.obj.f <- subset(expression.obj, subset = nCount_RNA >= 30 & nCount_RNA <= thres[1])
  #make plots
  Do_QC_plots(expression.obj,folder.1,paste(obj_name,"_before_filter.jpeg",sep=""))
  Do_QC_plots(expression.obj.f,folder.1,paste(obj_name,"_after_filter.jpeg",sep=""))
  return(expression.obj.f)
}