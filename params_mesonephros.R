save.dir='/home/labs/olenderlab/lvzvia/zelzer/Mor_grinshtein/analysis/limb_buds_aug25'
pca.num=20
res.num=0.8
normalization.method <- 'RNA'
folder.1 <- '1_QC'
folder.2<-'2_markers'
folder.3<-'3_automatic_annotation'
folder.4<-'4_comparison_across_conditions'
object.reduction="pca"
plot_users_genes=0
test_conditions=0
###############
csv_files_path="/home/labs/olenderlab/lvzvia/zelzer/Mor_grinshtein/data/limb_buds"
metadata_all<-read.csv(file="limb_buds_info.csv",header=T)
experiment_name="limb_buds" #for saving the object
#