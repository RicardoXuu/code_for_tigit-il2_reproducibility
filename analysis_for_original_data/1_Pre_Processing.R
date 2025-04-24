


# The raw data can be downloaded from the GSA database (https://ngdc.cncb.ac.cn/gsa/) with the accession number: CRA019549
# or contact the authors ( ustczxh@ustc.edu.cn )for the raw data


# Attention: The code works after getting the filtered data from the raw data by 10X Genomics Cell Ranger pipeline
# One need to run the Cell Ranger pipeline to get the filtered data first, then run the code below


# Code File 1 : Pre-processing the raw data

library(Seurat)
library(monocle3)
library(ggplot2)
library(stringr)


Cluster <- read.table("data/Analysis/4_Cell_clustering_analysis/4_1_all_cell_cluster/all_cells_features.csv",row.names =1,header = T,sep = ",")
length(unique(rownames(Cluster)))

Cell_ID_PBS <- gsub('_','-',rownames(Cluster)[Cluster$Samples%in%"PBS"])
Cell_ID_T2 <- gsub('_','-',rownames(Cluster)[Cluster$Samples%in%"T2"])


PBS_data <- Read10X(data.dir = "data/PBS/")
PBS_data <- PBS_data[,colnames(PBS_data)%in%Cell_ID_PBS]
dim(PBS_data)

T2_data <- Read10X(data.dir = "data/T2/") 
colnames(T2_data) <- gsub("-1","-2",colnames(T2_data))
T2_data <- T2_data[,colnames(T2_data)%in%Cell_ID_T2]
dim(T2_data)

rm(Cell_ID_PBS,Cell_ID_T2)


PBS <- CreateSeuratObject(counts = PBS_data, project = "PBS")
T2 <- CreateSeuratObject(counts = T2_data, project = "T2")
rm(PBS_data,T2_data)
gc()


seurat <- merge(PBS,T2)
rm(PBS,T2)

seurat@meta.data$Sample <- Cluster$Samples
seurat@meta.data$clusters <- Cluster$clusters

AnnoInfo <- read.table("data/Analysis/8_Reanno/all_cells_features.csv",row.names =1,header = T,sep = ",")
seurat@meta.data$anno <- AnnoInfo$clusters

rm(Cluster,AnnoInfo)






