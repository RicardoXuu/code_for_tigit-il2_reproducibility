

# Code File 3 : Analysis of Treg cells from CD4 cells

library(Seurat)
library(monocle3)
library(ggplot2)
library(stringr)


cds_treg <- cds[,colnames(cds)%in%treg_cell_id]
cds_treg <- cluster_cells(cds_treg, resolution=15e-4)
pdf("plot/Treg_Monocle3_anno.pdf",height = 3, width = 3)
plot_cells(cds_treg,
           group_label_size = 2,
           color_cells_by = "anno",
           group_cells_by = "anno")
dev.off()

cds_treg_PBS <- cds_treg[,pData(cds_treg)$Sample%in%"PBS"]
pdf("plot/Treg_Monocle3_anno_PBS.pdf",height = 3, width = 3)
plot_cells(cds_treg_PBS,
           group_label_size = 2,
           color_cells_by = "anno",
           group_cells_by = "anno")
dev.off()
rm(cds_treg_PBS)

cds_treg_T2 <- cds_treg[,pData(cds_treg)$Sample%in%"T2"]
pdf("plot/Treg_Monocle3_anno_T2.pdf",height = 3, width = 3)
plot_cells(cds_treg_T2,
           group_label_size = 2,
           color_cells_by = "anno",
           group_cells_by = "anno")
dev.off()
rm(cds_treg_T2)


pdf("plot/Treg_Monocle3_NewCluster.pdf",height = 3, width = 3)
plot_cells(cds_treg,
           group_label_size = 4,
           color_cells_by = "cluster")
dev.off()

cds_treg_PBS <- cds_treg[,pData(cds_treg)$Sample%in%"PBS"]
pdf("plot/Treg_Monocle3_NewCluster_PBS.pdf",height = 3, width = 3)
plot_cells(cds_treg_PBS,
           group_label_size = 4,
           color_cells_by = "cluster")
dev.off()
rm(cds_treg_PBS)

cds_treg_T2 <- cds_treg[,pData(cds_treg)$Sample%in%"T2"]
pdf("plot/Treg_Monocle3_NewCluster_T2.pdf",height = 3, width = 3)
plot_cells(cds_treg_T2,
           group_label_size = 4,
           color_cells_by = "cluster")
dev.off()
rm(cds_treg_T2)



pData(cds_treg)$highRes <- cds_treg@clusters@listData[["UMAP"]][["clusters"]]
treg_cell_cluster <- as.data.frame(pData(cds_treg))
treg_cell_cluster <- cbind(rownames(treg_cell_cluster),treg_cell_cluster)
rownames(treg_cell_cluster) <- 1:nrow(treg_cell_cluster)
colnames(treg_cell_cluster)[1] <- c("cell_id")
write.table(treg_cell_cluster, file = "data/treg_cell_cluster.txt",col.names = T,row.names = F,quote = F,sep = "\t")












