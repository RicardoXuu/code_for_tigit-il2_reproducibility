

# Code File 2 : Analysis of total CD4 cells 

# Monocle3 
library(Seurat)
library(monocle3)
library(ggplot2)
library(stringr)


seurat_cd4 <- seurat[,seurat$anno%in%unique(seurat@meta.data$anno)[1:2]]
seurat_cd4 <- JoinLayers(seurat_cd4)

matrix <- seurat_cd4@assays[["RNA"]]@layers[["counts"]]
matrix <- as.matrix(matrix)
cell_meta <- seurat_cd4@meta.data
gene_anno <- data.frame(gene_short_name = rownames(seurat_cd4))
rownames(gene_anno) <- gene_anno$gene_short_name
cds <- monocle3::new_cell_data_set(matrix,
                                   cell_metadata = cell_meta,
                                   gene_metadata = gene_anno)

rm(matrix,cell_meta,gene_anno)
rm(seurat_cd4)
gc()


cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds)

pdf("plot/CD4_Monocle3_anno.pdf",height = 4, width = 4)
plot_cells(cds, 
           group_label_size = 4,
           label_groups_by_cluster=FALSE,  
           color_cells_by = "anno")
dev.off()

cds_PBS <- cds[,pData(cds)$Sample%in%"PBS"]
pdf("plot/CD4_Monocle3_anno_PBS.pdf",height = 4, width = 4)
plot_cells(cds_PBS,
           group_label_size = 4,
           color_cells_by = "anno",
           group_cells_by = "anno")
dev.off()
rm(cds_PBS)

cds_T2 <- cds[,pData(cds)$Sample%in%"T2"]
pdf("plot/CD4_Monocle3_anno_T2.pdf",height = 4, width = 4)
plot_cells(cds_T2,
           group_label_size = 4,
           color_cells_by = "anno",
           group_cells_by = "anno")
dev.off()
rm(cds_T2)


cds <- cluster_cells(cds, resolution=1e-5)

pdf("plot/CD4_Monocle3_NewCluster.pdf",height = 4, width = 4)
plot_cells(cds,
           group_label_size = 5,
           color_cells_by = "cluster")
dev.off()


cds_PBS <- cds[,pData(cds)$Sample%in%"PBS"]
pdf("plot/CD4_Monocle3_NewCluster_PBS.pdf",height = 4, width = 4)
plot_cells(cds_PBS,
           group_label_size = 5,
           color_cells_by = "cluster")
dev.off()
rm(cds_PBS)

cds_T2 <- cds[,pData(cds)$Sample%in%"T2"]
pdf("plot/CD4_Monocle3_NewCluster_T2.pdf",height = 4, width = 4)
plot_cells(cds_T2,
           group_label_size = 5,
           color_cells_by = "cluster")
dev.off()
rm(cds_T2)



cds_umap <- reducedDims(cds)$UMAP
colnames(cds_umap) <- c("umap_1","umap_2")
class(cds_umap)
seurat_cd4 <- seurat[,seurat$anno%in%unique(seurat@meta.data$anno)[1:2]]
seurat_cd4 <- NormalizeData(seurat_cd4)
seurat_cd4 <- FindVariableFeatures(seurat_cd4)
seurat_cd4 <- ScaleData(seurat_cd4)
seurat_cd4 <- RunPCA(seurat_cd4)
seurat_cd4 <- FindNeighbors(seurat_cd4,reduction = "pca", dims = 1:10)
seurat_cd4 <- RunUMAP(seurat_cd4,reduction = "pca",dims = 1:10)
seurat_cd4@meta.data$new_anno <- cds@clusters@listData[["UMAP"]][["cluster_result"]][["optim_res"]][["membership"]]

seurat_cd4@reductions$umap@cell.embeddings <- cds_umap



# For Extended Data Fig.3c in Article
pdf("plot/CD4_Monocle3_anno_New.pdf",height = 7, width = 7)
DimPlot(seurat_cd4,reduction = "umap",group.by = "anno",pt.size = 1)+
  ggtitle("")+
  theme(legend.position = "right")+
  guides(colour=guide_legend(override.aes = list(size =6),ncol = 2)) +
  theme(panel.grid = element_blank())+
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_blank()) +
  theme(
    text = element_text(size = 20),
    legend.title = element_blank(), 
    legend.text = element_text(size = 14,hjust = 0))+
  NoLegend()
dev.off()


pdf("plot/CD4_Monocle3_anno_NewCluster_New.pdf",height = 7, width = 7)
DimPlot(seurat_cd4,reduction = "umap",group.by = "new_anno",pt.size = 1)+
  ggtitle("")+
  theme(legend.position = "right")+
  guides(colour=guide_legend(override.aes = list(size =6),ncol = 2)) +
  theme(panel.grid = element_blank())+
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_blank()) +
  theme(
    text = element_text(size = 20),
    legend.title = element_blank(), 
    legend.text = element_text(size = 14,hjust = 0))+
  NoLegend()
dev.off()


# For Extended Data Fig.3d in Article
pdf("plot/CD4_Monocle3_anno_NewCluster_DotPlot_New.pdf",height = 5, width = 4)
DotPlot(seurat_cd4,features = rev(c("Cd3e","Cd4","Foxp3","Ifng")),group.by = "new_anno",scale = F,
        cols = c("white", "#CD3333"))+
  xlab("")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(fill = NA,colour = "black",size = 1,linetype = "solid")) +
  theme(legend.position = "right",
        axis.title.x.bottom = element_blank(),
        axis.text.x = element_text(angle = 45,size = 12,hjust = 1,colour = "black"),
        axis.text.y = element_text(size = 12,colour = "black"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+coord_flip()
dev.off()



treg_cell_cluster <- read.table("data/treg_cell_cluster.txt",header = T,quote = "",sep = "\t")
seurat_treg <- seurat_cd4[,colnames(seurat_cd4)%in%treg_cell_cluster$cell_id]
seurat_treg@meta.data$new_anno_2 <- treg_cell_cluster$highRes
dim(treg_cell_cluster)
dim(seurat_treg)


# For Fig.3d in Article
pdf("plot/Treg_Monocle3_NewCluster_New.pdf",height = 7, width = 7)
DimPlot(seurat_treg,reduction = "umap",group.by = "new_anno_2",pt.size = 3)+
  ggtitle("")+
  theme(legend.position = "right")+
  guides(colour=guide_legend(override.aes = list(size =6),ncol = 2)) +
  theme(panel.grid = element_blank())+
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_blank()) +
  theme(
    text = element_text(size = 20),
    legend.title = element_blank(), 
    legend.text = element_text(size = 14,hjust = 0))+
  NoLegend()
dev.off()

pdf("plot/Treg_Monocle3_NewCluster_PBS_New.pdf",height = 7, width = 7)
DimPlot(seurat_treg[,seurat_treg$Sample == "PBS"],reduction = "umap",group.by = "new_anno_2",pt.size = 3)+
  ggtitle("")+
  theme(legend.position = "right")+
  guides(colour=guide_legend(override.aes = list(size =6),ncol = 2)) +
  theme(panel.grid = element_blank())+
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_blank()) +
  theme(
    text = element_text(size = 20),
    legend.title = element_blank(), 
    legend.text = element_text(size = 14,hjust = 0))+
  NoLegend()
dev.off()

pdf("plot/Treg_Monocle3_NewCluster_T2_New.pdf",height = 7, width = 7)
DimPlot(seurat_treg[,seurat_treg$Sample == "T2"],reduction = "umap",group.by = "new_anno_2",pt.size = 3)+
  ggtitle("")+
  theme(legend.position = "right")+
  guides(colour=guide_legend(override.aes = list(size =6),ncol = 2)) +
  theme(panel.grid = element_blank())+
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_blank()) +
  theme(
    text = element_text(size = 20),
    legend.title = element_blank(), 
    legend.text = element_text(size = 14,hjust = 0))+
  NoLegend()
dev.off()


seurat_treg_DEGs <- FindMarkers(seurat_treg,ident.1 = "T2",ident.2 = "PBS",logfc.threshold = 0.01)
seurat_treg_DEGs$gene <- rownames(seurat_treg_DEGs)
pdf(file = "plot/Treg_Monocle3_NewCluster_Ikzf2_Vln.pdf", width = 4, height = 4)
VlnPlot(seurat_treg,features = "Ikzf2",pt.size = 0,cols = c("grey","#FF9999"),)+
  xlab(seurat_treg_DEGs$p_val[seurat_treg_DEGs$gene=="Ikzf2"])+  
  ylab(seurat_treg_DEGs$avg_log2FC[seurat_treg_DEGs$gene=="Ikzf2"])+
  NoLegend()
dev.off()
rm(seurat_treg_DEGs)



treg_cell_id <- names(cds@clusters@listData[["UMAP"]][["cluster_result"]][["optim_res"]][["membership"]][cds@clusters@listData[["UMAP"]][["cluster_result"]][["optim_res"]][["membership"]] ==2])
write.table(as.data.frame(treg_cell_id), file = "data/treg_cell_id.txt",col.names = F,row.names = F,quote = F,sep = "\t")


