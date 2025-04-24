
# 

library(Seurat)
library(pheatmap)
library(openxlsx)
library(stringr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)

library(ggplot2)
library(scales)  
library(grid)    


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


seurat <- JoinLayers(seurat)
seurat <- NormalizeData(seurat)
seurat_macrophage <- seurat[,seurat$anno%in%c("Macrophages")]

seurat_macrophage_5 <- seurat_macrophage[,seurat_macrophage$clusters==5]
seurat_macrophage_13 <- seurat_macrophage[,seurat_macrophage$clusters==13]

DEGs_macrophage_5 <- FindMarkers(seurat_macrophage_5,ident.1 = "T2",ident.2 = "PBS")
DEGs_macrophage_5 <- DEGs_macrophage_5[DEGs_macrophage_5$avg_log2FC >= 1 & DEGs_macrophage_5$p_val <= 0.05, ]

DEGs_macrophage_13 <- FindMarkers(seurat_macrophage_13,ident.1 = "T2",ident.2 = "PBS")
DEGs_macrophage_13 <- DEGs_macrophage_13[DEGs_macrophage_13$avg_log2FC >= 1 & DEGs_macrophage_13$p_val <= 0.05, ]



# DEGs_macrophage_5
DEGs_macrophage_5$gene <- rownames(DEGs_macrophage_5)

DEGs_macrophage_5$ensembl <- "N"
for (i in 1:nrow(DEGs_macrophage_5)) {
  DEGs_macrophage_5$ensembl[i] <- gene2ensembl$ensembl[gene2ensembl$symbol%in%DEGs_macrophage_5$gene[i]]
}

ego_bp_macrophage_5 <- enrichGO(gene        = DEGs_macrophage_5$ensembl,
                                 OrgDb      = org.Mm.eg.db,
                                 keyType    = 'ENSEMBL',
                                 ont        = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.1,
                                 qvalueCutoff  = 0.1)


bp_used <- ego_bp_macrophage_5@result[1:30,]

bp_used$log_p_adjust <- -log10(bp_used$p.adjust)
break_points <- c(0.001, 0.003,  0.009)
labels <- c("1e-3", "3e-3", "9e-3")
log_breaks <- -log10(break_points)
log_breaks <- sort(log_breaks) 




pdf(file = "plot/GO_BP_Macrophage_5_Barplot.pdf", width = 12, height = 18)
ggplot(bp_used, aes(x = reorder(Description, Count), y = Count, fill = log_p_adjust)) +
  geom_bar(stat = "identity", width = 0.9) +
  scale_fill_gradientn(
    colors = colors,
    values = scales::rescale(c(min(bp_used$log_p_adjust), log_breaks, max(bp_used$log_p_adjust))),
    breaks = log_breaks,
    labels = labels, 
    guide = guide_colourbar(title = "p.adjust", barwidth = 2, barheight = 20)
  ) +
  
  labs(fill = "p.adjust") +
  xlab("") +
  ylab("") +
  ylim(0, max(bp_used$Count[1:30]) + 1) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 45)) +
  coord_flip() +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(size = 1.5),
    text = element_text(size = 22),
    axis.title = element_text(size = 22),
    axis.text.y = element_text(size = 22)
  )

dev.off()




# DEGs_macrophage_13
DEGs_macrophage_13$gene <- rownames(DEGs_macrophage_13)
DEGs_macrophage_13 <- DEGs_macrophage_13[DEGs_macrophage_13$gene%in%gene2ensembl$symbol,]
DEGs_macrophage_13$ensembl <- "N"
for (i in 1:nrow(DEGs_macrophage_13)) {
  DEGs_macrophage_13$ensembl[i] <- gene2ensembl$ensembl[gene2ensembl$symbol%in%DEGs_macrophage_13$gene[i]]
}

ego_bp_macrophage_13 <- enrichGO(gene      = DEGs_macrophage_13$ensembl,
                                OrgDb      = org.Mm.eg.db,
                                keyType    = 'ENSEMBL',
                                ont        = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.1,
                                qvalueCutoff  = 0.1)


bp_used <- ego_bp_macrophage_13@result[ego_bp_macrophage_13@result$p.adjust <= 0.05,]

bp_used[,c(2,6)]

bp_used$log_p_adjust <- -log10(bp_used$p.adjust)
break_points <- c(0.0042, 0.02, 0.048)
log_breaks <- -log10(break_points) 
log_breaks <- sort(log_breaks, decreasing = FALSE)
labels <- c("4.2e-3", "2.0e-2", "4.8e-2")



pdf(file = "plot/GO_BP_Macrophage_13_Barplot.pdf", width = 5, height = 2.6)
ggplot(bp_used, aes(x = reorder(Description, Count), y = Count, fill = log_p_adjust)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_gradientn(
    colors = colors,
    values = rescale(c(min(bp_used$log_p_adjust), log_breaks, max(bp_used$log_p_adjust))),
    breaks = log_breaks,
    labels = labels,
    guide = guide_colourbar(title = "p.adjust", barwidth = 2, barheight = 10)
  ) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 30)) +
  labs(fill = "p.adjust") +
  xlab("") +
  ylab("") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    text = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )
dev.off()


