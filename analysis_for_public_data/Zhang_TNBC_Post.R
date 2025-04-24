
library(Seurat)
library(openxlsx)
library(ggplot2)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(scales)  
library(grid)   
library(RColorBrewer)


Treg_Post_counts <- read.table("data/Treg_Post_counts.txt",row.names = 1,sep = "\t",quote = "")
Treg_Post <- CreateSeuratObject(counts = Treg_Post_counts, project = "TNBC_Treg_Post")


Treg_Post_meta <- read.table("data/Treg_Post_meta.txt",row.names = 1,sep = "\t",quote = "")
Treg_Post@meta.data <- cbind(Treg_Post@meta.data,Treg_Post_meta[,5:16])

Treg_Post <- NormalizeData(Treg_Post)


Treg_Post <- SCTransform(Treg_Post)
DefaultAssay(Treg_Post) <- "SCT"
Treg_Post <- NormalizeData(Treg_Post)

Treg_Post <- FindVariableFeatures(Treg_Post)
Treg_Post <- ScaleData(Treg_Post)


Treg_Post <- RunPCA(Treg_Post)
Treg_Post <- FindNeighbors(Treg_Post)
Treg_Post <- RunUMAP(Treg_Post,reduction = "pca",dims = 1:10)
Treg_Post <- RunTSNE(Treg_Post,reduction = "pca",dims = 1:10)


Treg_Post_1 <- Treg_Post
Treg_Post_1@reductions$umap@cell.embeddings[,1] <- (Treg_Post_1@reductions$umap@cell.embeddings[,2])
Treg_Post_1@reductions$umap@cell.embeddings[,2] <- -(Treg_Post@reductions$umap@cell.embeddings[,1])
Treg_Post <- Treg_Post_1
rm(Treg_Post_1)


Idents(Treg_Post) <- Treg_Post$Responder

DEGs_Post <- FindMarkers(Treg_Post,ident.1 = "Non_Responder",ident.2 = "Responder")

DEGs_Post <- DEGs_Post[DEGs_Post$avg_log2FC >= 0.32 & DEGs_Post$p_val <= 0.05, ]
DEGs_Post$gene <- rownames(DEGs_Post)
write.xlsx(DEGs_Post,file = "data/DEGs_Post_Zhang_TNBC.xlsx")



# Go enrich
ego_bp_1 <- enrichGO(gene     = DEGs_Post$gene,
                     OrgDb      = org.Hs.eg.db,
                     keyType    = 'SYMBOL',
                     ont        = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.1,
                     qvalueCutoff = 0.1)

bb <- ego_bp_1@result
write.xlsx(bb,file = "data/GO_DEGs_Post_Zhang_TNBC.xlsx")




DEGs_Post_Filter <- DEGs_Post[!(DEGs_Post$gene%in%DEGs$gene),]
# Go enrich
ego_bp <- enrichGO(gene       = DEGs_Post_Filter$gene,
                   OrgDb      = org.Hs.eg.db,
                   keyType    = 'SYMBOL',
                   ont        = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.1,
                   qvalueCutoff = 0.1)

aa <- ego_bp@result
write.xlsx(aa,file = "data/GO_DEGs_Post_Filter_Zhang_TNBC.xlsx")



bp_used <- read.xlsx("data/GO_Treg_0517.xlsx",colNames = F)
colnames(bp_used)[1] <- "term"

bp_used <- aa[aa$ID%in%bp_used$term,]

colnames(bp_used)



bp_used$log_p_adjust <- -log10(bp_used$p.adjust)

break_points <- c(max(bp_used$p.adjust), min(bp_used$p.adjust)*10)
break_points <- c(0.02,1e-3,1e-4,1e-5)
colors <- c("pink", "red")
log_breaks <- -log10(break_points)

labels <- c("0.02","1e-3","1e-4","1e-5")

pdf(file = "plot/GO_BP_Treg_Post_Barplot.pdf", width = 12, height = 12)
ggplot(bp_used, aes(x = reorder(Description, Count), y = Count, fill = log_p_adjust)) +
  geom_bar(stat = "identity", width = 0.9) +
  scale_fill_gradientn(colors = colors, 
                       values = rescale(c(log_breaks, max(bp_used$log_p_adjust))),
                       breaks = log_breaks,
                       labels = labels,  
                       guide = guide_colourbar(title = "p.adjust", barwidth = 2, barheight = 20)) +
  labs(fill = "p.adjust") +
  xlab("") +
  ylab("") +
  ylim(0, max(bp_used$Count[1:20]) + 1) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
  coord_flip() +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 1.5),
        text = element_text(size = 22),
        axis.title = element_text(size = 22),
        axis.text.y = element_text(size = 22),
  )  
dev.off()




FOXP3_Pre <- FindMarkers(Treg,ident.1 = "Non_Responder",ident.2 = "Responder",logfc.threshold = 0.01)
FOXP3_Pre$gene <- rownames(FOXP3_Pre)
pdf(file = "plot/Treg_FOXP3_Vln.pdf", width = 4, height = 6)
VlnPlot(Treg,features = "FOXP3",pt.size = 0)+
  xlab(FOXP3_Pre$p_val[FOXP3_Pre$gene=="FOXP3"])+
  ylab(FOXP3_Pre$avg_log2FC[FOXP3_Pre$gene=="FOXP3"])+
  NoLegend()
dev.off()



FOXP3_Post <- FindMarkers(Treg_Post,ident.1 = "Non_Responder",ident.2 = "Responder",logfc.threshold = 0.01)
FOXP3_Post$gene <- rownames(FOXP3_Post)
pdf(file = "plot/Treg_Post_FOXP3_Vln.pdf", width = 4, height = 6)
VlnPlot(Treg_Post,features = "FOXP3",pt.size = 0)+
  xlab(FOXP3_Post$p_val[FOXP3_Post$gene=="FOXP3"])+  
  ylab(FOXP3_Post$avg_log2FC[FOXP3_Post$gene=="FOXP3"])+
  NoLegend()
dev.off()







