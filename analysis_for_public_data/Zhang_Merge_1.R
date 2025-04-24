

rm(list =ls())


library(Seurat)
library(openxlsx)
library(ggplot2)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)

Treg_counts <- read.table("data/Treg_counts.txt",row.names = 1,sep = "\t",quote = "")
Treg <- CreateSeuratObject(counts = Treg_counts, project = "TNBC_Treg")


Treg_meta <- read.table("data/Treg_meta.txt",row.names = 1,sep = "\t",quote = "")
Treg@meta.data <- cbind(Treg@meta.data,Treg_meta[,5:16])


Treg_Post_counts <- read.table("data/Treg_Post_counts.txt",row.names = 1,sep = "\t",quote = "")
Treg_Post <- CreateSeuratObject(counts = Treg_Post_counts, project = "TNBC_Treg_Post")

Treg_Post_meta <- read.table("data/Treg_Post_meta.txt",row.names = 1,sep = "\t",quote = "")
Treg_Post@meta.data <- cbind(Treg_Post@meta.data,Treg_Post_meta[,5:16])




Treg_All <- merge(x = Treg, y = Treg_Post, project = "TNBC_Treg_All")
Treg_All <- JoinLayers(Treg_All)
Treg_All$SampleSum <- Idents(Treg_All)
table(Treg_All$SampleSum)


Treg_All <- SCTransform(Treg_All)

DefaultAssay(Treg_All) <- "SCT"


Treg_All <- RunPCA(Treg_All)
Treg_All <- FindNeighbors(Treg_All)
Treg_All <- RunUMAP(Treg_All,reduction = "pca",dims = 1:20)
Treg_All <- RunTSNE(Treg_All,reduction = "pca",dims = 1:10)



pdf("plot/TNBC_Treg_All_SCT_UMAP_Patient.pdf",width = 7.7,height = 7)
DimPlot(Treg_All,
        group.by = "Patient",
        pt.size = 2.5,
        reduction = "umap")
dev.off()


pdf("plot/TNBC_Treg_All_SCT_UMAP_Sample.pdf",width = 8.5,height = 7)
DimPlot(Treg_All,
        group.by = "SampleSum",
        cols = c("#6BAED6","#E69F00"),
        pt.size = 2.5,
        reduction = "umap")
dev.off()


pdf("plot/TNBC_Treg_All_SCT_UMAP_Responder.pdf",width = 8.6,height = 7)
DimPlot(Treg_All,
        group.by = "Responder",
        pt.size = 2.5,
        reduction = "umap")
dev.off()



pdf("plot/TNBC_Treg_Post_SCT_UMAP_Patient.pdf",width = 8.6,height = 7)
DimPlot(Treg_All[,Treg_All$SampleSum == "TNBC_Treg_Post"],
        group.by = "Patient",
        pt.size = 2.5,
        reduction = "umap")
dev.off()



Treg_All$Plot <- Treg_All$Responder
Treg_All$Plot[Treg_All$SampleSum == "TNBC_Treg"] <- "Others"
Treg_All$Plot <- factor(Treg_All$Plot,levels = c("Non_Responder","Responder","Others"))



plot_colors <- c("#F8766D","#00BFC4","white")
names(plot_colors) <- c("Non_Responder","Responder","Others")

table(Treg_All$Plot)
pdf("plot/TNBC_Treg_Post_SCT_UMAP_Responder.pdf",width = 7,height = 7)
plot <- DimPlot(Treg_All,
        group.by = "Plot",
        pt.size = 2.5,
        cols = plot_colors,
        order = c("Non_Responder","Responder","Others"),
        reduction = "umap")+NoLegend()
print(plot)
dev.off()




pdf("plot/TNBC_Treg_Pre_SCT_UMAP_Patient.pdf",width = 8.6,height = 7)
DimPlot(Treg_All[,Treg_All$SampleSum == "TNBC_Treg"],
        group.by = "Patient",
        pt.size = 2.5,
        reduction = "umap")
dev.off()



plot_colors <- c("#F8766D","#00BFC4","white")
names(plot_colors) <- c("Non_Responder","Responder","Others")

table(Treg_All$SampleSum)

Treg_All$Plot <- Treg_All$Responder
Treg_All$Plot[Treg_All$SampleSum == "TNBC_Treg_Post"] <- "Others"
Treg_All$Plot <- factor(Treg_All$Plot,levels = c("Non_Responder","Responder","Others"))
table(Treg_All$Plot )

pdf("plot/TNBC_Treg_Pre_SCT_UMAP_Responder.pdf",width = 7,height = 7)
plot <- DimPlot(Treg_All,
                group.by = "Plot",
                pt.size = 2.5,
                cols = plot_colors,
                order = c("Non_Responder","Responder","Others"),
                reduction = "umap")+NoLegend()
print(plot)
dev.off()



Treg_All$Plot  <- "NA"
save(Treg_All,file = "data/TNBC_Treg_All_SCT.RData")
load("data/TNBC_Treg_All_SCT.RData")

sct_counts <- GetAssayData(Treg_All, assay = "SCT", slot = "counts")
sct_counts_mat <- as.matrix(sct_counts)
write.table(sct_counts_mat, file = "data/TNBC_Treg_All_SCT_counts.txt",
            sep = "\t", quote = FALSE, 
            row.names = TRUE, col.names = TRUE)

meta_data <- Treg_All@meta.data
write.table(meta_data, file = "data/TNBC_Treg_All_SCT_meta.txt",
            sep = "\t", quote = FALSE, 
            row.names = TRUE, col.names = TRUE)


