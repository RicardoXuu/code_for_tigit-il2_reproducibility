
rm(list =ls())



library(Seurat)
library(openxlsx)
library(ggplot2)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)


load("data/TNBC_Treg_All_SCT.RData")


#########.
trj_genes <- read.table("data/TNBC_Treg_monocle2_trj_genes.txt", header = F)
trj_genes <- trj_genes$V1


table(rownames(Treg_All) %in% trj_genes) 
trj_genes_exp <- Treg_All@assays$SCT@data[rownames(Treg_All) %in% trj_genes,]
trj_genes_exp_mean <- as.data.frame(colMeans(trj_genes_exp))


Treg_All$trj_genes_exp_mean <- trj_genes_exp_mean$`colMeans(trj_genes_exp)`

new_gene_counts <- matrix(
  trj_genes_exp_mean$`colMeans(trj_genes_exp)`, 
  nrow = 1, 
  dimnames = list("trj_genes_exp_mean", colnames(Treg_All))
)

Treg_All@assays$SCT@counts <- rbind(Treg_All@assays$SCT@counts, new_gene_counts)
Treg_All@assays$SCT@data <- rbind(Treg_All@assays$SCT@data, new_gene_counts)


Idents(Treg_All) <- Treg_All$Responder
unique(Treg_All$Responder)

markers <- FindMarkers(
  object = Treg_All,
  ident.1 = "Non_Responder",
  ident.2 = "Responder",
  features = "trj_genes_exp_mean",
  logfc.threshold = 0,
  min.pct = 0,
  test.use = "wilcox"
)

print(markers)

plot_colors <- c("#F8766D","#00BFC4")
names(plot_colors) <- c("Non_Responder","Responder")
Treg_All$Responder <- factor(Treg_All$Responder, levels = c("Responder","Non_Responder"))

pdf("plot/TNBC_Treg_All_NR_VS_R_trj_genes_exp_mean_vlnplot.pdf")
VlnPlot(Treg_All, features = c("trj_genes_exp_mean"), 
        group.by = "Responder",
        cols = plot_colors,
        pt.size = 0) + NoLegend()
dev.off()


pdf("plot/TNBC_Treg_All_NR_VS_R_trj_genes_exp_mean_plot.pdf",width = 7.7, height = 7)
FeaturePlot(Treg_All, 
            features = "trj_genes_exp_mean", 
            order = TRUE,
            pt.size = 2.5)
dev.off()






#########.
Treg_Responder <- Treg_All[, Treg_All$Responder == "Responder"]

trj_genes <- read.table("data/TNBC_Treg_monocle2_trj_genes.txt", header = F)
trj_genes <- trj_genes$V1


table(rownames(Treg_Responder) %in% trj_genes) 
trj_genes_exp <- Treg_Responder@assays$SCT@data[rownames(Treg_Responder) %in% trj_genes,]
trj_genes_exp_mean <- as.data.frame(colMeans(trj_genes_exp))


Treg_Responder$trj_genes_exp_mean <- trj_genes_exp_mean$`colMeans(trj_genes_exp)`

new_gene_counts <- matrix(
  trj_genes_exp_mean$`colMeans(trj_genes_exp)`, 
  nrow = 1, 
  dimnames = list("trj_genes_exp_mean", colnames(Treg_Responder))
)

Treg_Responder@assays$SCT@counts <- rbind(Treg_Responder@assays$SCT@counts, new_gene_counts)
Treg_Responder@assays$SCT@data <- rbind(Treg_Responder@assays$SCT@data, new_gene_counts)


Idents(Treg_Responder) <- Treg_Responder$SampleSum
unique(Treg_Responder$SampleSum)

markers <- FindMarkers(
  object = Treg_Responder,
  ident.1 = "TNBC_Treg_Post",
  ident.2 = "TNBC_Treg",
  features = "trj_genes_exp_mean",
  logfc.threshold = 0,
  min.pct = 0,
  test.use = "wilcox"
)

print(markers)

plot_colors <- c("#6BAED6","#E69F00")
names(plot_colors) <- c("TNBC_Treg","TNBC_Treg_Post")

pdf("plot/TNBC_Treg_Responder_Post_VS_Pre_trj_genes_exp_mean_vlnplot.pdf")
VlnPlot(Treg_Responder, features = c("trj_genes_exp_mean"), 
        group.by = "SampleSum",
        cols = plot_colors,
        pt.size = 0) + NoLegend()
dev.off()


pdf("plot/TNBC_Treg_Responder_Post_VS_Pre_trj_genes_exp_mean_plot.pdf",width = 7.7, height = 7)
FeaturePlot(Treg_Responder, 
            features = "trj_genes_exp_mean", 
            order = TRUE,
            pt.size = 2.5)
dev.off()







#########.
Treg_Non_Responder <- Treg_All[, Treg_All$Responder == "Non_Responder"]

trj_genes <- read.table("data/TNBC_Treg_monocle2_trj_genes.txt", header = F)
trj_genes <- trj_genes$V1


table(rownames(Treg_Non_Responder) %in% trj_genes) 
trj_genes_exp <- Treg_Non_Responder@assays$SCT@data[rownames(Treg_Non_Responder) %in% trj_genes,]
trj_genes_exp_mean <- as.data.frame(colMeans(trj_genes_exp))


Treg_Non_Responder$trj_genes_exp_mean <- trj_genes_exp_mean$`colMeans(trj_genes_exp)`

new_gene_counts <- matrix(
  trj_genes_exp_mean$`colMeans(trj_genes_exp)`, 
  nrow = 1, 
  dimnames = list("trj_genes_exp_mean", colnames(Treg_Non_Responder))
)

Treg_Non_Responder@assays$SCT@counts <- rbind(Treg_Non_Responder@assays$SCT@counts, new_gene_counts)
Treg_Non_Responder@assays$SCT@data <- rbind(Treg_Non_Responder@assays$SCT@data, new_gene_counts)


Idents(Treg_Non_Responder) <- Treg_Non_Responder$SampleSum
unique(Treg_Non_Responder$SampleSum)

markers <- FindMarkers(
  object = Treg_Non_Responder,
  ident.1 = "TNBC_Treg_Post",
  ident.2 = "TNBC_Treg",
  features = "trj_genes_exp_mean",
  logfc.threshold = 0,
  min.pct = 0,
  test.use = "wilcox"
)

print(markers)

plot_colors <- c("#6BAED6","#E69F00")
names(plot_colors) <- c("TNBC_Treg","TNBC_Treg_Post")

pdf("plot/TNBC_Treg_Non_Responder_Post_VS_Pre_trj_genes_exp_mean_vlnplot.pdf")
VlnPlot(Treg_Non_Responder, features = c("trj_genes_exp_mean"), 
        group.by = "SampleSum",
        cols = plot_colors,
        pt.size = 0) + NoLegend()
dev.off()

pdf("plot/TNBC_Treg_Non_Responder_Post_VS_Pre_trj_genes_exp_mean_plot.pdf",width = 7.7, height = 7)
FeaturePlot(Treg_Non_Responder, 
            features = "trj_genes_exp_mean", 
            order = TRUE,
            pt.size = 2.5)
dev.off()






