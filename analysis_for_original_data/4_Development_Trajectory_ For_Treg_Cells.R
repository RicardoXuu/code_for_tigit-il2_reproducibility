

# Code File 4: Show Development Trajectory of Treg cells based on Monocle2

library(Seurat)
library(monocle)

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
gc()

seurat_cd4 <- seurat[,seurat$anno%in%unique(seurat@meta.data$anno)[1:2]]


treg_cell_id <- readLines("data/treg_cell_id.txt")

seurat_treg <- seurat_cd4[,colnames(seurat_cd4)%in%treg_cell_id]

 
treg_cell_cluster <- read.table("data/treg_cell_cluster.txt",header = T,quote = "",sep = "\t")
seurat_treg@meta.data$highRes <- as.character(treg_cell_cluster$highRes)

expr_matrix <- GetAssayData(seurat_treg, assay = "RNA", slot = "counts")
expr_matrix <- as(expr_matrix, "sparseMatrix")
cell_metadata <- seurat_treg@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expr_matrix), row.names = rownames(expr_matrix))

cds_new <- newCellDataSet(expr_matrix,
                          phenoData = new("AnnotatedDataFrame", data = cell_metadata),
                          featureData = new("AnnotatedDataFrame", data = gene_annotation))
rm(expr_matrix,cell_metadata,gene_annotation)

cds_new <- detectGenes(cds_new,min_expr = 0.1)  

cds_new <- estimateSizeFactors(cds_new)
cds_new <- estimateDispersions(cds_new,parallel = T)

cds_new_disptable <- dispersionTable(cds_new)
cds_new_reDim <- setOrderingFilter(cds_new,
                                   cds_new_disptable$gene_id)

rm(cds_new_disptable)

cds_new_reDim <- reduceDimension(cds_new_reDim,
                                 max_components = 2,
                                 method = 'DDRTree')

cds_new_reDim <- orderCells(cds_new_reDim, reverse = TRUE)

cds_new_reDim@reducedDimS[1,] <- -cds_new_reDim@reducedDimS[1,] 
cds_new_reDim@reducedDimK[1,] <- -cds_new_reDim@reducedDimK[1,]


pdf("plot/Treg_Monocle2_Sample.pdf",height = 4)
plot_cell_trajectory(cds_new_reDim, 
                     cell_size = 1,
                     show_branch_points =F,
                     color_by = "Sample")
dev.off()

# For Figure 3g in Article
pdf("plot/Treg_Monocle2_Sample_split.pdf",height = 6)
plot_cell_trajectory(cds_new_reDim, 
                     cell_size = 1,
                     show_branch_points =F,
                     color_by = "Sample")+
  facet_wrap(~Sample,nrow = 2)
dev.off()


pdf("plot/Treg_Monocle2_State.pdf",height = 4)
plot_cell_trajectory(cds_new_reDim,
                     cell_size = 1,
                     show_branch_points =F,
                     color_by = "State")
dev.off()


# For Figure 3g in Article
pdf("plot/Treg_Monocle2_Pseudotime.pdf",height = 4)
plot_cell_trajectory(cds_new_reDim, 
                     cell_size = 1,
                     show_state_number =F,
                     show_branch_points =F,
                     color_by = "Pseudotime")
dev.off()

pdf("plot/Treg_Monocle2_anno.pdf",height = 4)
plot_cell_trajectory(cds_new_reDim,
                     cell_size = 1,
                     show_state_number =F,
                     show_branch_points =F,
                     color_by = "anno")
dev.off()

pdf("plot/Treg_Monocle2_anno_split.pdf",height = 6)
plot_cell_trajectory(cds_new_reDim, 
                     cell_size = 1,
                     show_branch_points =F,
                     color_by = "anno")+
  facet_wrap(~anno,nrow = 2)
dev.off()


# For Figure 3g in Article
pdf("plot/Treg_Monocle2_NewCluster.pdf",height = 4)
plot_cell_trajectory(cds_new_reDim,
                     cell_size = 1,
                     show_branch_points =F,
                     color_by = "cluster_higRes")
dev.off()


pdf("plot/Treg_Monocle2_NewCluster_split.pdf",height = 10)
plot_cell_trajectory(cds_new_reDim,
                     cell_size = 1,
                     show_branch_points =F,
                     color_by = "cluster_higRes")+
  facet_wrap(~cluster_higRes,nrow = 6)
dev.off()


save(cds_new_reDim,file = "data/cds_new_reDim_monocle2.RData")


load("data/cds_new_reDim_monocle2.RData")


function_score <- read.xlsx("data/DEGs_0325.xlsx",sheet = 8,colNames = F)
colnames(function_score) <- c("gene","cell_type")

gene_list <- c(function_score$gene)
gene_anno <- as.data.frame(function_score$cell_type)
rownames(gene_anno) <- gene_list
colnames(gene_anno) <- c("Cell Type")

names(gene_list) <- function_score$cell_type
plot_cell_trajectory(cds_new_reDim, 
                     cell_size = 1,
                     show_state_number =F,
                     show_branch_points =F,
                     color_by = "Pseudotime")

diff_test_res <- differentialGeneTest(cds_new_reDim[Treg_function_score$gene,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))


# For Figure 3h in Article
pdf("plot/Function_Score_Treg_New_Heatmap.pdf",width = 4,height = 3)
plot_pseudotime_heatmap(cds_new_reDim[gene_list,],
                        num_clusters = 0,
                        cores = 1,
                        cluster_rows = FALSE,
                        # border_color = NA,
                        # add_annotation_row = gene_anno,
                        show_rownames = T)
dev.off()





pdf("plot/Function_Score_Treg_New_Exp.pdf",width = 7,height = 30)
plot_genes_in_pseudotime(cds_new_reDim[gene_list,])
dev.off()



library(reshape2)
library(dplyr)
Cell_Dist <-melt(table(pData(cds_new_reDim)[,c(4,21)]))


Cell_Dist <- Cell_Dist %>%
  group_by(cluster_higRes) %>%
  mutate(total = sum(value),  
         percent = value / total * 100)  



# For Figure 3E in Article
pdf("plot/Cell_Dist_Score_Treg_New_Exp.pdf",height = 5,width = 3)
p <- ggplot(Cell_Dist, aes(x = factor(cluster_higRes), y = percent, fill = Sample)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(  "PBS" = "grey", "T2" = "#FF9999")) +
  labs(x = "Cluster High Resolution", y = "Percentage (%)", fill = "Sample") + 
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.key.size = unit(0.8, "cm"), 
    legend.title = element_blank(), 
    legend.text = element_text(size = 12),
    axis.ticks.x = element_blank(), 
    axis.ticks.y = element_line(color = "black"), 
    axis.text.y = element_text(color = "black",size = 15), 
    axis.text.x = element_text(color = "black",size = 15), 
    axis.title = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 2)
    
  ) 
print(p)
dev.off()






