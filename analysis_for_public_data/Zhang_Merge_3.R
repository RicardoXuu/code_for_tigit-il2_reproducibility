

rm(list = ls())
gc()

library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(openxlsx)
library(pheatmap)
library(reshape2)

load("data/TNBC_Treg_All_SCT.RData")

unique(Treg_All@meta.data$Responder)

#########.
trj_genes <- read.table("data/TNBC_Treg_monocle2_trj_genes_Add_0327.txt", header = F)
trj_genes <- trj_genes$V1


expression_data <- Treg_All@assays$SCT@data
expression_data <- as.data.frame(expression_data)
rownames(expression_data) <- rownames(Treg_All)
colnames(expression_data) <- colnames(Treg_All)

responder_info <- Treg_All$Responder
genes_to_plot <- trj_genes[trj_genes%in%rownames(Treg_All)]

expression_subset <- expression_data[genes_to_plot, ]

expression_average <- as.data.frame(rowMeans(expression_subset[, responder_info == "Responder"])) 
colnames(expression_average) <- "Responder"

expression_average$Non_Responder <- rowMeans(expression_subset[, responder_info == "Non_Responder"])

expression_long <- melt(as.matrix(expression_average), 
                        varnames = c("Gene", "Group"), 
                        value.name = "Expression")


pdf("plot/TNBC_Treg_All_NR_VS_R_trj_genes_exp_0327.pdf", width = 3, height = 6)
ggplot(expression_long, aes(x = Group, y = Gene, fill = Expression)) +
  geom_tile(color = "white", size = 0.5) +  # 添加外框线
  scale_fill_gradient2(low = "white", mid = "pink", high = "red", midpoint = 0.4) + # 颜色映射
  labs(x = NULL, y = "", fill = "Expression") +
  ggtitle("") +
  theme_bw()+
  theme(
    axis.title.y = element_text(size = 20,color = "black",face = "bold"),  # 旋转 x 轴标签
    axis.text.x = element_text(size = 12, face = "bold",color = "black",angle = 45,hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold",color = "black"),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", size = 1.5),
    legend.title = element_text(angle = 90, vjust = 0.5, size = 16,color = "black")  # 旋转 legend 标题
  ) 
dev.off()



pdf("plot/TNBC_Treg_All_NR_VS_R_trj_genes_exp_0327_wide.pdf", width = 6, height = 3)
ggplot(expression_long, aes(x = Group, y = Gene, fill = Expression)) +
  geom_tile(color = "white", size = 0.5) +  # 添加外框线
  scale_fill_gradient2(low = "white", mid = "pink", high = "red", midpoint = 0.4) + # 颜色映射
  labs(x = NULL, y = "", fill = "Expression") +
  ggtitle("") +
  theme_bw()+
  theme(
    axis.title.y = element_text(size = 20,color = "black",face = "bold"),  # 旋转 x 轴标签
    axis.text.x = element_text(size = 12, face = "bold",color = "black",angle = 45,hjust = 1),
    axis.text.y = element_text(size = 12, face = "bold",color = "black"),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    legend.position = "top",
    panel.border = element_rect(color = "black", size = 1.5),
    legend.title = element_text( size = 16,color = "black")  # 旋转 legend 标题
  ) + coord_flip()
dev.off()


rm(list =ls())
gc()


###################. 

