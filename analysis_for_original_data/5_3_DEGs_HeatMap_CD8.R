library(Seurat)
library(pheatmap)
library(openxlsx)
library(stringr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)


gene2ensembl <- read.table("data/PBS/features.tsv")
gene2ensembl <- gene2ensembl[,1:2]
colnames(gene2ensembl) <- c("ensembl","symbol")

mouse_gene2go <- read.table("data/mouse_gene2go.txt")
colnames(mouse_gene2go) <- c("ensembl","go_term")


Gene_list <- read.xlsx("data/DEGs_0325.xlsx",sheet = 1)
capitalize_first <- function(x) {
  sapply(x, function(word) {
    if (!is.na(word)) {
      paste(toupper(substr(word, 1, 1)), tolower(substr(word, 2, nchar(word))), sep = "")
    } else {
      NA
    }
  })
}
Gene_list <- lapply(Gene_list, capitalize_first)
for (i in 1:length(Gene_list)) {
  names(Gene_list[[i]]) <- NULL
  Gene_list[[i]] <- Gene_list[[i]][!is.na(Gene_list[[i]])]
}
rm(i)



seurat <- JoinLayers(seurat)
seurat <- NormalizeData(seurat)
unique(seurat$anno)
seurat_cd8 <- seurat[,seurat@meta.data$anno == "CD8+T cells"]


DEGs <- FindMarkers(seurat_cd8,ident.1 = "T2",ident.2 = "PBS")

DEGs <- DEGs[DEGs$avg_log2FC >= 0.32 & DEGs$p_val <= 0.05, ]
DEGs$gene <- rownames(DEGs)


#########
# MHC.molecular

DEGs_1 <- DEGs[DEGs$gene%in%Gene_list[["Cytotoxicity"]],]

DEGs_1_count <- seurat_cd8@assays$RNA@layers$data[rownames(seurat_cd8)%in%DEGs_1$gene,]
DEGs_1_count <- as.data.frame(DEGs_1_count)
colnames(DEGs_1_count) <- colnames(seurat_cd8)
rownames(DEGs_1_count) <- rownames(seurat_cd8)[rownames(seurat_cd8)%in%DEGs_1$gene]

table(seurat_cd8$Sample)

DEGs_1_sum <- as.data.frame(matrix(0,nrow = nrow(DEGs_1_count),ncol = 2))
rownames(DEGs_1_sum) <- rownames(DEGs_1_count)
colnames(DEGs_1_sum) <- c("T2","PBS")
for (i in 1:nrow(DEGs_1_count)) {
  DEGs_1_sum$T2 <- rowMeans(DEGs_1_count[,colnames(DEGs_1_count)%in%colnames(seurat_cd8)[seurat_cd8$Sample == "T2"]])
  DEGs_1_sum$PBS <- rowMeans(DEGs_1_count[,colnames(DEGs_1_count)%in%colnames(seurat_cd8)[seurat_cd8$Sample == "PBS"]])
  
}



pdf(file = "plot/Cytotoxicity_CD8_Dotplot.pdf", width = 8, height = 4)
dp <- DotPlot(seurat_cd8,
              group.by = "Sample",
              #scale = F,
              features = Gene_list[["Cytotoxicity"]] ) +
  xlab("") +
  ylab("")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.border = element_rect(fill = NA,colour = "black",size = 2,linetype = "solid")) +
  theme(legend.position = "top",
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_blank(),
        legend.text = element_text(color = "black",face = "bold"),
        axis.text.x = element_text(angle = 45,size = 12,hjust = 1,color = "black",face ="bold"),
        axis.text.y = element_text(size = 12,color = "black",face = "bold"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
  )
max_value <- max(dp$data$avg.exp.scaled, na.rm = TRUE)
min_value <- min(dp$data$avg.exp.scaled, na.rm = TRUE)

colors <- c("skyblue", "white", "red")
dp <- dp + scale_color_gradientn(colors = colors, 
                                 breaks = c(min_value, 0,max_value), 
                                 labels = c(round(min_value,2),0, round(max_value,2)))

print(dp)
dev.off()




##########
library(RColorBrewer)

dp_scale <- dp[["data"]]
rownames(dp_scale) <- 1:nrow(dp_scale)
for (i in 1:(nrow(dp_scale)/2)) {
  dp_scale$avg.exp.new[i] <- dp_scale$avg.exp[i] / ( (dp_scale$avg.exp[i] + dp_scale$avg.exp[i+nrow(dp_scale)/2]) /2 )
}
for (i in (nrow(dp_scale)/2 +1 ):(nrow(dp_scale))) {
  dp_scale$avg.exp.new[i] <- dp_scale$avg.exp[i] / ( (dp_scale$avg.exp[i] + dp_scale$avg.exp[i-nrow(dp_scale)/2]) /2 )
}



########
# For Figure6h in Article
dp_scale$pct.exp.signed <- ifelse(dp_scale$id == "T2", dp_scale$pct.exp, -dp_scale$pct.exp)

x_limits <- max(abs(dp_scale$pct.exp.signed))
dp_scale$features.plot <- factor(dp_scale$features.plot,
                                 levels = as.character(dp_scale$features.plot[dp_scale$id == "T2"][order(dp_scale$avg.exp.new[dp_scale$id == "T2"],decreasing = F)]))

p <- ggplot(dp_scale, aes(x = features.plot, y = pct.exp.signed, fill = avg.exp.new)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  scale_fill_gradient2(low = "skyblue", mid = "white", high = "red", 
                       midpoint = 1, limits = c(min(dp_scale$avg.exp.new), max(dp_scale$avg.exp.new)), space = "Lab", 
                       name = "Exp Level") +
  labs(x = "", y = "") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(panel.background = element_rect(fill = "#F3F3F3")) +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 2, linetype = "solid")) +
  theme(legend.position = "top",
        axis.title.x.bottom = element_blank(),
        axis.title.y.left = element_blank(),
        legend.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(color = "black", face ="bold"),
        axis.text.y = element_text(size = 12, color = "black", face = "bold"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)
  ) +
  coord_flip() +
  scale_y_continuous(labels = function(x) paste0(abs(x), "%"), limits = c(-x_limits, x_limits))  

pdf(file = "plot/Cytotoxicity_CD8_Barplot.pdf", width = 4, height = 6)
print(p)
dev.off()






########
# GO enrich
DEGs$ensembl <- "N"
for (i in 1:nrow(DEGs)) {
  DEGs$ensembl[i] <- gene2ensembl$ensembl[gene2ensembl$symbol%in%DEGs$gene[i]]
}

ego_bp <- enrichGO(gene       = DEGs$ensembl,
                   OrgDb      = org.Mm.eg.db,
                   keyType    = 'ENSEMBL',
                   ont        = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.1,
                   qvalueCutoff = 0.1)

aa <- ego_bp@result
write.xlsx(aa,file = "data/GO_CD8.xlsx")



bp_used <- read.xlsx("data/DEGs_0325.xlsx",sheet = 5,colNames = F)
colnames(bp_used) <- "term"

bp_used <- aa[aa$ID%in%bp_used$term,]

colnames(bp_used)



library(ggplot2)
library(scales)  
library(grid)   


bp_used$log_p_adjust <- -log10(bp_used$p.adjust)

break_points <- c(max(bp_used$p.adjust), min(bp_used$p.adjust)*10)

break_points
break_points <- c(0.009,1e-3,1e-4,1e-5)
colors <- c("pink", "red")
log_breaks <- -log10(break_points)

labels <- c("0.009","1e-3","1e-4","1e-5")

pdf(file = "plot/GO_BP_CD8_Barplot.pdf", width = 11.5, height = 10)
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












