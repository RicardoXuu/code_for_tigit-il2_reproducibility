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


Gene_list <- read.xlsx("data/DEGs_0325.xlsx",sheet = 2)
capitalize_first <- function(x) {
  sapply(x, function(word) {
    if (!is.na(word)) {
      paste(toupper(substr(word, 1, 1)), substr(word, 2, nchar(word)), sep = "")
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

seurat_neu <- seurat[,seurat@meta.data$anno == "Neutrophils"]


DEGs <- FindMarkers(seurat_neu,ident.1 = "T2",ident.2 = "PBS")

DEGs <- DEGs[DEGs$avg_log2FC >= 0.32 & DEGs$p_val <= 0.05, ]
DEGs$gene <- rownames(DEGs)


#########
# MHC.molecular

DEGs_1 <- DEGs[DEGs$gene%in%Gene_list[["MHC.molecular"]],]

DEGs_1_count <- seurat_neu@assays$RNA@layers$data[rownames(seurat_neu)%in%DEGs_1$gene,]
DEGs_1_count <- as.data.frame(DEGs_1_count)
colnames(DEGs_1_count) <- colnames(seurat_neu)
rownames(DEGs_1_count) <- rownames(seurat_neu)[rownames(seurat_neu)%in%DEGs_1$gene]

table(seurat_neu$Sample)

DEGs_1_sum <- as.data.frame(matrix(0,nrow = nrow(DEGs_1_count),ncol = 2))
rownames(DEGs_1_sum) <- rownames(DEGs_1_count)
colnames(DEGs_1_sum) <- c("T2","PBS")
for (i in 1:nrow(DEGs_1_count)) {
  DEGs_1_sum$T2 <- rowMeans(DEGs_1_count[,colnames(DEGs_1_count)%in%colnames(seurat_neu)[seurat_neu$Sample == "T2"]])
  DEGs_1_sum$PBS <- rowMeans(DEGs_1_count[,colnames(DEGs_1_count)%in%colnames(seurat_neu)[seurat_neu$Sample == "PBS"]])
  
}


pdf(file = "plot/MHC_Neutrophils_Dotplot.pdf", width = 8, height = 4)
dp <- DotPlot(seurat_neu,
              group.by = "Sample",
              #scale = F,
              features = c(Gene_list[["MHC.molecular"]][1:5],"H2-Q4","H2-Q6","H2-Q7","H2-T22") ) +
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

pdf(file = "plot/MHC_Neutrophils_Barplot.pdf", width = 4, height = 6)
print(p)
dev.off()









########







########




########
# GO term 
go_term <- "GO:0034341"
GO_1_list <- mouse_gene2go$ensembl[mouse_gene2go$go_term %in% go_term]
GO_1_list <- gene2ensembl[gene2ensembl$ensembl%in%GO_1_list,]


DEGs_2 <- DEGs[DEGs$gene%in%GO_1_list$symbol,]

DEGs_2_count <- seurat_neu@assays$RNA@layers$data[rownames(seurat_neu)%in%DEGs_2$gene,]
DEGs_2_count <- as.data.frame(DEGs_2_count)
colnames(DEGs_2_count) <- colnames(seurat_neu)
rownames(DEGs_2_count) <- rownames(seurat_neu)[rownames(seurat_neu)%in%DEGs_2$gene]

table(seurat_neu$Sample)

DEGs_2_sum <- as.data.frame(matrix(0,nrow = nrow(DEGs_2_count),ncol = 4))
rownames(DEGs_2_sum) <- rownames(DEGs_2_count)
colnames(DEGs_2_sum) <- c("T2","PBS","All","Change")
for (i in 1:nrow(DEGs_2_count)) {
  DEGs_2_sum$T2 <- rowMeans(DEGs_2_count[,colnames(DEGs_2_count)%in%colnames(seurat_neu)[seurat_neu$Sample == "T2"]])
  DEGs_2_sum$PBS <- rowMeans(DEGs_2_count[,colnames(DEGs_2_count)%in%colnames(seurat_neu)[seurat_neu$Sample == "PBS"]])
  DEGs_2_sum$All <- rowMeans(DEGs_2_count)
  DEGs_2_sum$Change <- DEGs_2_sum$T2 / DEGs_2_sum$PBS
  
}

p1 <- pheatmap(DEGs_2_sum[,1:2],
         scale = "row")

########



########
# GO list
go_list <- read.xlsx("data/DEGs_0325.xlsx",sheet = 3,colNames = F)
colnames(go_list) <- "term"
go_list_add <- as.data.frame(c("GO:0034340","GO:0034341"))
colnames(go_list_add) <- "term"
go_list  <- rbind(go_list,go_list_add)
rm(go_list_add)
for (i in 1:nrow(go_list)) {
  go_list$num[i] <- nrow(gene2ensembl[gene2ensembl$ensembl%in%(mouse_gene2go$ensembl[mouse_gene2go$go_term %in% go_list$term[i]]),])
}
go_list <- go_list[-2,]


GO_Matrix <- as.data.frame(matrix(0,nrow = nrow(go_list),ncol = ncol(seurat)))
rownames(GO_Matrix) <- gsub(":","_",go_list$term)
colnames(GO_Matrix) <- colnames(seurat)
for (i in 1:nrow(go_list)) {
  go_term <- go_list$term[i]
  go_gene <- gene2ensembl$symbol[gene2ensembl$ensembl%in%(mouse_gene2go$ensembl[mouse_gene2go$go_term %in% go_term])]
  all_matrix <- as.matrix(seurat@assays$RNA@layers$data)
  go_exp <- all_matrix[rownames(seurat) %in% go_gene]
  GO_Matrix[i,] <- colMeans(seurat@assays$RNA@layers$data[rownames(seurat) %in% go_gene,])
}

go_term


GO_matrix <- as.matrix(GO_Matrix,"dgCMatrix")
dim(GO_matrix)

assay.v5 <- CreateAssay5Object(counts = GO_matrix)
seurat_go <- CreateSeuratObject(assay.v5)
seurat_go@meta.data <- cbind(seurat_go@meta.data,seurat@meta.data[,4:6])


VlnPlot(seurat_go[,seurat_go$anno == "Neutrophils"],
        features = gsub(":","-",go_list$term),
        group.by = "Sample")




# Go enrich
DEGs$gene[DEGs$gene == "Gbp6.1"] <- "Gbp6"
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

bp_used <- read.xlsx("data/DEGs_0325.xlsx",sheet = 4,colNames = F)
colnames(bp_used) <- "term"

bp_used$term <- sprintf("%07d", as.numeric(bp_used$term))
bp_used$term <- str_c("GO:",bp_used$term)
bp_used <- aa[aa$ID%in%bp_used$term,]

colnames(bp_used)



library(ggplot2)
library(scales)  
library(grid)   


bp_used$log_p_adjust <- -log10(bp_used$p.adjust)

break_points <- c(max(bp_used$p.adjust), 0.05, min(bp_used$p.adjust)*10)
colors <- c("skyblue", "pink", "red")
log_breaks <- -log10(break_points)

labels <- c("0.55", "0.05", "5.30e-8")


pdf(file = "plot/GO_BP_Neutrophils_Barplot.pdf", width = 12, height = 10)
ggplot(bp_used, aes(x = reorder(Description, Count), y = Count, fill = log_p_adjust)) +
  geom_bar(stat = "identity", width = 0.9) +
  scale_fill_gradientn(colors = colors, 
                       values = rescale(c(min(bp_used$log_p_adjust), log_breaks, max(bp_used$log_p_adjust))),
                       breaks = log_breaks,
                       labels = labels,  
                       guide = guide_colourbar(title = "p.adjust", barwidth = 2, barheight = 20)) +
  labs(fill = "p.adjust") +
  xlab("") +
  ylab("") +
  ylim(0, max(bp_used$Count[1:20]) + 1) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 50)) +
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



# Gene Set Enrichment Analysis（GSEA）
genelist <- DEGs$avg_log2FC
names(genelist) <- DEGs$gene
head(genelist)
genelist <- sort(genelist, decreasing = TRUE)
head(genelist)
gsemf <- gseGO(genelist,
               OrgDb = org.Mm.eg.db,
               keyType = "SYMBOL",
               ont="BP"
)
head(gsemf)
str(gsemf)
gseaplot(gsemf, geneSetID="GO:0034341")












