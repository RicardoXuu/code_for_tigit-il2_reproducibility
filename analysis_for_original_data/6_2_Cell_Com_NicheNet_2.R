
library(Seurat)
library(CellChat)
library(ggplot2)
library(stringr)
library(dplyr)


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


treg_cell_id <- read.table("data/treg_cell_id.txt",header = F,quote = "",sep = "\t")
treg_cell_id <- treg_cell_id$V1


seurat@meta.data$anno_new <- seurat@meta.data$anno


anno_new <- seurat@meta.data[,c(6:7)]
anno_new$cell_id <- rownames(anno_new)
anno_new$anno_new[anno_new$anno%in%c("CD4+Tconv cells","Treg cells")] <- c("CD4+Tconv cells")
anno_new$anno_new[anno_new$cell_id%in%treg_cell_id] <- c("Treg cells")

seurat@meta.data$anno_new <- anno_new$anno_new
table(seurat@meta.data$anno_new)
rm(anno_new,treg_cell_id)

seurat@meta.data$anno_sample <- str_c(seurat@meta.data$anno_new,"_",seurat@meta.data$Sample)
table(seurat@meta.data$anno_sample)


seurat <- JoinLayers(seurat)
seurat <- NormalizeData(seurat)

seurat@meta.data$samples <- seurat@meta.data$Sample


library(circlize)
library(nichenetr)



seurat_part <- seurat[,seurat@meta.data$anno_new%in%c("CD8+T cells","CD4+Tconv cells","Treg cells","Neutrophils")]

lr_network <- readRDS(url("https://zenodo.org/records/7074291/files/gr_network_mouse_21122021.rds"))
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/records/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))


# indicated cell types should be cell class identities
Idents(seurat_part) <- seurat_part@meta.data$anno_new
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = seurat_part, 
  receiver = "Neutrophils", 
  condition_colname = "Sample", 
  condition_oi = "T2", 
  condition_reference = "PBS", 
  sender = c("Treg cells"), 
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)



active_ligand_receptor_links_df <- nichenet_output$ligand_receptor_df
active_ligand_receptor_links_df <- active_ligand_receptor_links_df[order(active_ligand_receptor_links_df$weight,decreasing = T),]
active_ligand_receptor_links_df <- active_ligand_receptor_links_df[order(active_ligand_receptor_links_df$ligand),]
active_ligand_receptor_max <- active_ligand_receptor_links_df[!duplicated(active_ligand_receptor_links_df$ligand),]
active_ligand_receptor_max <- active_ligand_receptor_max[order(active_ligand_receptor_max$weight,decreasing = T),]

active_ligand_receptor_num <- table(active_ligand_receptor_links_df$ligand)[table(active_ligand_receptor_links_df$ligand) >= 0]

active_ligand_receptor_links_df <- active_ligand_receptor_links_df[active_ligand_receptor_links_df$ligand%in%active_ligand_receptor_max$ligand[active_ligand_receptor_max$weight >= 0.5] 
                                                                   & active_ligand_receptor_links_df$ligand%in%names(active_ligand_receptor_num) ,] 
unique(active_ligand_receptor_links_df$ligand)
unique(active_ligand_receptor_links_df$receptor)
rm(active_ligand_receptor_max,active_ligand_receptor_num)



active_ligand_target_links_df <- nichenet_output$ligand_target_df[nichenet_output$ligand_target_df$ligand%in%active_ligand_receptor_links_df$ligand,]
table(active_ligand_target_links_df$target)



#active_ligand_target_links_df = nichenet_output$ligand_target_df %>% mutate(target_type = "Neutrophils") %>% mutate(ligand_type = "Treg cells")

unique(active_ligand_receptor_links_df$receptor)
unique(active_ligand_target_links_df$target)

unique(active_ligand_receptor_links_df$receptor) %in% unique(active_ligand_target_links_df$target)




# 1 receptor gene
circos_links = active_ligand_receptor_links_df


circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% mutate(color_ligand_type = "orange") %>% mutate(color_target_type = "skyblue")
links_circle = circos_links %>% select(ligand,receptor, weight)


ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type
names(grid_ligand_color) <- ligand_color$ligand

target_color = circos_links %>% distinct(receptor,color_target_type)
grid_target_color = target_color$color_target_type
names(grid_target_color) <- target_color$receptor

grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 0.9 - (0.9*weight)) %>% .$transparency 


target_order = circos_links$receptor %>% unique()
ligand_order = circos_links$ligand %>% unique()
order = c(ligand_order,target_order)


width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_target = 15
width_same_cell_same_target_type = 0.5

gaps = c(
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(color_ligand_type == "orange") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(color_target_type == "skyblue") %>% distinct(receptor) %>% nrow() -1)),
  width_different_cell
)





circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = 0, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #
dev.off()

unique(Idents(seurat_part))

circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,
             link.sort = TRUE, link.decreasing = FALSE, 
             grid.col = grid_col,transparency = transparency, 
             diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow", 
             annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) 












# 2 target gene
circos_links = active_ligand_target_links_df


circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% mutate(color_ligand_type = "orange") %>% mutate(color_target_type = "skyblue")
links_circle = circos_links %>% select(ligand,target, weight)


ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type
names(grid_ligand_color) <- ligand_color$ligand

target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type
names(grid_target_color) <- target_color$target

grid_col =c(grid_ligand_color,grid_target_color)

transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 0.5 - (0.5*weight)) %>% .$transparency 


target_order = circos_links$target %>% unique()
ligand_order = circos_links$ligand %>% unique()
order = c(ligand_order,target_order)


width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_target = 15
width_same_cell_same_target_type = 0.5

gaps = c(
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(color_ligand_type == "orange") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(color_target_type == "skyblue") %>% distinct(target) %>% nrow() -1)),
  width_different_cell
)





circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = 0, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #
dev.off()



pdf("plot/NicheNet_Treg_Neutrophils_Target_circos.pdf",height = 9,width = 9)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,
             link.sort = TRUE, link.decreasing = FALSE, 
             grid.col = grid_col,transparency = transparency, 
             diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow", 
             annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) 
dev.off()


