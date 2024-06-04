library(Seurat)
library(ggplot2)
library(tibble)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggcorrplot)
library(RColorBrewer)
library(clusterProfiler)
library(enrichplot)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(enrichR)
library(gprofiler2)
library(STutility)
library(hdf5r)
library(writexl)
library(readxl)
set.seed(12345)


# load final object of Hudacove et al
experiment <- readRDS("T:/Erika/Single cell/Integrated/seurat objects/seurat_mesenchymal_final.rds")

## subsetting cluster with DC from other mesenchymal clusters
Idents(experiment) <- experiment@meta.data$seurat_clusters
experiment <- subset(experiment, idents = c(13), #where cluster 13 has been identified as dermal cells
                     invert = F)

## split the samples by condition
experiment.list <- SplitObject(experiment , split.by = "Condition")


## SCTransform normalization regrressed by mito genes, ribo genes, and cell cycle genes
for (i in 1:length(experiment.list)) {
  experiment.list[[i]] <- SCTransform(experiment.list[[i]], verbose = FALSE, 
                                      vars.to.regress = c("percent.mt", "G2M.Score", "S.Score", "percent.rib"))
}



# select the features for downstream integration 
experiment.features <- SelectIntegrationFeatures(object.list = experiment.list, nfeatures = 3000)
experiment <- PrepSCTIntegration(object.list = experiment.list, anchor.features = experiment.features, 
                                 verbose = TRUE)

experiment.anchors <- FindIntegrationAnchors(object.list = experiment, normalization.method = "SCT", 
                                             anchor.features = experiment.features, verbose = TRUE, dims = 1:30)
experiment.integrated <- IntegrateData(anchorset = experiment.anchors, normalization.method = "SCT", 
                                       verbose = TRUE, dims = 1:30)

# set the assay to an integrated analysis 
DefaultAssay(experiment.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
experiment <- RunPCA(experiment.integrated, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
experiment <- RunUMAP(experiment, reduction = "pca", dims = 1:25)
experiment <- FindNeighbors(experiment, reduction = "pca", dims = 1:25)

# after clustree analysis pick resolution
experiment <- FindClusters(experiment, resolution = 0.6)



## re-arrange your factors
experiment@meta.data$Genotype <- factor(experiment@meta.data$Genotype, levels = c("WT", "MUT"))
experiment@meta.data$Condition <- factor(experiment@meta.data$Condition, levels = c("WT_E12", "MUT_E12", "WT_E13", "MUT_E13"))


DefaultAssay(experiment) <- "RNA"

# normalize RNA data and scale them for heatmap and DGE
experiment <- NormalizeData(experiment, verbose = TRUE)
experiment <- ScaleData(experiment, verbose = TRUE)

#find all markers
all_markers <- FindAllMarkers(experiment, logfc.threshold = 0.25, min.pct = 0.20, only.pos = T)
write_xlsx(all_markers,"C:/Users/mehmet.kaplan/OneDrive - Ústav experimentální medicíny AV ČR, v. v. i/Desktop/all markers.xlsx")

# MAKE MODULE SCORES
DefaultAssay(experiment) <- "SCT"
# Fb module score: from https://doi.org/10.1016/j.devcel.2015.06.023
Fb_module_score <-list(c("Cd34", "Col1a1", "Dlk1", "Dpp4", "Dpt", "Ephb3", "Fbn1", "Lox", "Lrig1",
                         "Wnt11", "Lum", "Ccbe1", "Eln", "Sulf1", "Col1a2", "Twist1"))
experiment <- AddModuleScore(experiment, features = Fb_module_score, ctrl = 5, name = "Fb Module Score")

# Pre-DC module score: from https://doi.org/10.1016/j.devcel.2018.11.034
Pre_DC_module_Score <- list(c("Ccnd1", "Cpne8", "Kcnab3", "Kcnd3", "Lama1", "Malat1", "Pappa",
                              "Sdk2", "Slain1", "Slc16a9", "Zic3"))
experiment <- AddModuleScore(experiment, features = Pre_DC_module_Score, ctrl = 5, name = "Pre-DC Module Score")

# DC1 module score: from https://doi.org/10.1016/j.devcel.2018.11.034
DC1_module_Score <- list(c("1500009L16Rik", "Endou","Fam221a","Galnt3", "Grik3",
                           "H2-DMa", "Iqsec3", "Lynx1", "Olfr461", "Sgip1"))
experiment <- AddModuleScore(experiment, features = DC1_module_Score, ctrl = 5,
                             name = "DC1_module_Score")

# DC2 module score: from https://doi.org/10.1016/j.devcel.2018.11.034
DC2_module_Score <- list(c("1700019D03Rik", "2700046A07Rik", "2810055G20Rik","4930426D05Rik", "4933424G05Rik",
                           "Adamts20", "Adamts4", "Adarb2", "Ak3", "Anks1b", "Apod","Bmp3", "Cadps2","Cxcr4",
                           "D16Ertd472e", "Dclk2", "Dll1","Doc2g", "Dusp4","Ebf4", "Efcc1", "Epha3","Fam155a","Frzb",
                           "Gm16845", "Inhba", "Jak2","Kctd1", "Klhl14",  "Mir6912","Msra", "Necab1", "Ntng1",
                           "Pcyt1b","Pdlim3", "Plcb1","Prkg2", "Prl2c2","Prss48","Qpct","Rspo3","Rtn1","Scn3a",
                           "Setbp1", "Slitrk5","Smyd3", "Snrpn","Sobp","Sox2ot","Spock1", "Stac2","Syt6", "Thrb",
                           "Tmem132d", "Tspan7","Unc5c"))
experiment <- AddModuleScore(experiment, features = DC2_module_Score, ctrl = 5,
                             name = "DC2_module_Score")




# Cell adhesion module score
cell_adhesion_module_Score <- list(c("Chd2", "Cdh11", "Vcan", "Ncam1", "Actg1", "Fn1", "Tnc", "Alcam"))
experiment <- AddModuleScore(experiment, features = cell_adhesion_module_Score, ctrl = 5, name = "cell_adhesion_module_Score")

# ECM module score
ECM_module_score <-list(c("Bgn", "Cd34", "Col1a1", "Col1a2", "Col3a1", "Dcn", "Lama4", "P4hb", "Timp2", "Vcam1"))
experiment <- AddModuleScore(experiment, features = ECM_module_score, ctrl = 5, name = "ECM Module Score")


# make condiction cluster
experiment@meta.data$condition_cluster <- paste(experiment@meta.data$Condition ,
                                                experiment@meta.data$seurat_clusters ,
                                                sep = "_")
# make genotype cluster
experiment@meta.data$genotype_cluster <- paste(experiment@meta.data$Genotype,
                                               experiment@meta.data$seurat_clusters ,
                                               sep = "_")
#save your object
saveRDS(experiment, 'mesenchymal_dermalcondansate_cluster_subsetted.Rds')


# DEGs between cluster2 of control and mutant
DefaultAssay(experiment) <- 'RNA'
Idents(experiment) <- experiment@meta.data$genotype_cluster
control_vs_mutaant_DEGs_cluster2 <- FindMarkers(experiment, ident.1 = "WT_2", ident.2 = "MUT_2",
                                                logfc.threshold = 0.25,   min.pct = 0.25)

# write DEGs for cluster 2 markers to excel
names <- rownames(control_vs_mutaant_DEGs_cluster2)
rownames(control_vs_mutaant_DEGs_cluster2) <- NULL
control_vs_mutaant_DEGs_cluster2 <- cbind(names,control_vs_mutaant_DEGs_cluster2)
write_xlsx(control_vs_mutaant_DEGs_cluster2,"C:/Users/mehmet.kaplan/OneDrive - Ústav experimentální medicíny AV ČR, v. v. i/Desktop/control_vs_mutaant_DEGs_cluster2.xlsx")



# DEGs between  E12 control and mutant all clusters
DefaultAssay(experiment) <- 'RNA'
Idents(experiment) <- experiment@meta.data$Condition
control_vs_mutaant_DEGs_E12 <- FindMarkers(experiment, ident.1 = "WT_E12", ident.2 = "MUT_E12",
                                           logfc.threshold = 0.25,   min.pct = 0.25)
# write to excel
names <- rownames(control_vs_mutaant_DEGs_E12)
rownames(control_vs_mutaant_DEGs_E12) <- NULL
control_vs_mutaant_DEGs_E12 <- cbind(names,control_vs_mutaant_DEGs_E12)
write_xlsx(control_vs_mutaant_DEGs_E12,"C:/Users/mehmet.kaplan/OneDrive - Ústav experimentální medicíny AV ČR, v. v. i/Desktop/control_vs_mutaant_DEGs_E12.xlsx")


# DEGs between  E13 control and mutant all clusters
DefaultAssay(experiment) <- 'RNA'
Idents(experiment) <- experiment@meta.data$Condition
control_vs_mutaant_DEGs_E13 <- FindMarkers(experiment, ident.1 = "WT_E13", ident.2 = "MUT_E13",
                                           logfc.threshold = 0.25,   min.pct = 0.25)
# write to excel
names <- rownames(control_vs_mutaant_DEGs_E13)
rownames(control_vs_mutaant_DEGs_E13) <- NULL
control_vs_mutaant_DEGs_E13 <- cbind(names,control_vs_mutaant_DEGs_E13)
write_xlsx(control_vs_mutaant_DEGs_E13,"C:/Users/mehmet.kaplan/OneDrive - Ústav experimentální medicíny AV ČR, v. v. i/Desktop/control_vs_mutaant_DEGs_E13.xlsx")


## clusterProfiler GENE ONTOLOGY
cluster_0_markers <- read_excel("C:/Users/mehmet.kaplan/OneDrive - Ústav experimentální medicíny AV ČR, v. v. i/Desktop/book1.xlsx") #where top 100 markers of cluster 0 is added
# run GO 
GO <- enrichGO(cluster_0_markers$gene,OrgDb = "org.Mm.eg.db",  
               ont = "BP", 
               keyType = "SYMBOL") 
GO_data_frame <- as.data.frame(GO)
write_xlsx(GO_data_frame, "C:/Users/mehmet.kaplan/OneDrive - Ústav experimentální medicíny AV ČR, v. v. i/Desktop/book1.xlsx")











###PLOTS


setwd("C:/Users/mehmet.kaplan/OneDrive - Ústav experimentální medicíny AV ČR, v. v. i/Desktop/whisker paper/scRNA/")



## nonsplited dimplot
p <- DimPlot(experiment, reduction = "umap", label = TRUE, pt.size = 1, label.size = 0)
plot(p)
ggsave("Dimplot non-splitted.tiff", width = 4, height = 4, dpi = 600)


## splited dimplot by condition or genotype
p <- DimPlot(experiment, reduction = "umap", label = TRUE, pt.size = 1, label.size = 0, split.by = "Condition")
p <- DimPlot(experiment, reduction = "umap", label = TRUE, pt.size = 1, label.size = 0, split.by = "Genotype")
plot(p)
ggsave("Dimplot splitted by genotype.tiff", width = 10, height = 4, dpi = 600)



## nonsplited feature and vlnplots
DefaultAssay(experiment) <- 'SCT'
Idents(experiment) <- experiment@meta.data$seurat_clusters

p1 <- FeaturePlot(experiment, features = c("Hck"), min.cutoff = 'q9', max.cutoff = "q90",
                  pt.size = 2, ncol = 1)

p2 <- VlnPlot(experiment, features = c("Hck"))
plot(p1/p2)
ggsave("Zeb2 feature and vlnplot non-splitted.tiff", width = 4, height = 6, dpi = 600)


# HEATMAPS
all_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

DefaultAssay(experiment) <- 'RNA'
p <- DoHeatmap(experiment, features = top10$gene, raster=F) +  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) ) +  guides(color=FALSE)
plot(p)
ggsave("heatmap.tiff", width = 10, height = 7, dpi = 600)

## splited by condition feature and vlnplots
setwd("C:/Users/mehmet.kaplan/OneDrive - Ústav experimentální medicíny AV ČR, v. v. i/Desktop/whisker paper/scRNA/ECM bgn cd34 col1a1 col1a2 col3a1 dcn lama4 p4hb timp2 vcam1")
DefaultAssay(experiment) <- 'SCT'
Idents(experiment) <- experiment@meta.data$Condition

p1 <- FeaturePlot(experiment, features = c("Vcam1"), min.cutoff = 'q9', max.cutoff = "q90",
                  pt.size = 2, split.by = "Condition")

p2 <- VlnPlot(experiment, features = c("Timp2"))
plot(p1/p2)
ggsave("Timp2 feature and vlnplot splitted.tiff", width = 12, height = 8, dpi = 600)


## Fb preDC DC1 DC module scores nonsplitted
DefaultAssay(experiment) <- 'SCT'
Idents(experiment) <- experiment@meta.data$seurat_clusters
p1 <- FeaturePlot(experiment, features = c("Fb Module Score1","Pre-DC Module Score1",
                                           "DC1_module_Score1","DC2_module_Score1"), min.cutoff = 'q9', max.cutoff = "q90",
                  pt.size = 2, ncol = 2)
p2 <- VlnPlot(experiment, features = c("Fb Module Score1","Pre-DC Module Score1",
                                       "DC1_module_Score1","DC2_module_Score1"),
              ncol = 2)
plot(p1/p2)
ggsave("Fb preDC DC1 DC2 module scores feature and vlnplot splitted.tiff", width = 8, height = 16, dpi = 600)


## Cell adhesion module scores splitted
DefaultAssay(experiment) <- 'SCT'
Idents(experiment) <- experiment@meta.data$Condition

p1 <- FeaturePlot(experiment, features = c("cell_adhesion_module_Score1"), min.cutoff = 'q9', max.cutoff = "q90",
                  pt.size = 2, split.by = "Condition")

p2 <- VlnPlot(experiment, features = c("cell_adhesion_module_Score1"))
plot(p1/p2)
ggsave("cell_adhesion_module_Score1 feature and vlnplot splitted.tiff", width = 12, height = 8, dpi = 600)



## ECM adhesion module scores splitted
DefaultAssay(experiment) <- 'SCT'
Idents(experiment) <- experiment@meta.data$Condition

p1 <- FeaturePlot(experiment, features = c("ECM Module Score1"), min.cutoff = 'q9', max.cutoff = "q90",
                  pt.size = 2, split.by = "Condition")

p2 <- VlnPlot(experiment, features = c("ECM Module Score1"))
plot(p1/p2)
ggsave("ECM Module Score1 feature and vlnplot splitted.tiff", width = 12, height = 8, dpi = 600)




head(experiment)

## cluster proportion

## summarizing mouse proportions per cluster
clusters <- as.data.frame(experiment@meta.data) %>%
  group_by(Genotype, seurat_clusters, .drop = FALSE) %>%
  summarize(count = n()) %>%
  group_by(Genotype) %>%
  mutate(Proportion = count/sum(count))


p <- ggplot(clusters, aes(x = Genotype, y = Proportion)) +
  geom_col(aes(fill = seurat_clusters)) +
  labs(fill = "Seurat Clusters", x = "Sample") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(text = element_text(size=20))
plot(p)
ggsave("cluster proportion by genotype.tiff", width = 8, height = 4, dpi = 600)




ECM_module_score <-list(c("Bgn", "Cd34", "Col1a1", "Col1a2", "Col3a1", "Dcn", "Lama4", "P4hb", "Timp2", "Vcam1"))

experiment <- AddModuleScore(experiment, features = ECM_module_score, ctrl = 5, name = "ECM Module Score")
head(experiment)
VlnPlot(experiment, features = c("ECM Module Score1"))



