##====================Seurat============================================================
library(dplyr)
library("Seurat")
#----------Setup the Seurat objects
HeLa_1.data <- Read10X(data.dir = "./CL100121924_L01_HeLa_1/filtered_feature_bc_matrix/") 
SiHa_2.data <- Read10X(data.dir = "./CL100121924_L02_SiHa_2/filtered_feature_bc_matrix")
HeLa_1.data@Dimnames[[2]] <- paste("HeLa_1_", c(1:length(HeLa_1.data@Dimnames[[2]])), sep = "")
SiHa_2.data@Dimnames[[2]] <- paste("SiHa_2_", c(1:length(SiHa_2.data@Dimnames[[2]])), sep = "")
HeLa_1 <- CreateSeuratObject(counts = HeLa_1.data, project = "HeLa_1", min.cells = 3, min.features = 200)
SiHa_2 <- CreateSeuratObject(counts = SiHa_2.data, project = "SiHa_2", min.cells = 3, min.features = 200)

#----------Standard pre-processing workflow
HeLa_1[["percent.mt"]] <- PercentageFeatureSet(object = HeLa_1, pattern = "^MT-")
SiHa_2[["percent.mt"]] <- PercentageFeatureSet(object = SiHa_2, pattern = "^MT-")

HeLa_1 <- subset(x = HeLa_1, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 15)
SiHa_2 <- subset(x = SiHa_2, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 15)

#----------Normalizing the data
HeLa_1 <- NormalizeData(object = HeLa_1, normalization.method = "LogNormalize", scale.factor = 10000)
SiHa_2 <- NormalizeData(object = SiHa_2, normalization.method = "LogNormalize", scale.factor = 10000)

#----------Identification of highly variable features(feature selection)
HeLa_1 <- FindVariableFeatures(object = HeLa_1, selection.method = "vst", nfeatures = 2000)
SiHa_2 <- FindVariableFeatures(object = SiHa_2, selection.method = "vst", nfeatures = 2000)

#----------Perform integration
cervical.anchors <- FindIntegrationAnchors(object.list = list(HeLa_1, SiHa_2), dims = 1:20)
cervical.combined <- IntegrateData(anchorset = cervical.anchors, dims = 1:20)
#----------Perform an integrated analysis
DefaultAssay(object = cervical.combined) <- "integrated"
cervical.combined <- ScaleData(object = cervical.combined, verbose = F)
#----------Visualization the integrated data
library(rcolorbrewer)
plots <- VlnPlot(cervical.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident", pt.size = 0, combine = F, cols = brewer.pal(8, "Set1"))
CombinePlots(plots = plots, ncol = 3)

plot1 <- FeatureScatter(cervical.combined, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident", cols = brewer.pal(8, "Set1"))
plot2 <- FeatureScatter(cervical.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident", cols = brewer.pal(8, "Set1"))
CombinePlots(plots = list(plot1,plot2))

x <- c("passedcells","passedcells","passedfeatures","passedfeatures")
y <- c("HeLa","SiHa","HeLa","SiHa")
z <- c(7999,12963,24365,24897)
data <- data.frame(x = x,y=y,z=z)
plot1 <- ggplot(data = data, aes(x = factor(x), y = z, fill = y)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  scale_color_manual(name = "", values = brewer.pal(8, "Set1"), limits = c("HeLa", "SiHa")) + 
  xlab('')+ylab('')+theme_bw()+theme(panel.background=element_blank(),panel.grid=element_blank()) 
plot1

#----------Run PCA
cervical.combined <- RunPCA(object = cervical.combined, npcs = 30, verbose = F, features = VariableFeatures(object = cervical.combined))
print(x = cervical.combined[["pca"]], dims = 1:5, nfeatures = 5)

#----------determine the dims of the dataset
cervical.combined <- JackStraw(object = cervical.combined, num.replicate = 100)
cervical.combined <- ScoreJackStraw(object = cervical.combined, dims = 1:20)

JackStrawPlot(object = cervical.combined, dims = 1:20)
ElbowPlot(object = cervical.combined)

#----------Cluster and Visualization
cervical.combined <- FindNeighbors(object = cervical.combined, dims = 1:20)
cervical.combined <- FindClusters(object = cervical.combined, resolution = 0.5)

cervical.combined <- RunUMAP(object = cervical.combined, dims = 1:20, reduction = "pca")
#umap_batch
DimPlot(object = cervical.combined, reduction = "umap", group.by = "orig.ident")
#umap
DimPlot(object = cervical.combined, reduction = "umap", label = T)
#umap split by SiHa and HeLa
DimPlot(object = cervical.combined, reduction = "umap", split.by = "orig.ident", label = T)


#-----------Identify conserved cell type markers
DefaultAssay(object = cervical.combined) <- "RNA"
cervical.markers <- FindAllMarkers(object = cervical.combined, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)

#Visualization pan-epi markers expression
FeaturePlot(cervical.combined, features = c("CD24","CD44","TJP1","UBC","KRT7"), cols = c("grey","red"))

#-----------find diff exp genes
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
for(i in c(0:12)){
name <-i
cluster_i <- subset(cervical.combined, idents = i)
Idents(cluster_i) <- "orig.ident"
cluster_i_diffexp.genes <- FindMarkers(cluster_i, ident.1 = "HeLa",ident.2 = "SiHa",print.bar = F)
cluster_i_diffexp.genes <- cluster_i_diffexp.genes[cluster_i_diffexp.genes$p_val_adj < 0.01,] 
cluster_i_diffexp.genes$change <-  as.factor(ifelse(cluster_i_diffexp.genes$p_val_adj < 0.01 & abs(cluster_i_diffexp.genes$avg_logFC) > 0.5,ifelse(cluster_i_diffexp.genes$avg_logFC > 0.5,'UP','DOWN'),'NOT'))
cluster_i_diffexp.genes$sign <- ifelse(cluster_i_diffexp.genes$p_val_adj < 0.01 & abs(cluster_i_diffexp.genes$avg_logFC) > 2,rownames(cluster_i_diffexp.genes),NA)
p <- ggplot(data = cluster_i_diffexp.genes, aes(x = avg_logFC, y = -log10(p_val_adj), color = change)) +
  geom_point(alpha=0.8, size = 1) +
  theme_bw(base_size = 15) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  ylab("-log10P.adj") + 
  ggtitle(paste("cluster", i, "volcano",sep = "_")) +
  scale_color_manual(name = "", values = c("#BC3C28","#0072B5","grey"), limits = c("UP", "DOWN", "NOT"))
p
pdf(paste("cluster", i, "_volcano.pdf",sep = "_"),width = 8)
print(p)
dev.off()
}

#----------Analyze the cell cycle gene
cc.genes <- read.table("./Human_cell_cycle_gene.txt")
DoHeatmap(cervical.combined, features = cc.genes, group.by = "seurat_clusters") + 
  scale_fill_gradientn(colors = c("blue","white","red")) 

#----------Analyze the TF&Rec&Lig
newpalette<-colorRampPalette(brewer.pal(12,"Set3"))(13)
TF <- read.table("TF.txt")
for (i in TF[1:1665,]){
 try(plots <- VlnPlot(cervical.combined,features = i, pt.size = 0, cols = newpalette))
 name_TF <- i
 pdf(paste("TF",name_TF,"Vlnplot.pdf",sep = "_"),width = 14)
 plot(plots)
 dev.off()
}
lig <- read.table("ligend.txt")
for (i in lig[1:690,]){
 try(plots <- VlnPlot(cervical.combined,features = i, pt.size = 0, cols = newpalette))
 name_lig <- i
 pdf(paste("lig",name_lig,"Vlnplot.pdf",sep = "_"),width = 14)
 plot(plots)
 dev.off()
}
rec <- read.table("receptor.txt")
for (i in rec[1:708,]){
 try(plots <- VlnPlot(cervical.combined,features = i, pt.size = 0, cols = newpalette))
 name_rec <- i
 pdf(paste("rec",name_rec,"Vlnplot.pdf",sep = "_"),width = 14)
 plot(plots)
 dev.off()
}
