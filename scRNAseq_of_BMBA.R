library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(cowplot)

#Distibution-based background subtration,hashtag/sampletag demultiplexing, and Seurat Object creation for BD Rhapsody-derived WTA data was 
#conducted by using rDBEC package as previously described (Shichino, S. et al. Commun Biol. (2022). doi:10.1038/s42003-022-03536-0  Github: https://github.com/s-shichino1989/TASSeq-paper)

#===============================Quality Control and Selecting Cells===============================#
#perfomed by using R software package Seurat v2.3.4

#load Seurat Object pre-processed by rDBEC package (named as mBC)
load("MiyakeWTA2_Seurat.rda")



#add metadata
colnames(mBC@meta.data)[2]="nReads"

mito.genes = grep(pattern = "^mt.", x = rownames(x = mBC@data), value = TRUE)
percent.mito = Matrix::colSums(mBC@raw.data[mito.genes, ])/Matrix::colSums(mBC@raw.data)
mBC = AddMetaData(object = mBC, metadata = percent.mito, col.name = "percent.mito")

nReads_log = log10(Matrix::colSums(mBC@raw.data))
mBC = AddMetaData(object = mBC, metadata = nReads_log, col.name = "nReads.log")



#filter out doublets and mitochondrial gene-enriched cells
mBC = FilterCells(object = mBC, subset.names = c("percent.mito"),
                  low.thresholds = c(-Inf), high.thresholds = c(0.2))

doublets = mBC@meta.data$TagIDs %in% c("doublet", "not_detected")
mBC = SubsetData(object = mBC, cells.use=rownames(mBC@meta.data[!doublets,]),
                 subset.raw = TRUE)

#===============================Seurat Object Updating, Normalization, Variable Feature Detection===============================#
#perfomed by using R software package Seurat v4.0.4

#update Seurat Object
mBC_new = UpdateSeuratObject(mBC)



#data subsetting 
#We have performed multiplex assay of scRNA-seq analysis for Figure.1&2. We have extracted dataset for Figure.1.
#hashtag index A0306: Bone marrow cells cultured with IL-3 for 7 days
Idents(mBC_new) <- "TagIDs" 
mBC_sub = subset(mBC_new, idents = c("A0306"))



#normalization
mBC_sub = NormalizeData(mBC_sub, scale.factor = 1000000)



#find variable features
mBC_sub = FindVariableFeatures(mBC_sub, mean.function = ExpMean, 
                                dispersion.function = LogVMR, 
                                mean.cutoff = c(0.1,Inf),
                                dispersion.cutoff = c(0.5,Inf))

#===============================Data Scaling===============================#
#scaling data 
#regression of read counts
all.genes <- rownames(mBC_sub)
mBC_sub = ScaleData(mBC_sub, features = all.genes)



#Run PCA
mBC_sub = RunPCA(mBC_sub, features = VariableFeatures(mBC_sub),
                  npcs = 100)



#regressing out the difference between the G2M and S phase scores
#perfomed as described by Seurat provider (https://satijalab.org/seurat/articles/cell_cycle_vignette.html)

#segregate this list into markers of G2/M phase and markers of S phase
#lists of cell cycle markers below were from Tirosh et al, 2015, with some minor modifications.
s.genes <- read.csv("sgeneslist.csv")
s.genes = as.vector(s.genes$gene)
g2m.genes <- read.csv("g2mgeneslist.csv")
g2m.genes = as.vector(g2m.genes$gene)

mBC_sub <- CellCycleScoring(mBC_sub, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

mBC_sub$CC.Difference <- mBC_sub$S.Score - mBC_sub$G2M.Score
mBC_sub <- ScaleData(mBC_sub, vars.to.regress = c("nReads", "CC.Difference"), features = rownames(mBC_sub))


#===============================PCA, Clustering, Dimensional Reduction===============================#
#principal component analysis
mBC_sub = RunPCA(mBC_sub, features = VariableFeatures(mBC_sub),
                 npcs = 100)



#JackStraw procedure
mBC_sub = JackStraw(mBC_sub, num.replicate = 100, dims = 100)
mBC_sub = ScoreJackStraw(mBC_sub, dims = 1:100)



#selecting principal components with their p-value less than 0.05
tmp = as.data.frame(mBC_sub@reductions$pca@jackstraw@overall.p.values)
tmp1 = tmp[tmp$Score>0.05, 1]
dims = c(1:(min(tmp1)-1))



#clustering and UMAP dimensional reduction
mBC_sub = FindNeighbors(mBC_sub, dims = dims)
mBC_sub = FindClusters(mBC_sub, resolution = 0.2)
mBC_sub = RunUMAP(mBC_sub, dims = dims)

#===============================Data Subsetting===============================#
#subset baso-lineage clusters
baso_lin = subset(mBC_sub, idents=c(0,2,10))

#change order of clusters
my_levels = c(2,0,10)
baso_lin@active.ident = factor(x=baso_lin@active.ident, levels = my_levels)

#change names of clusters 
current.cluster.ids <- c(2,0,10)
new.cluster.ids <- c('BMBA1', 'BMBA2', 'BMBA3')
baso_lin@active.ident <- plyr::mapvalues(x = baso_lin@active.ident, from = current.cluster.ids, to = new.cluster.ids)

#change colors of clusters
cols = c('BMBA1'="#0000FFFF",'BMBA2'="#FF0000",'BMBA3'="#FFCC33FF")

#===============================Detecting Differentially Expressed Genes===============================#
DEGs <- FindMarkers(mBC_sub, ident.1 = c("0","10"), ident.2 = c("2"), 
                      logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)

#===============================Visualization===============================#
#DimPlot
DimPlot(mBC_sub, reduction = "umap", pt.size = 1, label=T)+coord_fixed()
DimPlot(baso_lin, reduction = "umap", pt.size = 1, cols=cols, label=F)+xlim(-8,4)+ylim(-7, 2)+coord_fixed()



#FeaturePlot
FeaturePlot(mBC_sub, features = c("Mcpt8", "Cd200r3", "Kit", "Ly6c1", "Cd68", "Adgre1", "S100a8", "Tpsb2", "Elane", "Prg2", "Il5ra", "Flt3", "Cd79a"))

FeaturePlot(baso_lin, features = c("Fcer1a"), min.cutoff=5, pt.size=1)+xlim(-8,4)+ylim(-7, 2)+coord_fixed()
FeaturePlot(baso_lin, features = c("Cd9"), min.cutoff=5, pt.size=1)+xlim(-8,4)+ylim(-7, 2)+coord_fixed()
FeaturePlot(baso_lin, features = c("Itga2"), pt.size=1)+xlim(-8,4)+ylim(-7, 2)+coord_fixed()
FeaturePlot(baso_lin, features = c("Clec12a"), pt.size=1)+xlim(-8,4)+ylim(-7, 2)+coord_fixed()



#VlnPlot
VlnPlot(baso_lin, features = c("Fcer1a"), cols = cols, pt.size = 0)+NoLegend()
VlnPlot(baso_lin, features = c("Itga2"), cols = cols, pt.size = 0)+NoLegend()
VlnPlot(baso_lin, features = c("Clec12a"), cols = cols, pt.size = 0)+NoLegend()
VlnPlot(baso_lin, features = c("Cd9"), cols = cols, pt.size = 0)+NoLegend()



#Volcano Plot of DEGs
interested = c("Clec12a", "Cd9")
interested2 = c("Fcer1a", "Itga2", "Mcpt8", "Prss34")

DEGs$diffexpressed <- "NO"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
DEGs$diffexpressed[DEGs$avg_log2FC > 0.5 & DEGs$p_val_adj < 0.05] <- "UP"
# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
DEGs$diffexpressed[DEGs$avg_log2FC < -0.5 & DEGs$p_val_adj < 0.05] <- "DOWN"

DEGs$diffexpressed[DEGs$gene %in% interested] <- "INTEREST"
DEGs$diffexpressed[DEGs$gene %in% interested2] <- "INTEREST2"
DEGs$diffexpressed = factor(DEGs$diffexpressed, levels = c("INTEREST", "INTEREST2", "NO", "UP", "DOWN"))

p = ggplot(data=DEGs, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, size=diffexpressed)) +
  geom_point() + 
  scale_color_manual(values=c("black", "black", "gray", "deeppink1", "deepskyblue1"))+
  scale_size_manual(values=c(4,2.5, 1.5,1.5,1.5))+
  xlim(-6,6)+
  theme_classic()+
  theme(legend.position = "none")
p = LabelPoints(plot = p, points = c(interested, interested2), repel = TRUE)
p
