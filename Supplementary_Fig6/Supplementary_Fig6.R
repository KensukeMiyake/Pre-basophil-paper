library(ggplot2)
library(Seurat)
library(dplyr)
library(stringr)
library(viridis)
library(future)

#===============================loading data from GSE140802===============================#
#create seurat Object
#Dataset is downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4185643. 
expression_matrix_stateFate_inVivo <- ReadMtx(mtx = "GSM4185643_stateFate_inVivo_normed_counts.mtx.gz",
                                    features = "GSM4185643_stateFate_inVivo_cell_barcodes.txt.gz",
                                    cells = "GSM4185643_stateFate_inVivo_gene_names.txt.gz",
                                    feature.column = 1)

expression_matrix_stateFate_inVivo =t(expression_matrix_stateFate_inVivo)
seurat_object_stateFate_inVivo <- CreateSeuratObject(counts = expression_matrix_stateFate_inVivo, project = "stateFate_inVivo")

metadata = seurat_object_stateFate_inVivo@meta.data
metadata_inVivo = read.table("GSM4185643_stateFate_inVivo_metadata.txt.gz", header = T, sep = "\t")
metadata_new = cbind(metadata, metadata_inVivo)
seurat_object_stateFate_inVivo@meta.data = metadata_new



#===============================Normalization, Variable Feature Detection, Data Scaling===============================#
mBC = NormalizeData(object = seurat_object_stateFate_inVivo, scale.factor=1000000)
all.genes <- rownames(mBC)
mBC = ScaleData(object = mBC, vars.to.regress = c("nCount_RNA"), features = all.genes, verbose = TRUE)
mBC = FindVariableFeatures(mBC, selection.method = "vst", nfeatures = 2000)



#===============================PCA, Dimensional Reduction, Clustering===============================#
#principal component analysis
mBC <- RunPCA(mBC, npcs = 100, verbose = FALSE)

#JackStraw procedure
mBC = JackStraw(mBC, num.replicate = 100, dims = 100) 
mBC = ScoreJackStraw(mBC, dims = 1:100)

#Dimensional Reduction, Clustering
tmp = as.data.frame(mBC@reductions$pca@jackstraw@overall.p.values)
tmp1 = tmp[tmp$Score>0.05, 1] #selecting principal components with their p-value less than 0.05
dims = c(1:(min(tmp1)-1))

mBC <- RunUMAP(mBC, reduction = "pca", dims = dims)
mBC <- FindNeighbors(mBC, reduction = "pca", dims = dims)
mBC <- FindClusters(mBC, resolution = 1.0)

# Visualization
DimPlot(mBC, reduction = "umap", label=T)+coord_fixed()
DimPlot(mBC, reduction = "umap", split.by = "Cell.type.annotation", label = T,ncol=5)+coord_fixed()

#FindAllMarkers
markers <- FindAllMarkers(mBC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


#===============================Reclustering of clusters including basophil-lineage cells===============================#
Baso_lin = subset (mBC, idents=c(12,1,35,16,25,20))
Baso_lin = NormalizeData(Baso_lin, scale.factor = 1000000)
Baso_lin = FindVariableFeatures(Baso_lin, mean.function = ExpMean, 
                               dispersion.function = LogVMR, 
                               mean.cutoff = c(0.1,Inf),
                               dispersion.cutoff = c(0.5,Inf))
all.genes <- rownames(Baso_lin)
Baso_lin = ScaleData(Baso_lin, vars.to.regress = c("nCount_RNA"), features = all.genes)
Baso_lin = RunPCA(Baso_lin, features = VariableFeatures(Baso_lin),npcs = 100)
Baso_lin = JackStraw(Baso_lin, num.replicate = 100, dims = 100)
Baso_lin = ScoreJackStraw(Baso_lin, dims = 1:100)
JackStrawPlot(Baso_lin, dims = 1:40)
ElbowPlot(object = Baso_lin)
tmp = as.data.frame(Baso_lin@reductions$pca@jackstraw@overall.p.values)
tmp1 = tmp[tmp$Score>0.05, 1]
dims = c(1:(min(tmp1)-1))
Baso_lin = FindNeighbors(Baso_lin, dims = dims)
Baso_lin = FindClusters(Baso_lin, resolution = 0.5)
Baso_lin = RunUMAP(Baso_lin, dims = dims)
DimPlot(Baso_lin, reduction = "umap", pt.size = 1, label=T)+coord_fixed()
DimPlot(Baso_lin, reduction = "umap", split.by="Cell.type.annotation", pt.size = 1,ncol=5)+coord_fixed()

DefaultAssay(Baso_lin) <- "RNA"
FeaturePlot(Baso_lin, features = c("Mcpt8","Clec12a","Fcer1a","Cd34","Kit","Cma1"),  ncol = 3, coord.fixed = T)
VlnPlot(Baso_lin, features = c("Mcpt8","Clec12a","Fcer1a","Cd34","Kit","Cma1"),  ncol = 3)
save(Baso_lin,dims,file="221124_Weinreb_baso_recluster.rda")

Baso = subset (Baso_lin, idents=c(13,4))
Baso = NormalizeData(Baso, scale.factor = 1000000)
Baso = FindVariableFeatures(Baso, mean.function = ExpMean, 
                               dispersion.function = LogVMR, 
                               mean.cutoff = c(0.1,Inf),
                               dispersion.cutoff = c(0.5,Inf))
all.genes <- rownames(Baso)
Baso = ScaleData(Baso, vars.to.regress = c("nCount_RNA"), features = all.genes)
Baso = RunPCA(Baso, features = VariableFeatures(Baso),npcs = 100)
Baso = JackStraw(Baso, num.replicate = 100, dims = 100)
Baso = ScoreJackStraw(Baso, dims = 1:100)
JackStrawPlot(Baso, dims = 1:40)
ElbowPlot(object = Baso)
tmp = as.data.frame(Baso@reductions$pca@jackstraw@overall.p.values)
tmp1 = tmp[tmp$Score>0.05, 1]
dims = c(1:(min(tmp1)-1))
Baso = FindNeighbors(Baso, dims = dims)
Baso = FindClusters(Baso, resolution = 0.6)
Baso = RunUMAP(Baso, dims = dims)

require(scales)
# Create vector with levels of object@ident
identities <- levels(Baso@active.ident)
# Create vector of default ggplot2 colors
my_color_palette <- hue_pal()(length(identities))
my_color_palette
cols = c('0'="#A58AFF", "1"="#00C094", '2'="#C49A00", '3'="#53B400",'4'="#00B6EB", '5'="#F8766D")
Bas = subset (Baso, idents=c(0,1,2,3,4,5))
my_levels = c(5,2,3,1,4,0)
Baso@active.ident = factor(x=Baso@active.ident, levels = my_levels)

DimPlot(Bas, reduction = "umap", pt.size = 1, label=T, cols=cols)+coord_fixed()
DimPlot(Baso, reduction = "umap", split.by="Cell.type.annotation", pt.size = 1,ncol=3)+coord_fixed()
DimPlot(Baso, reduction = "umap", split.by="TagIDs", pt.size = 1,ncol=1)
VlnPlot(Baso, features = "Cd34")
DefaultAssay(Baso) <- "RNA"
FeaturePlot(Bas, features = c("Kit","Cd34","Clec12a","Cd9"),  ncol = 2, coord.fixed = T)

VlnPlot(Bas, features = c("Kit","Cd34","Clec12a","Cd9"), cols=cols, ncol = 2, pt.size = 0)
a = DotPlot(Baso,  features = c("Kit","Il1rl1","Cd34","Clec12a","Mcpt8","Gata1"))
table(Idents(Baso), Baso$orig.ident)
View(a$data)
save(Baso,dims,file="221124_integrated_baso_recluster.rda")
load("Weinreb et al 2020 Science/221124_integrated_baso_recluster.rda")


#Monocle3 pseudotime analysis
current.cluster.ids <- c(5,2,3,1,4,0,6)
new.cluster.ids <- c('pre-BMP-like', 'Immature1', 'Immature2', 'Mature1', 'Mature2','Mature3', 'non-Baso')
Baso@active.ident <- plyr::mapvalues(x = Baso@active.ident, from = current.cluster.ids, to = new.cluster.ids)
