library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(monocle3)

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
#We have performed multiplex assay of scRNA-seq analysis for Figure.1&2. We have extracted dataset for Figure.2.
#hashtag index A0301: Bone marrow basophils from Mcpt8-GFP Tg mice
#hashtag index A0302: Spleen basophils from Mcpt8-GFP Tg mice
Idents(mBC_new) <- "TagIDs" 
mBC_sub = subset(mBC_new, idents = c("A0301", "A0302"))



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
mBC_sub = FindClusters(mBC_sub, resolution = 1.0)
mBC_sub = RunUMAP(mBC_sub, dims = dims)

#===============================Data Subsetting and Reclustering===============================#
#subset baso-lineage clusters
baso_lin = subset(mBC_sub, idents=c(0,1,4))

#recluster baso-lineage clusters
baso_lin = NormalizeData(baso_lin, scale.factor = 1000000)
baso_lin = FindVariableFeatures(baso_lin, mean.function = ExpMean, 
                               dispersion.function = LogVMR, 
                               mean.cutoff = c(0.1,Inf),
                               dispersion.cutoff = c(0.5,Inf))
baso_lin <- ScaleData(baso_lin, features = rownames(baso_lin)) 
baso_lin = RunPCA(baso_lin, features = VariableFeatures(baso_lin),
                 npcs = 100)
baso_lin = JackStraw(baso_lin, num.replicate = 100, dims = 100)
baso_lin = ScoreJackStraw(baso_lin, dims = 1:100)

tmp = as.data.frame(baso_lin@reductions$pca@jackstraw@overall.p.values)
tmp1 = tmp[tmp$Score>0.05, 1]
dims = c(1:(min(tmp1)-1))

baso_lin = FindNeighbors(baso_lin, dims = dims)
baso_lin = FindClusters(baso_lin, resolution = 0.2)
baso_lin = RunUMAP(baso_lin, dims = dims)

#DimPlot
DimPlot(baso_lin, reduction = "umap", pt.size = 1, label=T)

#change order of clusters
my_levels = c(3,0,1,2)
baso_lin@active.ident = factor(x=baso_lin@active.ident, levels = my_levels)

#change names of clusters 
current.cluster.ids <- c(3,0,1,2)
new.cluster.ids <- c('pre-BMP-like', 'Baso1', 'Baso2', 'Baso3')
baso_lin@active.ident <- plyr::mapvalues(x = baso_lin@active.ident, from = current.cluster.ids, to = new.cluster.ids)

#change colors of clusters
cols = c('pre-BMP-like'='#BAD3FFFF','Baso1'='#0000FFFF','Baso2'='#FF0000','Baso3'='#FFCC33FF')




#subset baso1 cluster
prebaso = subset(baso_lin, idents=c("Baso1"))

#recluster baso-lineage clusters
prebaso = NormalizeData(prebaso, scale.factor = 1000000)
prebaso = FindVariableFeatures(prebaso, mean.function = ExpMean, 
                                dispersion.function = LogVMR, 
                                mean.cutoff = c(0.1,Inf),
                                dispersion.cutoff = c(0.5,Inf))
prebaso <- ScaleData(prebaso, features = rownames(prebaso)) 
prebaso = RunPCA(prebaso, features = VariableFeatures(prebaso),
                  npcs = 100)
prebaso = JackStraw(prebaso, num.replicate = 100, dims = 100)
prebaso = ScoreJackStraw(prebaso, dims = 1:100)

tmp = as.data.frame(prebaso@reductions$pca@jackstraw@overall.p.values)
tmp1 = tmp[tmp$Score>0.05, 1]
dims = c(1:(min(tmp1)-1))

prebaso = FindNeighbors(prebaso, dims = dims)
prebaso = FindClusters(prebaso, resolution = 0.2)
prebaso = RunUMAP(prebaso, dims = dims)

#DimPlot
DimPlot(prebaso, reduction = "umap", pt.size = 1, label=F)+coord_fixed()

#===============================Monocle3 Pseudotime===============================#
#Script is from user rwo012 on Github: https://github.com/satijalab/seurat/issues/1658 
#and user tlusardi on Github: https://github.com/cole-trapnell-lab/monocle-release/issues/388, 
#as previously described (de Jong MME, et al. Nat Immunol. (2021) doi:10.1038/s41590-021-00931-3, Github: https://github.com/MyelomaRotterdam/Microenvironment) 

# Create an expression matrix
expression_matrix <- baso_lin@assays$RNA@counts

# Get cell metadata
cell_metadata <- baso_lin@meta.data
if (all.equal(colnames(expression_matrix), rownames(cell_metadata))) {
  print(sprintf("Cell identifiers match"))
} else {
  print(sprintf("Cell identifier mismatch - %i cells in expression matrix, %i cells in metadata",
                ncol(expression_matrix), nrow(cell_metadata)))
  print("If the counts are equal, sort differences will throw this error")
}

# get gene annotations
gene_annotation <- data.frame(gene_short_name = rownames(baso_lin@assays$RNA), row.names = rownames(baso_lin@assays$RNA))
if (all.equal(rownames(expression_matrix), rownames(gene_annotation))) {
  print(sprintf("Gene identifiers all match"))
} else {
  print(sprintf("Gene identifier mismatch - %i genes in expression matrix, %i gene in gene annotation",
                nrow(expression_matrix), nrow(gene_annotation)))
  print("If the counts are equal, sort differences will throw this error")
}

# Seurat-derived CDS
my.cds <- new_cell_data_set(expression_matrix,
                            cell_metadata = cell_metadata,
                            gene_metadata = gene_annotation)


# Transfer Seurat embeddings
reducedDim(my.cds, type = "PCA") <- baso_lin@reductions$pca@cell.embeddings 
my.cds@preprocess_aux$prop_var_expl <- baso_lin@reductions$pca@stdev
plot_pc_variance_explained(my.cds)

# Transfer Seurat UMAP embeddings
my.cds@int_colData@listData$reducedDims$UMAP <- baso_lin@reductions$umap@cell.embeddings
plot_cells(my.cds)


# Copy cluster info from Seurat
my.cds@clusters$UMAP_so$clusters <- baso_lin@meta.data$seurat_clusters

my.cds <- cluster_cells(my.cds, reduction_method = "UMAP", resolution = 1e-3)

# Fix from https://gitmemory.com/cole-trapnell-lab
rownames(my.cds@principal_graph_aux$UMAP$dp_mst) <- NULL
colnames(my.cds@int_colData@listData$reducedDims$UMAP) <- NULL

DimPlot(baso_lin, reduction = "umap")
plot_cells(my.cds, color_cells_by = "partition", group_label_size = 3.5)
plot_cells(my.cds, color_cells_by = "seurat_clusters", show_trajectory_graph = F, group_label_size = 3.5)
my.cds = learn_graph(my.cds)
plot_cells(my.cds,
           color_cells_by = "seurat_clusters",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)


my.cds = order_cells(my.cds, reduction_method = "UMAP")
plot_cells(my.cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           cell_size = 1,
           trajectory_graph_color = "#000000")+coord_fixed()






genes <- c("Clec12a", "Mcpt8", "Fcer1a", "Prss34")

lineage_cds <- my.cds[rowData(my.cds)$gene_short_name %in% genes]

cluster_cols = c('#0000FFFF', '#FF0000', '#FFCC33FF','#BAD3FFFF')
plot_genes_in_pseudotime(lineage_cds,
                         color_cells_by="seurat_clusters",
                         cell_size=0.5, 
                         panel_order = c("Clec12a", "Mcpt8", "Fcer1a", "Prss34"),
                         min_expr=0.5, ncol=2)+
scale_color_manual(values = cluster_cols, name = "cluster")

#===============================Visualization===============================#
#DimPlot
DimPlot(mBC_sub, reduction = "umap", pt.size = 0.5, split.by="TagIDs", label=T)+coord_fixed()
DimPlot(baso_lin, reduction = "umap", pt.size = 0.5, split.by="TagIDs", ncol=1, cols = cols)+NoLegend()
DimPlot(prebaso, reduction = "umap", pt.size = 1, label=F)+coord_fixed()


#FeaturePlot
FeaturePlot(mBC_sub, features = c("Mcpt8", "Cd200r3", "Elane", "S100a8", "Lyz2", "Il5ra", "Cd34", "Kit", "Cd68", "Cd19", "Cma1"))

FeaturePlot(baso_lin, features = c("Itga2"), pt.size=1)+coord_fixed()
FeaturePlot(baso_lin, features = c("Clec12a"), pt.size=1)+coord_fixed()
FeaturePlot(baso_lin, features = c("Mcpt8"), pt.size=1)+coord_fixed()
FeaturePlot(baso_lin, features = c("Prss34"), pt.size=1)+coord_fixed()
FeaturePlot(baso_lin, features = c("Cd34"), pt.size=1)+coord_fixed()
FeaturePlot(baso_lin, features = c("Kit"), pt.size=1)+coord_fixed()
FeaturePlot(baso_lin, features = c("Cd9"), min.cutoff=5, pt.size=1)+coord_fixed()
FeaturePlot(baso_lin, features = c("Fcer1a"), min.cutoff=6, pt.size=1)+coord_fixed()

#FeaturePlot for S phase score
baso_lin = AddModuleScore(baso_lin, features = s_phase_genes, name = "S_phase_score")
FeaturePlot(baso_lin, features = c("S_phase_score1"), min.cutoff = 0, pt.size = 1)+coord_fixed()

#VlnPlot
VlnPlot(baso_lin, features = c("Clec12a"), cols = cols, pt.size = 0)+NoLegend()
VlnPlot(baso_lin, features = c("Cd9"), cols = cols, pt.size = 0)+NoLegend()
VlnPlot(baso_lin, features = c("Fcer1a"), cols = cols, pt.size = 0)+NoLegend()
VlnPlot(baso_lin, features = c("Cd34"), cols = cols, pt.size = 0)+NoLegend()
VlnPlot(baso_lin, features = c("Kit"), cols = cols, pt.size = 0)+NoLegend()

VlnPlot(prebaso, features = c("Cd34"), pt.size=0)+NoLegend()