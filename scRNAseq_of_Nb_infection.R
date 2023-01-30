library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)


#Distibution-based background subtration,hashtag/sampletag demultiplexing, and Seurat Object creation for BD Rhapsody-derived WTA data was 
#conducted by using rDBEC package as previously described (Shichino, S. et al. Commun Biol. (2022). doi:10.1038/s42003-022-03536-0  Github: https://github.com/s-shichino1989/TASSeq-paper)

#===============================Quality Control and Selecting Cells===============================#
#perfomed by using R software package Seurat v2.3.4

#load Seurat Object pre-processed by rDBEC package (named as mBC)
load("MiyakeWTA5_Seurat.rda")



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
#We have performed multiplex assay of scRNA-seq analysis. We have extracted dataset for Figure.5.
#hashtag index A0312: skin basophils during second Nb infection
#hashtag index A0313: bone marrow basophils during second Nb infection
Idents(mBC_new) <- "TagIDs" 
mBC_sub = subset(mBC_new, idents = c("A0312", "A0313"))


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

DimPlot(mBC_sub, reduction = "umap", pt.size = 1, label=T)


#===============================Data Subsetting===============================#
#subset baso-lineage clusters
baso_lin = subset(mBC_sub, idents=c(0,1,2))

#change order of clusters
my_levels = c(1,0,2)
baso_lin@active.ident = factor(x=baso_lin@active.ident, levels = my_levels)

#change names of clusters 
current.cluster.ids <- c(1,0,2)
new.cluster.ids <- c('Pre-basophil', 'Mature1', 'Mature2')
baso_lin@active.ident <- plyr::mapvalues(x = baso_lin@active.ident, from = current.cluster.ids, to = new.cluster.ids)

#change colors of clusters
cols = c('Pre-basophil'='#0000FFFF','Mature1'='#FF0000FF','Mature2'='#FF6C00FF')


#subset Seurat object of the infected skin
Idents(baso_lin) <- "TagIDs" 
Skin = subset(baso_lin, idents = c("A0312"))

Idents(baso_lin) <- "seurat_clusters"
Idents(Skin) <- "seurat_clusters"

#change order of clusters
my_levels = c(1,0,2)
Skin@active.ident = factor(x=Skin@active.ident, levels = my_levels)

#change names of clusters 
current.cluster.ids <- c(1,0,2)
new.cluster.ids <- c('Pre-basophil', 'Mature1', 'Mature2')
Skin@active.ident <- plyr::mapvalues(x = Skin@active.ident, from = current.cluster.ids, to = new.cluster.ids)


#===============================calculate the proportion of each cluster among basophil-lineage cells===============================#
a = prop.table(table(Idents(baso_lin), baso_lin$TagIDs), margin = 2)
a


#===============================Detecting Differentially Expressed Genes===============================#
#detect differentially expressed genes
DEGs <- FindMarkers(Skin, ident.1 = c("Mature1", "Mature2"), ident.2 = c("Pre-basophil"), 
                    logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)
DEGs$gene <- rownames(DEGs)


#==============================GO Enrichment Analysis (Pre-baso, Mature)==============================#

#extract genesets upregulated/downregulated in mature basophils compared to pre-basophils
upDEGs = DEGs %>% filter(avg_log2FC> 1) %>% filter(p_val_adj< 0.05) %>% dplyr::select(avg_log2FC) %>% arrange(desc(avg_log2FC))
downDEGs = DEGs %>% filter(avg_log2FC< -1) %>% filter(p_val_adj< 0.05) %>% dplyr::select(avg_log2FC) %>% arrange(desc(avg_log2FC))

upgenes <- rownames(upDEGs)
downgenes <- rownames(downDEGs)



#convert gene symbols into ENTREZID
upgenes <- bitr(upgenes, 
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Mm.eg.db)

downgenes <- bitr(downgenes, 
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Mm.eg.db)



#GO BP enrichment analysis
ego_up <- enrichGO(gene = upgenes$ENTREZID,
                   OrgDb         = org.Mm.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05, 
                   readable      = TRUE)
head(as.data.frame(ego_up))
ego_up_simplified <- simplify(ego_up)
head(as.data.frame(ego_up_simplified))

ego_down <- enrichGO(gene = downgenes$ENTREZID,
                     OrgDb         = org.Mm.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05, 
                     readable      = TRUE)
head(as.data.frame(ego_down))
ego_down_simplified <- simplify(ego_down)
head(as.data.frame(ego_down_simplified))


#Visualization (dotplot)
clusterProfiler::dotplot(ego_up, showCategory = 10)
clusterProfiler::dotplot(ego_down_simplified, showCategory = 10)


#===============================Gene Set Enrichment Analysis===============================#
DEGlist = DEGs[,c("gene", "avg_log2FC")] %>% arrange(desc(avg_log2FC))
DEGgenes = DEGlist$gene
DEGgenes = bitr(DEGgenes, 
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Mm.eg.db)
DEGgenes = mutate(DEGgenes, gene=SYMBOL)
DEGgenes = DEGgenes[,c("gene", "ENTREZID")]

DEGlist = merge(DEGlist, DEGgenes, by="gene")
DEGname = DEGlist$ENTREZID
DEG = as.vector(DEGlist$avg_log2FC)
names(DEG) <- as.vector(DEGname)
DEG = sort(DEG, decreasing = T)
gseGO_diff <- gseGO(geneList     = DEG,
                    OrgDb        = org.Mm.eg.db,
                    ont          = "BP",@
                    pvalueCutoff = 0.05,
                    verbose      = FALSE)

#visualization of GSEA result ("chromosome segregation")
res = gseaplot2(gseGO_diff, geneSetID = "GO:0007059")
res


#===============================Visualization===============================#
#DimPlot
DimPlot(mBC_sub, reduction = "umap", pt.size = 1, split.by="TagIDs", label=T)+coord_fixed()
DimPlot(baso_lin, reduction = "umap", pt.size = 0.5, split.by="TagIDs",ncol=3, cols = cols)+ylim(-8,3)+coord_fixed()



#FeaturePlot
FeaturePlot(mBC_sub, features = c("Mcpt8", "Cd200r3", "Mrc1", "Arg1", "Cd3e", "Gzma", "Cd19", "S100a8", "Col1a1", "Ccr3"))

FeaturePlot(mBC_sub, features = c("Cd200r3"), pt.size=1)+coord_fixed()
FeaturePlot(mBC_sub, features = c("Mcpt8"), pt.size=1)+coord_fixed()

FeaturePlot(baso_lin, features = c("Fcer1a"), min.cutoff=6, pt.size=1)+ylim(-8,3)+coord_fixed()
FeaturePlot(baso_lin, features = c("Cd9"), min.cutoff=5, pt.size=1)+ylim(-8,3)+coord_fixed()
FeaturePlot(baso_lin, features = c("Itga2"), pt.size=1)+ylim(-8,3)+coord_fixed()
FeaturePlot(baso_lin, features = c("Clec12a"), pt.size=1)+ylim(-8,3)+coord_fixed()
FeaturePlot(baso_lin, features = c("Mcpt8"), pt.size=1)+ylim(-8,3)+coord_fixed()
FeaturePlot(baso_lin, features = c("Prss34"), pt.size=1)+ylim(-8,3)+coord_fixed()


#VlnPlot
VlnPlot(baso_lin, features = c("Fcer1a"), cols = cols, pt.size = 0)+NoLegend()
VlnPlot(baso_lin, features = c("Itga2"), cols = cols, pt.size = 0)+NoLegend()
VlnPlot(baso_lin, features = c("Clec12a"), cols = cols, pt.size = 0)+NoLegend()
VlnPlot(baso_lin, features = c("Cd9"), cols = cols, pt.size = 0)+NoLegend()
VlnPlot(baso_lin, features = c("Mcpt8"), cols = cols, pt.size = 0)+NoLegend()
VlnPlot(baso_lin, features = c("Prss34"), cols = cols, pt.size = 0)+NoLegend()
