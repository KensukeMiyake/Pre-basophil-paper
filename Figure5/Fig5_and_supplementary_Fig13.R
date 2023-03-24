library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

# installation of rDBEC package from Github: https://github.com/s-shichino1989/TASSeq-paper)
# is required for the reproduction of results.
library(rDBEC)

#Distibution-based background subtration,hashtag/sampletag demultiplexing, and Seurat Object creation for 
# BD Rhapsody-derived WTA data was conducted by using rDBEC package as previously described 
# Shichino, S. et al. Commun Biol. (2022). doi:10.1038/s42003-022-03536-0
# Github: https://github.com/s-shichino1989/TASSeq-paper)

#=================Background subtraction by DEBC filtering and create Seurat object===============================#
setwd("/media/owner/Data/scRNAseq_Nb_second/GSE206590_Nb_infection_scRNA-seq/processed file")
# read processed count data matrix ("matrix_inflection_MiyakeWTA5.txt.gz")
tablelist = lapply(as.vector("matrix_inflection_MiyakeWTA5.txt.gz"), tableread_fast_sparse)
fnames = "MiyakeWTA5"
names(tablelist)=fnames

# DBEC filtering by using rDBEC package
DBEC_filter = background_subtraction_Biex(tablelist, min.event=100, minimum.max.expr=7,
                                          species="mmu", min.ave=6.0, 
                                          min.diff=5.5, modelnames="E",
                                          uncert.thre=1, nthreads=12, sample.name=fnames)
DBEC_res = apply_DBEC_filter(tablelist, DBEC_filter=DBEC_filter, nthreads=24, sample.name = fnames)
names(DBEC_res) = fnames

# create Seurat object (as list)
for (i in 1:length(DBEC_res)){
  colnames(DBEC_res[[i]]) = paste(names(DBEC_res)[i], colnames(DBEC_res[[i]]), sep='_')
}
seu = list(CreateSeuratObject(count=DBEC_res[[1]], min.cells = 5))

#========================Add hashtag annotation metadata to Seurat object===============================#
tableread_fast_hashtag = function(x, sep="\t", header=TRUE){
  tmp = data.table::fread(x, header=header, sep=sep, quote="")
  tmp = as.data.frame(tmp)
  return(tmp)
}
hashtag_data=lapply(as.vector("Hashtag_top1M_MiyakeTag5.txt.gz"), tableread_fast_hashtag)
for (i in 1:length(hashtag_data)){
  hashtag_data[[i]][,1] = paste(names(tablelist)[i], hashtag_data[[i]][,1], sep='_')
}

for (i in 1:length(hashtag_data)){
  seu[[i]] = Demultiplex_DNAtags(hashtag_data[[i]], seu[[i]],
                                 nTags=4, nCells_in_cartridge=20000, 
                                 sample.name=fnames[[i]], scale.factor = 4001613,
                                 dir.name="./")
}

mBC = seu[[1]][[1]]
dim(mBC@meta.data)
mBC$percent.mito <- PercentageFeatureSet(object = mBC, pattern = "mt-")
mBC$percent.mito <- mBC@meta.data$percent.mito / 100

metadata = mBC@meta.data
metadata <- metadata %>%
  dplyr::rename(nReads = nCount_RNA, nGene = nFeature_RNA)
mBC@meta.data = metadata 
View(metadata)

#add hashtag ID to cell_barcodes and export expression matrix 
raw.data = mBC@assays$RNA@counts
tmp = paste(mBC@meta.data$TagIDs, rownames(mBC@meta.data), sep="_")
colnames(raw.data)=tmp
colnames(raw.data)=gsub("not_detected", "not-detected", colnames(raw.data))
raw.data=as.matrix(raw.data)
fwrite(as.data.frame(raw.data), file = "matrix_inflection_demulti_DBEC_miyakeWTA5.txt.gz", 
       row.names=T, col.names=T, sep="\t", eol="\n", quote=F, compress="gzip")

#===============================Quality Control and Selecting Cells===============================#
#perfomed by using R software package Seurat v4.0.4
expression_matrix = read.table("matrix_inflection_demulti_DBEC_miyakeWTA5.txt.gz")
cellnames = str_split(colnames(expression_matrix), pattern = "_", n=2, simplify = T)
TagIDs = cellnames[,1]
colnames(expression_matrix) = cellnames[,2]
View(expression_matrix)

mBC <- CreateSeuratObject(counts = expression_matrix, project = "miyake_WTA5")

mBC$percent.mito <- PercentageFeatureSet(object = mBC, pattern = "mt-")
mBC$percent.mito <- mBC@meta.data$percent.mito / 100

metadata = mBC@meta.data
metadata <- metadata %>%
  dplyr::rename(nReads = nCount_RNA, nGene = nFeature_RNA)
metadata$TagIDs = TagIDs
mBC@meta.data = metadata 
View(metadata)

# QC Filtering
# Filter out low quality reads using selected thresholds
mBC <- subset(x = mBC, (nReads >= 500) & (nGene >= 500) )
mBC <- subset(x = mBC, (TagIDs == "A0312") | (TagIDs == "A0313") | (TagIDs == "A0314") |  (TagIDs == "A0315") )
table(mBC$TagIDs)
dim(mBC@meta.data)

#===============================Seurat Object Updating, Normalization, Variable Feature Detection===============================#
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
mBC_sub = RunPCA(mBC_sub, features = VariableFeatures(mBC_sub),npcs = 100)

#regressing out the difference between the G2M and S phase scores
#perfomed as described by Seurat provider (https://satijalab.org/seurat/articles/cell_cycle_vignette.html)

]

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
                    ont          = "BP",?@
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
