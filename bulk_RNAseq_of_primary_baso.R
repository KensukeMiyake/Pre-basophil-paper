library(TCC)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(org.Mm.eg.db)

#==============================Load Count Table (BaP, Pre-baso, Mature, Spl baso)==============================#
#load count table
data = read.csv("mouse_5Seq_table.csv", header=T)

rownames(data) = data[,c(1)]
data = data[, c(2:13)]

#==============================TCC Object Construction and Gene Filtering (BaP, Pre-baso, Mature, Spl baso)==============================#
#construct TCC cluss object
group = c(1,1,1,2,2,2,3,3,3,4,4,4)
tcc = new("TCC", data, group)



#filter low-count genes
tcc <- filterLowCountGenes(tcc)

#==============================Normalization and DEG Detection (BaP, Pre-baso, Mature, Spl baso)==============================#
#normalization of multi-group count data
tcc = calcNormFactors(tcc, norm.method = "tmm", test.method = "edger", iteration = 3)

normalized.count <- getNormalizedData(tcc)



#DE analysis
tcc <- estimateDE(tcc, test.method = "edger", FDR = 0.1)
result <- getResult(tcc, sort = TRUE)

#==============================Heatmap (BaP, Pre-baso, Mature, Spl baso)==============================#
#create the normalized count data of DEGs
normalized.count = as.data.frame(normalized.count)
normalized.count = mutate(normalized.count, gene_id=rownames(normalized.count))
norm.count = merge(normalized.count, result, by="gene_id")
norm.count = filter(norm.count, estimatedDEG==1)
rownames(norm.count) = norm.count$gene_id
norm.count = norm.count[, c(2:13)]



#Z-score calculation
zscore <- t(apply(norm.count, 1, scale)) 
colnames(zscore) = colnames(norm.count)
head(zscore)
zscore = as.matrix(zscore)



#Heatmap annotation
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

Group = c(rep(1,"BaP", 3), rep(2,"Pre-baso", 3), rep(3,"Mature", 3), rep(4,"Spleen", 3))
cols = c('1'='#777777', '2'='blue', '3'='red', '4'='purple')

ha = HeatmapAnnotation(Group = Group,
                       col = list(Group = cols))

Gene_list=c("Clec12a", "Itga2", "Fcer1a", "Cd34", "Cd9", "Mcpt8", "Prss34")
Gene_position=which(rownames(zscore) %in% Gene_list)


ha1 = rowAnnotation(link = anno_mark(at = Gene_position, labels = rownames(zscore)[Gene_position],
                                     labels_gp = gpar(fontsize = 10), padding = unit(0.01, "mm")))



#Depicting Heatmap
ht_list = Heatmap(zscore, name = "zscore", col=col_fun,
        show_row_dend  = T, cluster_columns = F, 
        clustering_method_rows = "ward.D2",
        show_row_names = F,
        bottom_annotation = ha,
        right_annotation = ha1)
draw(ht_list, heatmap_legend_side="left", annotation_legend_side = "right")









#==============================Load Count Table (Pre-baso, Mature)==============================#
#load count table
data = data[, c(4:9)]

#==============================TCC Object Construction and Gene Filtering (Pre-baso, Mature)==============================#
#construct TCC cluss object
group = c(2,2,2,3,3,3)
tcc = new("TCC", data, group)

tcc <- filterLowCountGenes(tcc)

#==============================Normalization and DEG Detection (Pre-baso, Mature)==============================#
#normalization of two-group count data
tcc = calcNormFactors(tcc, norm.method = "deseq2", test.method = "deseq2", iteration = 3)

normalized.count <- getNormalizedData(tcc)



#DE analysis
tcc <- estimateDE(tcc, test.method = "deseq2", FDR = 0.1)
result <- getResult(tcc, sort = TRUE)

#==============================GO Enrichment Analysis (Pre-baso, Mature)==============================#

#extract genesets upregulated/downregulated in mature basophils compared to pre-basophils
result$diffexpressed <- "NO"
result$diffexpressed[result$m.value > 0 & result$estimatedDEG==1] <- "UP"
result$diffexpressed[result$m.value < 0 & result$estimatedDEG==1] <- "DOWN"

downgenes <- result$gene_id[result$diffexpressed=="DOWN"]
upgenes <- result$gene_id[result$diffexpressed=="UP"]



#convert gene symbols into ENTREZID
downgenes <- bitr(downgenes, 
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)

upgenes <- bitr(upgenes, 
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Mm.eg.db)



#GO BP enrichment analysis
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



#Visualization (dotplot)
clusterProfiler::dotplot(ego_down_simplified, showCategory = 10)
clusterProfiler::dotplot(ego_up_simplified, showCategory = 10)

