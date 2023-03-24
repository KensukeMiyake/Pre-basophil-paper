library(TCC)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(dendextend)

#Set working directory
setwd("/media/owner/Data/GEO submission fiiles/geo_submission_IL3_vs_TSLP_basophils/")

#Read expression table data
data = read.table("count_table_IL3_vs_TSLP_rev.txt", header=T)
colnames(data)
head(data)
colnames(data)
rownames(data) = data[,1]
data = data[, c(c(5:7),c(2:4),c(11:13),c(8:10))]
head(data)

#construct TCC cluss object
group = c(rep(1,"IL3_pre-baso", 3), rep(2,"IL3_mature", 3), 
          rep(3,"TSLP_pre-baso", 3),rep(4,"TSLP_mature", 3))
tcc = new("TCC", data, group)
head(tcc)

#filter low-count genes
dim(tcc$count)
tcc <- filterLowCountGenes(tcc)
dim(tcc$count)
#normalization of multi-group count data
tcc <- calcNormFactors(tcc, norm.method="tmm", test.method="edger",
                       iteration=3, FDR=0.1, floorPDEG=0.05)
tcc$norm.factors
normalized.count <- getNormalizedData(tcc)
write.csv(normalized.count, file = "Bulk_IL3_vs_TSLP_normalized_count.csv")
head (normalized.count)
dim(normalized.count)
normalized.count[rownames(normalized.count)=="Il6",]

#DE analysis
designall <- model.matrix(~as.factor(group))
designall
tcc <- estimateDE(tcc, test.method="edger", FDR=0.1,
                  design=designall,coef=c(2:4))
result <- getResult(tcc, sort=FALSE)
head(result)
sum(tcc$stat$q.value < 0.1) 
tmp <- cbind( normalized.count, result)
head(tmp)

#Fold Change
mean_G1 <- log2(apply(sweep(as.matrix(normalized.count[,group==1]), 1, 0.125, FUN="+"),1,mean))
mean_G2 <- log2(apply(sweep(as.matrix(normalized.count[,group==2]), 1, 0.125, FUN="+"),1,mean))
mean_G3 <- log2(apply(sweep(as.matrix(normalized.count[,group==3]), 1, 0.125, FUN="+"),1,mean))
mean_G4 <- log2(apply(sweep(as.matrix(normalized.count[,group==4]), 1, 0.125, FUN="+"),1,mean))

FoldChange <-cbind (mean_G1,mean_G2,mean_G3,mean_G4)
head(FoldChange)
MAX_FC<- apply(FoldChange,1,max)
MIN_FC<- apply(FoldChange,1,min)
FoldChange = cbind( FoldChange,MAX_FC,MIN_FC)
tmp<- cbind(tmp,FoldChange)
colnames(tmp)
tmp =tmp[,-c(14:15)]
head(tmp)
write.csv(tmp, file="221027_DEG_TCC_all.csv", quote=F, row.names=T)

#Heatmap
#create the normalized count data of DEGs
norm.count = filter(tmp, tmp$estimatedDEG==1)
head(norm.count)
norm.count = norm.count[, c(1:12)]
head(norm.count)

# get log normalized counts and change to z-score
zscore <- t(apply(norm.count, 1, scale))
colnames(zscore) = colnames(norm.count)
head(zscore)
zscore <- as.matrix(zscore)

#heatmap colors
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

cols = c('1'="blue", '2'="red", '3'="turquoise",'4'="pink")


ha <- HeatmapAnnotation(Group = group,
                        col = list(Group = cols))
ht_list = Heatmap(zscore, name = "Z-score", col=col_fun,
        show_row_dend  = T, cluster_columns =  F, 
        clustering_method_rows = "ward.D2",
        show_row_names = F,
        bottom_annotation = ha,
        use_raster = F)
draw(ht_list, heatmap_legend_side="left", annotation_legend_side = "left")

