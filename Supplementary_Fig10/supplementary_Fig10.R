library(TCC)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)

#Set working directory
setwd("/media/owner/Data/221026 Bulk RNA-seq immature and mature BMBA stimulation")

#Read expression table data
data1 = read.table("/media/owner/Data/221026 Bulk RNA-seq immature and mature BMBA stimulation/Miyake1121/Miyake1121_RNASeq_expression_table.txt", header=T)
data2 = read.table("/media/owner/Data/221026 Bulk RNA-seq immature and mature BMBA stimulation/Miyake1122/Miyake1122_RNASeq_expression_table.txt", header=T)
colnames(data1)
colnames(data1) = c("Symbol","Control_prebaso_1","Control_prebaso_2","Control_prebaso_3","Control_mature_1","Control_mature_2","Control_mature_3",
                    "TNP_prebaso_1","TNP_prebaso_2","TNP_prebaso_3")
head(data1)

colnames(data2)
colnames(data2) = c("Symbol","TNP_mature_1","TNP_mature_2","TNP_mature_3","IL3_prebaso_1","IL3_prebaso_2","IL3_prebaso_3",
                    "IL3_mature_1","IL3_mature_2","IL3_mature_3")
head(data2)
data_integrated = merge(data1, data2, by="Symbol", all = F)

#Treat NA as 0
data_integrated <- mutate_all(data_integrated, ~replace(., is.na(.), 0))
head (data_integrated)
dim(data_integrated)
#save data_integrated as tab-delineated text file for GEO submission
write.table (data_integrated, file="Miyake_1121_1122_RNAseq_expression_table.txt", quote= F, sep = "\t")


data_integrated = data_integrated[, c(2:19)]
head (data_integrated)

write.table (data_integrated, file="Miyake_1121_1122_RNAseq_expression_table.txt", quote= F, sep = "\t")

#construct TCC cluss object
group = c(rep(1,"OVA_immature", 3), rep(2,"OVA_mature", 3), rep(3,"TNP_immature", 3), 
          rep(4,"TNP_mature", 3),rep(5,"IL3_immature", 3),rep(6,"IL3_mature", 3))
tcc = new("TCC", data_integrated, group)
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
write.csv(normalized.count, file = "221027_Bulk_mature_immature_stim_normalized_count.csv")
head (normalized.count)
dim(normalized.count)
normalized.count[rownames(normalized.count)=="Il6",]

#DE analysis
designall <- model.matrix(~as.factor(group))
designall
tcc <- estimateDE(tcc, test.method="edger", FDR=0.1,
                  design=designall,coef=c(2:6))
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
mean_G5 <- log2(apply(sweep(as.matrix(normalized.count[,group==5]), 1, 0.125, FUN="+"),1,mean))
mean_G6 <- log2(apply(sweep(as.matrix(normalized.count[,group==6]), 1, 0.125, FUN="+"),1,mean))

FoldChange <-cbind (mean_G1,mean_G2,mean_G3,mean_G4,mean_G5,mean_G6)
head(FoldChange)
MAX_FC<- apply(FoldChange,1,max)
MIN_FC<- apply(FoldChange,1,min)
FoldChange = cbind( FoldChange,MAX_FC,MIN_FC)
tmp<- cbind(tmp,FoldChange)
colnames(tmp)
tmp =tmp[,-c(19:21)]
head(tmp)
write.csv(tmp, file="221027_DEG_TCC_all.csv", quote=F, row.names=T)

#Heatmap
#create the normalized count data of DEGs
norm.count = filter(tmp, tmp$estimatedDEG==1)
head(norm.count)
norm.count = norm.count[, c(1:18)]
head(norm.count)

# get log normalized counts and change to z-score
zscore <- t(apply(norm.count, 1, scale))
colnames(zscore) = colnames(norm.count)
head(zscore)
zscore <- as.matrix(zscore)

#heatmap colors
library(circlize)
library(dendextend)

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
group = c(rep(1,"OVA_immature", 3), rep(2,"OVA_mature", 3), rep(3,"TNP_immature", 3), 
          rep(4,"TNP_mature", 3),rep(5,"IL3_immature", 3),rep(6,"IL3_mature", 3))
cols = c('1'="blue", '2'="red", '3'="turquoise",'4'="pink",'5'="purple",'6'="green")

Gene_list=c("Il4","Il6","Il10","Il13","Csf1","Ccl3","Tnf")
Gene_position=which(rownames(zscore) %in% Gene_list)
Gene_position
ha1 = rowAnnotation(link = anno_mark(at = Gene_position, labels = rownames(zscore)[Gene_position],
                                     labels_gp = gpar(fontsize = 10), padding = unit(0.01, "mm")))

ha <- HeatmapAnnotation(Group = group,
                        col = list(Group = cols))
ht_list = Heatmap(zscore, name = "Z-score", col=col_fun,
        show_row_dend  = T, cluster_columns =  F, 
        clustering_method_rows = "ward.D2",
        show_row_names = F,
        bottom_annotation = ha,
        right_annotation = ha1,
        column_order = c("OVA_immature_1","OVA_immature_2","OVA_immature_3","TNP_immature_1","TNP_immature_2","TNP_immature_3","IL3_immature_1","IL3_immature_2","IL3_immature_3",
                         "OVA_mature_1","OVA_mature_2","OVA_mature_3","TNP_mature_1","TNP_mature_2","TNP_mature_3","IL3_mature_1","IL3_mature_2","IL3_mature_3"),
        use_raster = F)
draw(ht_list, heatmap_legend_side="left", annotation_legend_side = "left")

