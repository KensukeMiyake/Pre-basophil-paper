library(TCC)
library(dplyr)
library(ggplot2)
library(ggrepel)


#==============================Load Count Table==============================#
#load count table
data = read.csv("mouse_5Seq_table_BMBA.csv", header=T)
rownames(data) = data[,c(1)]
data = data[, c(2:7)]


#==============================TCC Object Construction and Gene Filtering==============================#
#construct TCC cluss object
group = c(1,1,1,2,2,2)
tcc = new("TCC", data, group)



#filter low-count genes
tcc <- filterLowCountGenes(tcc)



#==============================Normalization and DEG Detection==============================#
#normalization of two-group count data
tcc = calcNormFactors(tcc, norm.method = "tmm", test.method = "edger", iteration = 3)



#DE analysis
tcc <- estimateDE(tcc, test.method = "edger", FDR = 0.1)
result <- getResult(tcc, sort = TRUE)

#==============================Visualization (MA plot)==============================#
result$diffexpressed <- "NO"
result$diffexpressed[result$m.value > 1 & result$estimatedDEG==1] <- "UP"
result$diffexpressed[result$m.value < -1 & result$estimatedDEG==1] <- "DOWN"
result$diffexpressed[result$gene_id == "Cxcr4"] <- "INTEREST"
result$diffexpressed = factor(result$diffexpressed, levels = c("NO", "UP", "DOWN", "INTEREST"))


result$delabel <- NA
result$delabel[result$gene_id == "Cxcr4"] <- result$gene_id[result$gene_id == "Cxcr4"]


#MA plot
ggplot(data=result, aes(x=a.value, y=m.value, col=diffexpressed, size=diffexpressed, label=delabel))+
  geom_point() +
  geom_text_repel() +
  scale_color_manual(values=c("gray", "gray40", "gray40", "red")) +
  scale_size_manual(values=c(1,1,1,2))+
  theme_classic()+
  theme(legend.position = "none")
