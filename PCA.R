# 载入程序包
library(actuar)
library(TeachingDemos) 
library(VennDiagram)
library(openxlsx)
library(dplyr)
library(plyr)
library(ggplot2)
library(gtable)
library(grid)
library(cowplot)
library(RColorBrewer)
suppressMessages(library(foreach))
suppressMessages(library(doSNOW))
suppressMessages(library(doParallel))
library(xml2)
library(httr)
library(jsonlite)
library(rvest)
library(RCurl)
library(FactoMineR)
library(car)
library(WGCNA)
library(getopt)
library(plotrix)
library(lemon)
library(impute)
library(myplot)
library(Mfuzz)
library(stringr)
library(plotly)
library(htmlwidgets)
library(ggpubr)
library(ggrepel)
library(igraph)


#PCA绘图
# ******************************************************************************************
#颜色
#####FF99CC00表示透明
colpalette <- c("#1d953f", "#102b6a", "#c77eb5", "#2585a6", "purple", "#e0861a", "#d71345", "#6b473c", "#ffce7b", "#78a355", "#fdb933", "#5e7c85", "#411445", "#c37e00", "#bed742", "#009ad6", "#9d9087", "#aa2116", "#225a1f", "#f3715c", "#7bbfea", "#dbce8f", "#f8aba6", "#778899", "#2F4F4F", "#2E8B57", "#8B6969", "#EE2C2C", "#6495ED", "#87CEEB", "#0000CD", "#63B8FF", "#006400", "#00FF7F", "#9400D3", "#FFAEB9", "#8B5A2B", "#FFFF00", "#BDB76B", "#FF1493")
colpalette1 <- c("#1d953f","#102b6a")
pca.plot <- function(pca.path, ploMN0, samDF, res.pca){
  xy <- as.data.frame(dataEllipse(ploMN0[,1],ploMN0[,2],levels = 0.95,draw = FALSE))
  xmax <- ifelse(max(xy[,1])>=max(ploMN0[,1]), max(xy[,1]), max(ploMN0[,1]))
  ymax <- ifelse(max(xy[,2])>=max(ploMN0[,2]), max(xy[,2]), max(ploMN0[,2]))
  
  par(mar=c(12, 6, 4, 3)+0.1,font.lab=2,font.axis=1,cex.lab=1.5,cex.axis=1,las=2,xpd=TRUE) #family="RMN",
  png(filename = paste0(pca.path, "/PCA.png"),width = 18,height = 14,units = "cm",res = 300)
  plot(ploMN0,main = "Scores (PCA)",
       xlab = paste0("t", "1", " (",as.character(round(res.pca[["eig"]][1,2],1)),"%)"),
       ylab = paste0("t", "2", " (",as.character(round(res.pca[["eig"]][2,2],1)),"%)"),
       xlim = c(-xmax,xmax),
       ylim = c(-ymax,ymax),
       col = colpalette[as.numeric(as.factor(samDF[,2]))],
       pch = as.numeric(as.factor(samDF[,2])) + ifelse(length(unique(as.numeric(as.factor(samDF[,2])))) > 6, 0, 14),
       cex = 1.2)
  abline(v = 0)
  abline(h = 0)
  dataEllipse(ploMN0[,1],ploMN0[,2],add=T,levels = 0.95,plot.points = F,col="black",lwd=1,center.pch=FALSE)
  legend("topright", legend = unique(as.character(as.factor(samDF[,2]))), col = unique(colpalette[as.numeric(as.factor(samDF[,2]))]), 
         pch = unique(as.numeric(as.factor(samDF[,2]))) + ifelse(length(unique(as.numeric(as.factor(samDF[,2])))) > 6, 0, 14))
  dev.off()
  pdf(file = paste0(pca.path, "/PCA.pdf"),width = 9,height = 7)
  plot(ploMN0,main = "Scores (PCA)",
       xlab = paste0("t", "1", " (",as.character(round(res.pca[["eig"]][1,2],1)),"%)"),
       ylab = paste0("t", "2", " (",as.character(round(res.pca[["eig"]][2,2],1)),"%)"),
       xlim = c(-xmax,xmax),
       ylim = c(-ymax,ymax),
       col = colpalette[as.numeric(as.factor(samDF[,2]))],
       pch = as.numeric(as.factor(samDF[,2])) + ifelse(length(unique(as.numeric(as.factor(samDF[,2])))) > 6, 0, 14),
       cex = 1.2)
  abline(v = 0)
  abline(h = 0)
  dataEllipse(ploMN0[,1],ploMN0[,2],add=T,levels = 0.95,plot.points = F,col="black",lwd=1,center.pch=FALSE)
  legend("topright", legend = unique(as.character(as.factor(samDF[,2]))), col = unique(colpalette[as.numeric(as.factor(samDF[,2]))]), 
         pch = unique(as.numeric(as.factor(samDF[,2]))) + ifelse(length(unique(as.numeric(as.factor(samDF[,2])))) > 6, 0, 14))
  dev.off()
  
  
  resultmax <- max(max(xy),xmax,ymax)
  png(filename = paste0(pca.path, "/PCA_samples.png"),width = 18,height = 14,units = "cm",res = 300)
  plot(ploMN0,main="Scores (PCA)",
       xlab = paste0("t", 1, " (",round(res.pca[["eig"]][1,2],1),"%)"),
       ylab = paste0("t", 2, " (",round(res.pca[["eig"]][2,2],1),"%)"),
       xlim = c(-resultmax,resultmax),
       ylim = c(-resultmax,resultmax),
       col = "white",
  )
  text(ploMN0[,1],ploMN0[,2],labels = rownames(ploMN0),col = colpalette[as.numeric(as.factor(samDF[,2]))],cex = 0.8)
  abline(v = 0)
  abline(h = 0)
  dataEllipse(ploMN0[,1],ploMN0[,2],add=T,levels = 0.95,plot.points = F,col="black",lwd=1,center.pch=FALSE)
  legend("topright", legend = unique(as.character(as.factor(samDF[,2]))), col = unique(colpalette[as.numeric(as.factor(samDF[,2]))]), 
         pch = unique(as.numeric(as.factor(samDF[,2]))) + ifelse(length(unique(as.numeric(as.factor(samDF[,2])))) > 6, 0, 14))
  dev.off()
  pdf(file = paste0(pca.path, "/PCA_samples.pdf"),width = 9,height = 7)
  plot(ploMN0,main="Scores (PCA)",
       xlab = paste0("t", 1, " (",round(res.pca[["eig"]][1,2],1),"%)"),
       ylab = paste0("t", 2, " (",round(res.pca[["eig"]][2,2],1),"%)"),
       xlim = c(-resultmax,resultmax),
       ylim = c(-resultmax,resultmax),
       col = "white",
  )
  text(ploMN0[,1],ploMN0[,2],labels = rownames(ploMN0),col = colpalette[as.numeric(as.factor(samDF[,2]))],cex = 0.8)
  abline(v = 0)
  abline(h = 0)
  dataEllipse(ploMN0[,1],ploMN0[,2],add=T,levels = 0.95,plot.points = F,col="black",lwd=1,center.pch=FALSE)
  legend("topright", legend = unique(as.character(as.factor(samDF[,2]))), col = unique(colpalette[as.numeric(as.factor(samDF[,2]))]), 
         pch = unique(as.numeric(as.factor(samDF[,2]))) + ifelse(length(unique(as.numeric(as.factor(samDF[,2])))) > 6, 0, 14))
  dev.off()
}

#包含QC的PCA图
data <- log2(temp.qc[,-1] + 1)
data <- na.omit(data)
data_t <- data.frame(t(data))
res.pca <- PCA(data_t, graph = FALSE)
ploMN0 <- as.matrix(res.pca[["ind"]]$coord)
samDF <- data.frame(name = rownames(data_t), group = samples[rownames(data_t), 3])
pca.path <- createdir(paste0(getwd(), '/PCA_QC'))
pca.plot(pca.path, ploMN0, samDF, res.pca)

#不包含QC的PCA图
temp <- temp.qc[,-grep("(^QC-\\d+$)|(^QC\\d+$)",colnames(temp.qc))]
data <- log2(temp[,-1] + 1)
data <- na.omit(data)
data_t <- data.frame(t(data))
res.pca <- PCA(data_t, graph = FALSE)
ploMN0 <- as.matrix(res.pca[["ind"]]$coord)
samDF <- data.frame(name = rownames(data_t), group = samples[rownames(data_t), 3])
pca.path <- createdir(paste0(getwd(), '/PCA'))
pca.plot(pca.path, ploMN0, samDF, res.pca)

