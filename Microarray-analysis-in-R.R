#Microarray analysis in R

BiocManager::install("AnnotationDbi")
BiocManager::valid()
BiocManager::version()
BiocManager::available()
BiocManager::install("GEOquery")

#download processed data
library(GEOquery)

eList<- getGEO("GSE11675")
class(eList)
length(eList)
names(eList)
eData<-eList[[1]]
eData

#download raw data
elist2<-getGEOSuppFiles("GSE11675")
elist2
class(elist2)
rownames(elist2)<-basename(rownames(elist2))
elist2
tarArchive<-rownames(elist2[1])
tarArchive

#downloading data
library(GEOquery)

getGEOSuppFiles("GSE38792")
list.files("GSE38792")
untar("GSE38792/GSE38792_RAW.tar", exdir = "GSE38792/CEL")
list.files("GSE38792/CEL")

#reading raw data
BiocManager::install("oligo")
library(oligo)
celfiles<-list.files("GSE38792/CEL", full=TRUE)
celfiles
rawData<-read.celfiles(celfiles)
rawData

#genefeatureset object
getClass("GeneFeatureSet")
exprs(rawData[1:4,1:3])

#the unit of measurement is integer
#on a 16 bit scanner
2^16
max(exprs(rawData))

#clean up the phenotype information
filename<-sampleNames(rawData)
filename
pData(rawData)$filename<-filename
sampleNames<-sub(".*_", "", filename)
sampleNames<-sub(".CEL.gz$", "", sampleNames)
sampleNames(rawData)<- sampleNames
pData(rawData)$group<-ifelse(grepl("^OSA", sampleNames(rawData)), "OSA", "Control")
pData(rawData)

#normalization
boxplot(rawData, target="core")

normData<-rma(rawData)
normData
exprs(normData)
boxplot(normData)

#expressionset object
BiocManager::install("ALL")
library(ALL)
data(ALL)
ALL
library(Biobase)
experimentData(ALL)
str(exprs(ALL)) 
class(exprs(ALL))
exprs(ALL[1:4,1:4])
sampleNames(ALL)
head(sampleNames(ALL))
head(featureNames(ALL))
head(pData(ALL))
head(pData(ALL)$sex)
head(ALL$sex)
ALL$sex
ALL[,1:5]
ALL[1:10,]
ALL[,c(3,2,1)]
ALL$sex[c(1,2,3)]
ALL[,c(1,2,3)]$sex

#featuredata
featureData(ALL)

#annotation
ids<-featureNames(ALL)[1:5]
ids
BiocManager::install("hgu95av2.db")
library(hgu95av2.db)
hgu95av2ENTREZID
ids
as.list(hgu95av2ENTREZID[ids])

#load the data
library(leukemiasEset)
data("leukemiasEset")
leukemiasEset
leukemiasEset$LeukemiaType
table(leukemiasEset$LeukemiaType)

#subset
ourData<-leukemiasEset[,leukemiasEset$LeukemiaType%in% c("ALL", "NoL")]
ourData$LeukemiaType<-factor(ourData$LeukemiaType)
ourData$LeukemiaType

#limma
library(limma)
design <- model.matrix(~ ourData$LeukemiaType)
fit<-lmFit(ourData, design)
fit<-eBayes(fit)
topTable(fit)

#level
ourData$LeukemiaType

#Fc by hand
topTable(fit, n=1)
genename<-rownames(topTable(fit, n=1))
genename
exprs(ourData)[genename,]
typemean<- tapply(exprs(ourData)[genename,], ourData$LeukemiaType, mean)
typemean
typemean["NoL"] - typemean["ALL"]

#design2
design2<-model.matrix(~ourData$LeukemiaType - 1)
design2
colnames(design2)
class(design2)
colnames(design2)<-c("ALL", "NoL")
fit2<-lmFit(ourData,design2)
fit2
contrast.matrix<-makeContrasts("ALL-NoL", levels=design2)
contrast.matrix
fit2c<-contrasts.fit(fit2, contrast.matrix)
fit2c<-eBayes(fit2c)
topTable(fit2c)
browseVignettes("limma")
tbl<-topTable(fit2c)
tbl
tbl<-topTable(fit2c, n=200)
tbl
write.table(tbl, file="R/DEGS.txt", quote=FALSE, sep= "\t")
list.files("R")
system("open R/DEGS.txt")
write.csv

#exploitary data analysis
#PCA
library(leukemiasEset)
data("leukemiasEset")
leukemiasEset
gset<-leukemiasEset
gset
exprs(gset)
t(exprs(gset))
PCA<-prcomp(t(exprs(gset)), scale=FALSE)
PCA
percentvar<-round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio<-sqrt(percentvar[2]/percentvar[1])
dataGG<-data.frame(PC1=PCA$x[,1], PC2=PCA$x[,2], phenotype=pData(gset)$LeukemiaType)
dataGG
library(ggplot2)
g<-ggplot(dataGG, aes(PC1,PC2))+geom_point(aes(shape=phenotype, colour=phenotype))+ ggtitle("PCA Plot")+ xlab(paste0("PC1,VarExp:",percentvar[1],"%"))+ ylab(paste0("PC2,VarExp:",percentvar[2],"%"))+ theme(plot.title = element_text(hjust = 0.5))+ coord_fixed(ratio = sd_ratio)
g
ggsave(filename="R/PCA.png", g, width=180,height=100, units="mm", dpi=300)


#heatmap
exp<-exprs(gset)
exp
annot<-data.frame(type=gset$LeukemiaType,row.names = sampleNames(gset))
dists<-as.matrix(dist(t(exp), method = "manhattan"))
row.names(dists)<-row.names(pData(gset))

library(RColorBrewer)
cols<- sort(colorRampPalette(brewer.pal(9,"YlOrRd"))(255))
colnames(dists)<-NULL
diag(dists)<- NA
annot_colors<-list(phenotype=c(ALL="charteuse4", AML="burlywood3", CLL="blue4", CML="cadetblue2", NoL="snow3"))

#pheatmap
library(pheatmap)
pheatmap(dists, col=cols, annotation_row = annot, annotation_colors = annot_colors, legend = TRUE, treeheight_row = 0, legend_breaks = c(min(dists, na.rm = TRUE), max (dists, na.rm = TRUE) ), legend_labels = (c("small distance", "large distance")), main = "LeukemiasEset")

#heatmap on genes with the lowest p.val
eset<-leukemiasEset[,leukemiasEset$LeukemiaType%in% c("ALL", "Nol",)]
eset$LeukemiaType
eset$LeukemiaType<-factor(eset$LeukemiaType, levels=c("ALL", "NoL"))
eset$LeukemiaType


#limma
design<-model.matrix(~eset$LeukemiaType)
library(limma)
fit<-lmFit(eset, design)
fit<- eBayes(fit)
tbl<-topTable(fit, number=200)
ind<- rownames(tbl)
mat<-exprs(gset[ind,])

library(pheatmap)
pheatmap(mat)

#colors
cols<-colorRampPalette(c("green", "black", "red"))(n=255)
cols
pheatmap(mat=mat, color=cols)

#clustering the columns
pheatmap(mat=mat, color=cols, cluster_cols=FALSE)

#clustering distance rows
pheatmap(mat=mat, color=cols, cluster_cols=FALSE, clustering_distance_rows="manhattan")

#cutree columns
pheatmap(mat=mat, color=cols, cutree_cols=5)

#annotation col
pheatmap(mat=mat, color=cols, cutree_cols=5, annotation_col=annot)

#show row names
pheatmap(mat=mat, color=cols, cutree_cols=5, annotation_col=annot, show_rownames=FALSE)

#fontsize of columns
pheatmap(mat=mat, color=cols, cutree_cols=5, annotation_col=annot, show_rownames=FALSE, fontsize_col=5)

#rename
colnames(mat)<-sub("\\.CEL$","",colnames(mat))
pheatmap(mat=mat, color=cols, cutree_cols=5, annotation_col=annot, show_rownames=FALSE, fontsize_col=5)

#angle
pheatmap(mat=mat, color=cols, cutree_cols=5, annotation_col=annot, show_rownames=FALSE, fontsize_col=5, angle_col="45")

#main
pheatmap(mat=mat, color=cols, cutree_cols=5, annotation_col=annot, show_rownames=FALSE, fontsize_col=5, angle_col="45", main="leukemiasEset")

#save
pheatmap(mat=mat, color=cols, cutree_cols=5, annotation_col=annot, show_rownames=FALSE, fontsize_col=5, angle_col="45", main="leukemiasEset", filename= "R/heatmap.png")
dev.off()

#enhanced volcano
#DEGs
library(leukemiasEset)
eset<-leukemiasEset[,leukemiasEset$LeukemiaType%in% c("ALL", "Nol",)]
eset$LeukemiaType
eset$LeukemiaType<-factor(eset$LeukemiaType, levels=c("ALL", "NoL"))
eset$LeukemiaType
design <- model.matrix(~ ourData$LeukemiaType)
#limma
library(limma)

fit<-lmFit(eset, design)
fit<-eBayes(fit)
tbl<-topTable(fit, n=2000)

library(EnhancedVolcano)
EnhancedVolcano(toptable=tbl, lab=rownames(tbl), x="logFC",y="adj.P.Val")

#title
EnhancedVolcano(toptable=tbl, lab=rownames(tbl), x="logFC",y="adj.P.Val", title="leukemiasEset")

#pCutOff
EnhancedVolcano(toptable=tbl, lab=rownames(tbl), x="logFC",y="adj.P.Val", title="leukemiasEset", pCutoff=10e-8)

#FCcutoff
EnhancedVolcano(toptable=tbl, lab=rownames(tbl), x="logFC",y="adj.P.Val", title="leukemiasEset", pCutoff=10e-8, FCcutoff=4)

#pointsize
EnhancedVolcano(toptable=tbl, lab=rownames(tbl), x="logFC",y="adj.P.Val", title="leukemiasEset", pCutoff=10e-8, FCcutoff=4, pointSize=0.6)



