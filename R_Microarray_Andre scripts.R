##################################################
# R code for analysis of Andre's microarray data #
# A. P. Levine, 02/11/2014			  			 #
##################################################


#Lumi
#http://www.bioconductor.org/packages/release/bioc/vignettes/lumi/inst/doc/lumi.pdf

#Limma, page 40 for differential
#http://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf

.libPaths(c('C:/Rlibrary',.libPaths()))

#Install Bioconductor packages, only do this the very first time on your computer
source("http://bioconductor.org/biocLite.R")
biocLite("lumi")
biocLite("lumiHumanIDMapping")
biocLite("limma")
biocLite("illuminaHumanv4.db")

#Load packages, do every time
library(lumi)
library(lumiHumanIDMapping)
library(limma)
library(illuminaHumanv4.db)

#Set working directory
#######CHANGE THIS TO THE DIRECTORY WHERE THE DATA ARE SAVED
#setwd("E:/Documents/APLevine/THP1_benzoate")


#Load data
files <- c('SmithSampleProbe301014.txt')
#expr = lumiR.batch(files,lib.mapping='lumiHumanIDMapping')
expr = lumiR.batch(files)

##Add control probes
#controlData <- getControlData('SmithControlProbe301014.txt')
#plotControlData(controlData, type='NEGATIVE')
#expr.c <- addControlData2lumi(controlData, expr)

#log2 transform
exprT = lumiT(expr,method="log2")

#Plot distributions
plot(exprT)

#Quantile normalisation
exprTN <- lumiN(exprT,method='quantile')

#Plot distributions
plot(exprTN)

#Plot relations
plot(exprTN, what='sampleRelation')

#Extract data
exprTNd <- exprs(exprTN)

#Detection summary
thr = detectionCall(exprTN, Th = 0.01, type = c('probe'))
barplot(table(thr))
keep <- which(thr>1)
exprTNdk = exprTNd[keep,]

#Run PCA
exprTNdk.center = exprTNdk - apply(exprTNdk,1,mean)
pc = prcomp(exprTNdk.center)

#Find group that each sample belongs to
id <- row.names(pc$rotation)
group=c()
for (i in id){
	this=substr(i,1,nchar(i)-1)
	group=c(group,this)
}

#Give each group a colour
cols = rep(NA,length(group))
count=0
for (i in unique(group)){
	count=count+1
	cols[which(group==i)] <- count
}

#Plot PCA
plot(pc$rotation[,1],pc$rotation[,2],col=cols,xlab="PC1",ylab="PC2",cex=2)
legend(x="topleft",legend=unique(group),col=1:length(unique(group)),pch=1,cex=1.5)

#Gene names
ls("package:illuminaHumanv4.db")
genes.symbol <- as.data.frame(illuminaHumanv4SYMBOL)

#Differential analysis
#Subset 'Un' and 'NaB'
test.groups <- which(group %in% c("Un","NaB"))
test.data <- exprTNdk[,test.groups]
groups.groups <- group[test.groups]
#Prepare design layout
left <- row.names(as.data.frame(test.data[1,]))
g1 = rep(0,length(left))
g1[which(groups.groups=="Un")] <- 1
g2 = rep(0,length(left))
g2[which(groups.groups=="NaB")] <- 1
design <- cbind(Un=g1,NaB=g2)
row.names(design) <- left

eset <- test.data


#Add gene names
top.g = merge(top,genes.symbol,by.x="ID",by.y="probe_id")
top = top.g[order(top.g$adj.P.Val),]

write.csv(top,file="differential.csv",row.names=F)

