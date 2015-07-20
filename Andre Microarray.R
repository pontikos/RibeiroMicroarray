##################################################
# R code for analysis of Andre's microarray data #
# A. P. Levine, 02/11/2014			  			 #
##################################################

load ("Microarray results")

#Lumi
#http://www.bioconductor.org/packages/release/bioc/vignettes/lumi/inst/doc/lumi.pdf

#Limma, page 40 for differential
#http://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf

source("http://bioconductor.org/biocLite.R")
biocLite("ALL")
updateR()
#biocLite("BiocUpgrade")
biocLite()


#Install Bioconductor packages, only do this the very first time on your computer
#source("http://bioconductor.org/biocLite.R")
#biocLite("lumi")
#biocLite("lumiHumanIDMapping")
#biocLite("limma")
#biocLite("illuminaHumanv4.db")

biocLite("GenomicRanges")
biocLite("GenomeInfoDb")
biocLite("RSQLite")
biocLite("registry")
biocLite("checkmate")
biocLite("locfit")
biocLite("beanplot")
biocLite("preprocessCore")
biocLite("methylumi")
biocLite("XML")
biocLite("rb")
biocLite("")
biocLite("")




#Load packages, do every time
library(lumi)
library(lumiHumanIDMapping)
library(limma)
library(illuminaHumanv4.db)

#Set working directory
#######CHANGE THIS TO THE DIRECTORY WHERE THE DATA ARE SAVED
#setwd("E:/Documents/APLevine/THP1_benzoate")

setwd("//ad.ucl.ac.uk/slms/home3/rmhvlrr/Desktop/Microarray results")

#Load data
files <- c('SmithSampleProbe301014.txt')
expr = lumiR.batch(files,lib.mapping='lumiHumanIDMapping')
expr = lumiR.batch(files)

#Andre
#expr
#files
#lumiR.batch(files)
read.table("SmithSampleProbe301014.txt")
d = read.table("SmithSampleProbe301014.txt", 
               sep="\t", 
               col.names=c("all"), 
               fill=FALSE, 
               strip.white=TRUE)

write.csv(expr,file="non.normalization.csv",row.names=F)


##Add control probes
controlData <- getControlData('SmithControlProbe301014.txt')
plotControlData(controlData, type='NEGATIVE')
expr.c <- addControlData2lumi(controlData, expr)

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

#exprTNd

#head(exprTNdk)

#Gene names
ls("package:illuminaHumanv4.db")
genes.symbol <- as.data.frame(illuminaHumanv4SYMBOL)


#exprTNdk
exprTNdk = as.data.frame(exprTNd)
exprTNdk$probe_id = row.names(exprTNdk)
out = merge(exprTNdk,genes.symbol,by.x="probe_id",by.y="probe_id")
write.csv(exprTNdk,file="post_normalisation3.csv",row.names=F)

#write.csv(exprTNdk,file="post_normalisation.csv",row.names=T)
#write.csv(exprTNd ,file="normalization.csv",row.names=F)

exprTN
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

#Differential analysis UN x NaB
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

fit <- lmFit(eset, design)
cont.matrix <- makeContrasts(UnVNaB=Un-NaB, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
top <- topTable(fit2, adjust="BH", p.value=0.05,number=500)
top$ID = row.names(top)

#Add gene names
top.g = merge(top,genes.symbol,by.x="ID",by.y="probe_id")
top = top.g[order(top.g$adj.P.Val),]

write.csv(top,file="differential UNxNaB.csv",row.names=F)


#Differential analysis EC x ENaB
#Subset 'EC' and 'ENaB'
test.groups <- which(group %in% c("EC","ENaB"))
test.groups
test.data <- exprTNdk[,test.groups]
groups.groups <- group[test.groups]
#Prepare design layout
left <- row.names(as.data.frame(test.data[1,]))
g1 = rep(0,length(left))
g1[which(groups.groups=="EC")] <- 1
g2 = rep(0,length(left))
g2[which(groups.groups=="ENaB")] <- 1
design <- cbind(EC=g1,ENaB=g2)
row.names(design) <- left
design
eset <- test.data
fit <- lmFit(eset, design)
cont.matrix <- makeContrasts(ECVENaB=EC-ENaB, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
top <- topTable(fit2, adjust="BH", p.value=0.05,number=500)
top$ID = row.names(top)
top
#Add gene names
top.g = merge(top,genes.symbol,by.x="ID",by.y="probe_id")
top = top.g[order(top.g$adj.P.Val),]

write.csv(top,file="differential ECxENaB.csv",row.names=F)

#Differential analysis PG x PGNaB
#Subset 'PG' and 'PGNaB'
test.groups <- which(group %in% c("PG","PGNaB"))
test.groups
test.data <- exprTNdk[,test.groups]
groups.groups <- group[test.groups]
#Prepare design layout
left <- row.names(as.data.frame(test.data[1,]))
g1 = rep(0,length(left))
g1[which(groups.groups=="PG")] <- 1
g2 = rep(0,length(left))
g2[which(groups.groups=="PGNaB")] <- 1
design <- cbind(PG=g1,PGNaB=g2)
row.names(design) <- left
design
eset <- test.data
fit <- lmFit(eset, design)
cont.matrix <- makeContrasts(PGVPGNaB=PG-PGNaB, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
top <- topTable(fit2, adjust="BH", p.value=0.05,number=500)
top$ID = row.names(top)
top
#Add gene names
top.g = merge(top,genes.symbol,by.x="ID",by.y="probe_id")
top = top.g[order(top.g$adj.P.Val),]

write.csv(top,file="differential PGxPGNaB.csv",row.names=F)

save.image ("Microarray results")

