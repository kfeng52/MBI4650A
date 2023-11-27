#Lecture 4 Live Coding Exercise

#Clear your environment
rm(list = ls())

#Load packages
library(factoextra)
library(devtools)
library(Biobase)
library(AnnotationDbi)
library(stringr)
library(ggfortify)

#We will use the methylation set from lecture 3
mdata <- read.csv("/Users/Christina/Desktop/IDATfiles/GSE68777/idat/GSE68777_RAW/Beta.csv", row.names=1)
pdata <- read.csv("/Users/Christina/Desktop/IDATfiles/GSE68777/idat/GSE68777_RAW/Pheno.csv", row.names=1)

#We may wish to pull out the variables of slide and well ID from the sample name
pdata
pdata2 <- str_split_fixed(pdata$title, '_', 2)
colnames(pdata2) <- c("slide", "well")
pdata <- cbind(pdata, pdata2)

#Confirm correct classes
pdata$sex <- as.factor(pdata$sex)
pdata$group <- as.factor(pdata$group)
pdata$slide <- as.factor(pdata$slide)
pdata$well <- as.factor(pdata$well)
str(pdata)

#Run PCA on your dataframe of interest (all CpGs/Samples), here is some sample code. 
#In mDataNAT, you want rows to be samples and columns to be CpGs
#We also must remove all NAs prior to running PCA
mdataNA <- na.omit(mdata) 
mDataNAT <- t(mdataNA) 
#Check that rows are samples and columns are features*
pcaMdata <- prcomp(mDataNAT, center=TRUE, scale=FALSE) 
attributes(pcaMdata) 
summary(pcaMdata) 
PCs = pcaMdata$x 
head(PCs)

#Plot the PCA results 

#Visualize eigenvalues (scree plot). Show the percentage of variances explained by each principal component.
fviz_eig(pcaMdata)

#Make the same plot by hand and without the histogram
SumPCA <- summary(pcaMdata)
#Pull out the proportion of variance, row 2 of the table, for the first 10 PCs
plot(SumPCA$importance[2,1:10]*100, ylab="% Variance Explained", xlab="PC", type="b")

#Pull Out First 5 PCs
PC1<-pcaMdata$x[,1]
PC2<-pcaMdata$x[,2]
PC3<-pcaMdata$x[,3]
PC4<-pcaMdata$x[,4]
PC5<-pcaMdata$x[,5]

#Visualize PC1 versus PC2
#autoplot in package ggfortify will do this but a ggplot scatter plot works well also
autoplot(pcaMdata, data = pdata, colour = 'sex')
#Which PC is picking up sex?
autoplot(pcaMdata, data = pdata, colour = 'group')
autoplot(pcaMdata, data = pdata, colour = 'slide')
autoplot(pcaMdata, data = pdata, colour = 'well')
#if there is no clustering it doesn't mean that things are not related, it just means that PC1 and 2 are not picking up this variance

#Colour by phenotypes using ggplot
ggplot(PCs, aes(y=PC2, x=PC1, colour=pdata$sex)) + geom_point(shape=19, size=1) + stat_ellipse()
ggplot(PCs, aes(y=PC4, x=PC3, colour=pdata$group)) + geom_point(shape=19, size=1) + stat_ellipse()
ggplot(PCs, aes(y=PC4, x=PC3, colour=pdata$well)) + geom_point(shape=19, size=1)
ggplot(PCs, aes(y=PC4, x=PC3, colour=pdata$slide)) + geom_point(shape=19, size=1)

#Look for PC outliers
densityplot(PC1, pch=19)
densityplot(PC2, pch=19)
densityplot(PC2, pch=19, col=pdata$sex)
densityplot(PC3, pch=19)

#Drop PC outliers
sample_outliers=c()
alloutliers=c()
for(i in 1:10){
  a<-subset(rownames(pcaMdata$x), pcaMdata$x[,i] > (mean(pcaMdata$x[,i])+3*sd(pcaMdata$x[,i])))
  b<-subset(rownames(pcaMdata$x), pcaMdata$x[,i] < (mean(pcaMdata$x[,i])-3*sd(pcaMdata$x[,i])))
  out<-c(a,b)
  sample_outliers <- c(sample_outliers,out)
  print(paste("outliers in PCA",i,":",sep=""))
  print(sample_outliers)
  alloutliers=c(alloutliers,sample_outliers)
  sample_outliers=c()
}	
outlier<-unique(alloutliers)
outlier

densityplot(PC5, pch=19)

#You may wish to remove this sample from your data before looking at the correlation matrix
#We are going to leave all samples in for now


# Make a correlation matrix
#There are many different ways to do this, here is one method

#First we will merge the PCs with our pdata
pdata$title <- sub("^", "X", pdata$title )
pdataPCs <- merge(pdata, PCs, by.x="title", by.y="row.names")

#Plot correlation matrix, we will also send this to a file since it is so large
setwd("/Users/Christina/Desktop/IDATfiles/GSE68777/idat/GSE68777_RAW/")
pdf("modTraitCor_Group.pdf")
twolines = function(x,y) {
  points(x,y,pch=20, col=pdataPCs$group)
  abline(lm(y~x),col="red")
  #legend("bottomright", paste("R=",prettyNum(cor(x,y), digits=3)),bty="n" ,cex=1.5)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use="complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
mydiag.panel <- function( x,  labels, ...){
  ll <- par("usr")
  rect(ll[1], ll[3], ll[2], ll[4], col="darkolivegreen1")
}

diag.labels=c("group", "sex", "slide", "well", "PC1", "PC2","PC3","PC4")

plot.formula=as.formula(~group+sex+slide+well+PC1+PC2+PC3+PC4)

pairs(plot.formula, data=pdataPCs, upper.panel=twolines, labels=diag.labels, diag.panel=mydiag.panel, lower.panel=panel.cor, label.pos=0.5, main="Correlation Between Variables")
dev.off()

#We notice that slide and group are correlated

# TOPIC 5: Heatmaps
betast <- as.data.frame(t(mdata))
betas_pheno <- merge(pdata, betast, by.x="title", by.y="row.names")
row.names(betas_pheno) <- betas_pheno$title

#A common type of plot for 'omics data is a heatmap. They are usually used for visualizing matrices. For example we can look at all CpGs that meet a certain variance threshold: 
betas_phenot <- as.data.frame(t(betas_pheno))
betas_phenot$var <- apply(betas_phenot, 1, var)
#choose a reasonable cutoff
betasHighVar <- subset(betas_phenot, betas_phenot$var > 0.04)

#2384

betasHighVar$var <- NULL
betasHighVar$age <- NULL
ematrix = data.matrix(betasHighVar, rownames.force = NA)

pdf("heatmap.pdf")
heatmap(ematrix)
dev.off()
#We can change the coloring since this one is a little hard to see. To do this you have to set up a color palette. 
colramp = colorRampPalette(c(3,"white",2))(9)
heatmap(ematrix,col=colramp)

#You might have noticed some automatic clustering here, you can turn that off 
heatmap(ematrix,col=colramp,Rowv=NA,Colv=NA)

#If you load the `gplots` package you can add a color scale with the `heatmap.2` package. Here we have to add some options to make one dendogram disappear, scale the data by rows, and remove a tracing plot

pdf("heatmap2.pdf")
par(mar=c(5.1,4.1,4.1,2.1))
heatmap.2(ematrix,col=colramp,Rowv=NA,
          dendrogram="column", scale="row",trace="none")
dev.off()

#What is this clustering picking up? Try group, then try sex
colnames(ematrix) <- betas_pheno$sex
heatmap(ematrix)

#END