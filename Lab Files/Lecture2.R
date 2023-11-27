#Lecture 2 Live Coding Exercise

#To put files online for them, first transfer into my folder then sudo transfer them into /home/shared/
#Read in a sample methylation dataset for data exploration
rm(list = ls())
library(lattice)
library(ggplot2)
library(gplots)
#cd /home/ccastel/MBI4650F
betas <- read.csv("/Users/Christina/Desktop/IDATfiles/GSE68777/idat/GSE68777_RAW/Beta.csv", row.names=1)
pdata <- read.csv("/Users/Christina/Desktop/IDATfiles/GSE68777/idat/GSE68777_RAW/Pheno.csv", row.names=1)

# TOPIC 1: Merging

#Find a common barcode, in this case the row.names of betas and the column called title in pdata
#But wait, they don't actually match

#Let's add an X before the ID in pdata$title
#The regex pattern "^" (outside any character-class brackets) represents the point just before the first character of a "character"-class item
pdata$title <- sub("^", "X", pdata$title )

#The beta table needs to be transposed for this and other tasks later
betast <- t(betas)
betast <- as.data.frame(betast)

#Now we can merge the two files
betas_pheno <- merge(pdata, betast, by.x="title", by.y="row.names")
row.names(betas_pheno) <- betas_pheno$title
#Check that it merged correctly

#Why is merging important?

# TOPIC 2: Exploration

##Look at the properties of betas and pdata

##Dimensions first, what are we making sure of?
#1. Numbers are as expected
#2. The files "match up"
dim(betas)
dim(pdata)
betas[1:5,1:5]
pdata[1:5,]

#Quick tables
table(pdata$sex)
#Cross tables, For example sex by group
table(pdata$group,pdata$sex)
#What do you notice? Does anything look unbalanced?
myTable <- table(pdata$sex)
myTable
prop.table(myTable)
#Can also do two-way tables using this format
myTable2 <- table(pdata$group, pdata$sex)
myTable2
prop.table(myTable2)

#Check a summary of the distribution to look for scale, this is also one way to check for `NA` values
summary(betas)
summary(pdata)
#str is useful for looking at classes of data
str(betas)
str(pdata)

#The class of group and sex is incorrect, for our statistical analysis we want these to be factors not character vetors
pdata$sex <- as.factor(pdata$sex)
pdata$group <- as.factor(pdata$group)
str(pdata)

#Look for missing values
#`NA` is the most common character for missing values, but sometimes they are coded as spaces, 999, -1 or "missing". Check for missing values in a variety of ways
#Use option useNA to include NA's in table
table(pdata$sex,useNA="ifany")
pdata$group
#is.na checks for NA values
table(is.na(pdata$group))

#Check genomic data for NAs
sum(is.na(betas))

# TOPIC 3: QC Plots
#Density plot to check for outliers
#Let's start with a boxplot
boxplot(betas,col=2,range=0) #whiskers extend to data extremes

#We can also look at this sample by sample with histograms
par(mfrow=c(1,2))
hist(betas[,1],col=2)
hist(betas[,2],col=2)
#What do you notice?

setwd("/Users/Christina/Desktop/IDATFiles/")
#Or with density plots
pdf("Example1.pdf")
densityplot(betas[,1], pch=19)
dev.off()
#What does the distribution look like by sample?
dev.off()
plot(density(betas[,1]),col=2)
lines(density(betas[,2]),col=3)

#density plot with ggplot
ggplot(betas,aes(x=X5958091020_R02C02)) + geom_density(alpha=0.25) + geom_rug() + xlab("X5958091020_R02C02")

##This is a good time to check for obvious outliers, do we have any?

#Density plot of one CpG coloured by sex
ggplot(betas_pheno,aes(x=cg18478105, color=as.factor(sex), fill=as.factor(sex))) + geom_density(alpha=0.25) + geom_rug() + xlab("cg18478105")
#What do you conclude?

#Density plot of one CpG coloured by sample group
ggplot(betas_pheno,aes(x=cg18478105, color=as.factor(group), fill=as.factor(group))) + geom_density(alpha=0.25) + geom_rug() + xlab("cg18478105")
#What do you conclude?

#Keep in mind this is only 1 CpG
 
# TOPIC 4: Subsetting, correlations and getting setup for regressions

#What parameters can we subset data based on?
#Subset to just include males

BP_MOnly <- subset(betas_pheno, betas_pheno$sex == "Male")
dim(BP_MOnly)

#Subset to only include samples with an age above 20
betas_pheno$age <- sample(10:70, nrow(betas_pheno), replace=FALSE)
densityplot(betas_pheno$age)
BP_OldOnly <- subset(betas_pheno, betas_pheno$age > 20)
dim(BP_OldOnly)

#simple correlations
cor(betas_pheno$age,  betas_pheno$cg20826792)
cor(betas_pheno$age,  betas_pheno$cg20826792, method="spearman")

#simple plot of the same data, think about who is dependent and who is independent variable
plot(betas_pheno$age,  betas_pheno$cg20826792, pch=19)
abline(lm(betas_pheno$cg20826792~betas_pheno$age), col=3)

# TOPIC 5: Heatmaps

#A common type of plot for 'omics data is a heatmap. They are usually used for visualizing matrices. For example we can look at all CpGs that meet a certain variance threshold: 
betas_phenot <- as.data.frame(t(betas_pheno))
betas_phenot$var <- apply(betas_phenot, 1, var)
#look at the plot to determine something reasonable as a cutoff
betasHighVar <- subset(betas_phenot, betas_phenot$var > 0.04)
#2384
#colnames(betasHighVar) <- betas_pheno$group
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

#END
