#Lecture 6 Live Coding Exercise

#Clear your environment
rm(list = ls())

#Load packages
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(FlowSorted.CordBlood.450k)

#Load data into minfi
#Get the 450k annotation data
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)
#Now we need to import the data
#Set up a path to your data directory
dataDirectory <-  "/Users/Christina/Desktop/Western Assistant Professor/Teaching/MBI4650F/IDATs/"
#Read in the sample sheet for the experiment
targets <- read.metharray.sheet(dataDirectory, pattern="sample_sheet.csv")
targets$Basename <- paste0("/Users/Christina/Desktop/Western Assistant Professor/Teaching/MBI4650F/IDATs/", targets$Sample_Name)
targets
#Ensure you have 22 samples

#Read in the raw data from the IDAT files
rgSet <- read.metharray.exp(targets=targets)
rgSet
#Give the samples descriptive names
targets$ID <- paste(targets$Group,targets$Sample,sep=".")
sampleNames(rgSet) <- targets$ID
rgSet

#Calculate the detection p-values
detP <- detectionP(rgSet)
head(detP)
#Remove poor quality samples
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]
rgSet
#Remove poor quality samples from detection p-value table
detP <- detP[,keep]
dim(detP)

#Normalize the data; this results in a GenomicRatioSet object
mSetSq <- preprocessQuantile(rgSet)

#Ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]
#Remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
table(keep)
mSetSqFlt <- mSetSq[keep,]
mSetSqFlt
dim(mSetSqFlt)

#Remove probes with SNPs at CpG site
#You can either remove all probes affected by SNPs (default), or only those with minor allele frequencies greater than a specified value. We will use the default for today
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
dim(mSetSqFlt)

#We will also filter out probes that have shown to be cross-reactive, that is, probes that have been demonstrated to map to multiple places in the genome. This list was originally published by @Chen2013 and can be obtained from the authors [website](http://www.sickkids.ca/MS-Office-Files/Research/Weksberg%20Lab/48639-non-specific-probes-Illumina450k.xlsx).

#Exclude cross reactive probes
xReactiveProbes <- read.csv(file=paste(dataDirectory,
                                       "ChenEtAlList.csv",
                                       sep="/"), stringsAsFactors=FALSE)
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]
dim(mSetSqFlt)

#Calculate methylation beta values
bVals <- getBeta(mSetSqFlt)
head(bVals[,1:5])

#Cell Type Composition
CC <- estimateCellCounts(rgSet, compositeCellType = "CordBlood",
                   processMethod = "auto", probeSelect = "auto",
                   cellTypes = c("Bcell", "CD4T", "CD8T", "Gran", "Mono", 
                                 "nRBC"),
                   referencePlatform = "IlluminaHumanMethylation450k",
                   returnAll = FALSE, meanPlot = FALSE, verbose = TRUE)

#We should run PCA now, but we learned about that last week so we will skip that step for today

#The biological question of interest for this particular dataset is to discover differentially methylated probes between preterm and term births. There are a lot of factors that you might imagine need to be taken into account when we perform the statistical analysis. For example, the weight of the child at birth may affect the outcome of methylation and therefore needs to be accounted for.

#To account for any required covariates we can run a regression model such that each methylation site is the dependent variable and the birth term (preterm or term) status is the independet variable, adjusted for the other relevant covariates.
  
#This is the factor of interest
levels(targets$Group) <- c("Preterm", "Term")
birth <- factor(targets$Group)
#These are some of the individual effects that we can account for, there are plenty of other factors we can/should include but we will focus on just these sex for today, keep in mind that cell types would be adjusted for here
sex_child <- as.factor(targets$Sex_Child)

#Try the formula out with my favourite CpG
summary(lm(bVals[grep("cg26094004", rownames(bVals)),] ~ birth +  sex_child))

#Running all CpGs would look like this in minfi (or we can raw code a loop)
design <- model.matrix(~0+birth+sex_child, data=targets)
colnames(design) <- c("Preterm", "Term", "sex_childmale")

#Fit the actual linear model to the data
fit <- lmFit(bVals, design)

#Create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(Term-Preterm,
                            levels=design)
contMatrix

#Fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
#Rank genes
fit2 <- eBayes(fit2)
#Get the table of results
DMPs <- topTable(fit2, num=Inf, coef=1)
head(DMPs)

## look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))

#qqplot of the results
pvals <- DMPs$P.Value
observed <- sort(pvals)
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
m="qqplot"
plot(c(0,7), c(0,7), col="red", lwd=4, type="l", xlab="Expected (-logP)",
     ylab="Observed (-logP)", xlim=c(0,7), ylim=c(0,7), las=1, xaxs="i", yaxs="i", bty="l", main=m)
points(lexp, lobs, pch=23, cex=.5, col="black", bg="black")

#Calculate lambda
library(QCEWAS)
P_lambda(pvals)

#Adjust the pvalues, convert them to FDRs
#Because we are testing thousands of hypotheses, we need to correct our p-values for this multiplicity of tests, **p.adjust()** returns adjusted p-values for different methods:

pvals_fdr <- p.adjust(pvals, method = "fdr")
summary(pvals_fdr)
head(sort(pvals_fdr))

#END