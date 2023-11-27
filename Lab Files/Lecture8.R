#Lecture 8 Live Coding Exercise

#Today we will work with matched RNA sequencing and methylation array (EPIC) data. First we will use the EPIC data to call DMRs using the package DMRcate, then we will look at how to integrate two types of "omics" data. As we mentioned in class, there are many different ways to integrate multiomics and a long list of considerations for the integration of multiple datatypes. We will explore one available method using the packacge Elmer in R.

#Clear your environment
rm(list = ls())

#Load libraries
library(ELMER)
library(data.table)
library(dplyr)
library(biomaRt)
library(matrixStats)
library(sesameData)
library(DMRcate)

setwd("/home/shared/Epigenomics/")
Meth <- read.csv("Meth.csv", row.names=1)
RNA <- read.csv("RNA.csv", row.names=1)

dim(Meth)
dim(RNA)

#Let's talk about what the data represents
Meth[1:5,1:5]
RNA[1:5,1:5]

########################
########################
######DMR Calling#######
#######################
#######################

type <- factor(c("KO", "NC", "KO", "NC", "KO", "NC", "NC", "NC", "NC", "KO", "KO", "KO"))
design <- model.matrix(~type)

#coef tells the program which column of model matrix to use
Meth <- as.matrix(Meth)
myannotation <- cpg.annotate("array", Meth, arraytype="EPIC", analysis.type="differential", design=design, coef=2, what="Beta", fdr=0.01)

#C=scaling factor for bandwidth (usually 2), lambda = 1000 nucleotides
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2, min.cpgs=10)
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
results <- as.data.frame(results.ranges)

head(results)

groups <- c(KO="magenta", NC="forestgreen")
cols <- groups[as.character(type)]
cols

#Visualization
#This command is currently not working as it is a known issue with DMRcate at the moment
pdf('DMR_1.pdf', height=12, width=11)
DMR.plot(ranges=results.ranges, dmr=1, CpGs=Meth, what="Beta",
arraytype = "EPIC", phen.col=cols, genome="hg19")
dev.off()

#######################
#######################
#ELMER Omics Integration#
#######################
#######################

#A quick crash course on count measures in gene expression data
#https://www.reneshbedre.com/blog/expression_units.html

#Create a dataframe to specify sample group Experiment (KO), Control (NC)
primary <- colnames(RNA)
GroupLabel <- c("Experiment", "Control", "Experiment", "Control", "Experiment", "Control", "Control", "Control", "Control", "Experiment", "Experiment", "Experiment")
Sample <- cbind(primary, GroupLabel)
Sample <- as.data.frame(Sample)
row.names(Sample) <- Sample$primary
Sample

#For more information about data input in ELMER, see: https://www.bioconductor.org/packages/release/bioc/vignettes/ELMER/inst/doc/input.html#Input_data

#A Multi Assay Experiment object from the MultiAssayExperiment package is the input for multiple main functions of ELMER.

#Create MAE object for ELMER
data <- createMAE(exp = RNA, 
                  met = Meth,
                  met.platform = "EPIC",
                  genome = "hg19",
                  save = FALSE,
                  TCGA = FALSE, 
                  colData=Sample,
                  met.na.cut=0.05
)

#exp	 (RNA): a matrix or path of file containing the expression data. Rownames should be either Ensembl gene id (ensembl_gene_id) or gene symbol (external_gene_name)
#met	(Meth): a matrix or path of file containing only the methylation data
#colData (Sample): a DataFrame or data.frame of the phenotype data for all participants. Must have column primary (sample ID).
#met.na.cut (0.05): defines the percentage of NA that the line should have to remove the probes for human methylation platforms.
#met.platform (EPIC): DNA methylation platform "450K" or "EPIC"
#genome (hg19): Which is the default genome to make gene information. Options hg19 and hg38

#Find the nearest 20 genes (10 upstream, 10 downstream) of significant CpGs
#To make this run faster we will select 5 instead
nearGenes <- GetNearGenes(data = data, 
                          probes = row.names(Meth), 
                          numFlankingGenes = 5)
                          
#save(nearGenes, file="nearGenes.rda")
load(file="nearGenes.rda")

#This takes a long time to run so we will not run it today
pairs <- suppressMessages(get.pair(data = data,
                  group.col = "GroupLabel",
                  group1 =  "Control",
                  group2 = "Experiment",
                  nearGenes = nearGenes,
                  mode = "supervised",
                  minSubgroupFrac = 1,
                  permu.size = 100, 
                  raw.pvalue = 1E-4,   
                  Pe = 1E-4, 
                  diffExp = TRUE,
                  filter.probes = FALSE, 
                  filter.percentage = 0.05,
                  filter.portion = 0.3,
                  cores = 1,
                  diff.dir="both",
                  dir.out = '../ELMER'))

#filter for p-value and max distance of 1Mb
pairs <- subset(pairs, pairs$top.dose.vs.NC.diff.pvalue < 0.001)
pairs <- subset(pairs, abs(Distance) < 1e6)

#save plots
CpG <- pairs$Probe
gene_id <- pairs$GeneID

scatter.plot(data = data,
                     byPair = list(probe = c(CpG[1]), gene = gene_id[1]), 
                     category = "GroupLabel", save = TRUE, lm_line = TRUE, save = TRUE) 

schematic.plot(pair = pairs, 
        data = data_Acet,  group.col = "Top_Dose",  byGene = pairs$GeneID[1], save = TRUE  )


#END


