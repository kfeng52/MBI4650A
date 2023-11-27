#Lecture 9 Live Coding Exercise

#The file needed for today's exercise is the server at /home/shared/Epigenomics ('EdgeR_RNA_all_genes.csv')

#Clear your environment
rm(list = ls())

#RNAseq example with a gene list
library(limma)
library(biomaRt)
library(dplyr)
library(AnnotationDbi)
library(enrichplot)
library(clusterProfiler)
library(DOSE)
library(annotables)
library(TimeSeriesExperiment)
library(org.Hs.eg.db)

#We will start with the "summary statistics" from an experiment, or the results of a differential analysis. 

#First step when working with a microarray or similar is annotating your result to an annotation or manifest file to get your gene names, we will talk about more about this on Monday. For today we don't have to worry about this.

setwd("/Users/Christina/Desktop/Western Assistant Professor/Teaching/MBI4650F")

#We will use a common example for RNA sequencing where your edgeR results have ENSG identifiers for genes:
Seq <- read.csv("EdgeR_RNA_all_genes.csv", header=T) #edgeR results
names(Seq)[names(Seq) == "X"] <- "ensembl_gene_id" #rename first column

#Use Bonferroni correction for p-value
Ntests = nrow(Seq)
pval = 0.05/Ntests
pval
#3.441867e-06

#Filtering and joining with entrez IDs
Seq_sig <- merge(Seq, grch38, by.x="ensembl_gene_id", by.y="ensgene", sort=F)
Seq_sig <- filter(Seq_sig, PValue < pval)
Seq_sig <- Seq_sig[!is.na(Seq_sig$entrez),]
Seq_sig <- Seq_sig[order(Seq_sig$PValue),]
sig_genes <- Seq_sig$entrez

Seq_all <- merge(Seq, grch38, by.x="ensembl_gene_id", by.y="ensgene", sort=F)

length(sig_genes)
dim(Seq_all)

#The Importance of a 'universe'

#The background list or universe in a GO or KEGG analysis (or otherwise), is the list of all genes that you tested in your study. In the example of RNA sequencing, the universe would comprise all genes that were tested in your final model (in this case ~14000 genes)

#GO analysis

#The Gene Ontology (GO) produces a bird’s-eye view of biological systems by building a tree of terms related to biological functions. This is particularly helpful when dealing with results from genome-wide experiments (e.g. transcriptomics) since classifying genes into groups of related functions can assist in the interpretation of results. Rather than focusing on each gene, one by one, the researcher gets access to metabolic pathways, functions related to development, etc.

#The GO resource is divided into 3 main subdomains:
  
#Biological Process (BP): a series of molecular events with a defined beginning and end relevant for the function of an organism, a cell, etc.
#Cellular Component (CC): the part of a cell.
#Molecular Function (MF): the enzymatic activities of a gene product.

#GO Enrichment
go_rna <- goana(sig_genes, universe = Seq_all$entrez, species="Hs")
go_rna <- go_rna %>% arrange(P.DE)
go_rna$GO_label <- rownames(go_rna)

#Can pull top results with topGO
GO_enrich <- topGO(go_rna)
GO_enrich
#How do we interpret the results?

#Visualizing GO results
#Simple method from TimeSeriesExperiment package
pdf("GO_enrich.pdf")
plotEnrichment(GO_enrich, n_max = 15)
dev.off()

#Barplot
d=Seq_sig
geneList = d[,2]
names(geneList) = as.character(d[,7])

de <- names(geneList)[abs(geneList) > 2]

ego <- enrichGO(gene=de, 
                universe=as.character(Seq_all$entrez),
                OrgDb=org.Hs.eg.db,
                ont="ALL",
                pAdjustMethod="BH",
                pvalueCutoff=0.01,
                qvalueCutoff=0.05,
                readable=TRUE)

pdf("GO_bar.pdf")
barplot(ego, showCategory=15)
dev.off()

#Dotplot
deg = names(geneList)[abs(geneList) > 1]
do = enrichGO(deg, OrgDb='org.Hs.eg.db', 
              universe= as.character(Seq_all$entrez), 
              ont           = "ALL",
              pAdjustMethod = "BH",pvalueCutoff  = 0.01,
              qvalueCutoff  = 0.05,
              readable      = TRUE)

pdf("GO_dot.pdf")
dotplot(do, showCategory=20)
dev.off()

#END
