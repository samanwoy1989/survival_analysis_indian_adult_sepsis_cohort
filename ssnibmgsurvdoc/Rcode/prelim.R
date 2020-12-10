###############################################################
# This is a code to load initial libreries required for the downstream analysis
###############################################################
options(width=40, digits=3)

###################################################
# Load libraries
###################################################
# Load the libraries
library("Biobase")
library("GEOquery")
library("limma")
library("genefilter")
library("org.Hs.eg.db")
library("hgu133plus2.db")
library("illuminaHumanv2.db")
library("GSEABase")
library("Category")
library("SPIA")
library("KEGGREST")
library("graph")
library("rgl")
library("pca3d")
library("annotate")
library("gplots")
library("stringr")
library("pathview")
library("KEGG.db")
library("ggplot2")
###################################################
### Function Definitions
###################################################
source("Rcode/funcDefs/ids2colors.R") # get color code for given ids
source("Rcode/funcDefs/getORApvals.R") # ORA
source("Rcode/funcDefs/box.plot.KEGG.R") # drawing box plot of pathway gexp
#source("Rcode/funcDefs/getPathview.R") # ORA
source("Rcode/funcDefs/drawPermutHist.R") # GSEA
source("Rcode/funcDefs/drawPermutHist_ggplot.R") # GSEA
source("Rcode/funcDefs/Fisher.test.R")
source("Rcode/funcDefs/getpmeta.R")
#source("Rcode/funcDefs/getpGSEA.R")
#source("Rcode/funcDefs/my.nsFilter.R")
#source("Rcode/funcDefs/my.getP2.R")
source("Rcode/funcDefs/getPval.R")

###################################################
# Load meta data
###################################################
# Gene-by-pathway matrix (KEGG)

# genes-by-pathway list
gbyp.fn <- "metadata/gbyp.rda"
if(file.exists(file=gbyp.fn)) {
  #cat("Loading the gbyp file ...")
  load("metadata/gbyp.rda")
  load(file=gbyp.fn)
  cat("\n")
} else {
  # Uncomment and run the code to get the gene list for KEGG patwhays
  
  # created named list, eg:  path:map00010: 
  # "Glycolysis /   Gluconeogenesis" 
  #  pathways.list <- keggList("pathway", "sav")
  
  # make them into KEGG-style human pathway identifiers
  #  pathway.codes <- sub("path:", "", names(pathways.list))
  
  # subsetting by c(TRUE, FALSE) -- which repeats
  # as many times as needed, sorts through some
  # unexpected packaging of geneIDs in the GENE element
  # of each pw[[n]]
  #  genes.by.pathway <- sapply(pathway.codes,
  #                           function(pwid){
  #                               pw <- keggGet(pwid)
  #                               pw[[1]]$GENE[c(TRUE, FALSE)]
  #                               })
  #  save(genes.by.pathway, pathways.list, file=gbyp.fn)
}
