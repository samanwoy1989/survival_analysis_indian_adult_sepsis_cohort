###################################################
# GSEAperm: Gene Set Enrichment Analysis
#				permutation t-test for
#				up- and down-regulation
###################################################
# Create the incidence matrix
pathway.names <- names(genes.by.pathway)
pathway.genes <- unique(unlist(genes.by.pathway))
numPathways <- length(pathway.names)
numGenes <- length(pathway.genes)
Am <- matrix(0, nrow=numPathways, ncol=numGenes)
rownames(Am) <- pathway.names
colnames(Am) <- pathway.genes
for(i in 1:length(genes.by.pathway)) {
  Am[i,genes.by.pathway[[i]]] <- 1
}

# Keep the incidence matrix for genes available in our expression set
Am<- Am[,intersect(featureNames(esetm), colnames(Am))]

# Reduce the incidence matrix by removing all gene sets
# that have 10 genes or fewer
selectedRows = (rowSums(Am)>10)
Am2 = Am[selectedRows, ]

gseap.fn <- 
  paste("Results/gseap_nperm.10000_2017.11.30.rda", sep="")
if(file.exists(file=gseap.fn)) {
  cat("\n Loading GSEA permutation t.test result from file ...\n ")
  load(file=gseap.fn)
  cat(" done!\n")
} else {
  cat("\n Calculating GSEA permutation t.test result from data ...\n")
  gseap <- list()
  set.seed(123)
  NPERM <- 1000
  gseap <- gseattperm(esetm[colnames(Am),], factor(esetm$Group), Am, NPERM)
  save(gseap, file=gseap.fn)
  cat(" done!\n")
}

pGSEAup   <- p.adjust(gseap[, "Upper"], method="bonferroni")
pGSEAdown <- p.adjust(gseap[, "Lower"], method="bonferroni")

################################################
require("GSEABase")
nsF <- esetm
nsF <- nsF[intersect(featureNames(nsF), colnames(Am2)),]
rtt <- rowttests(exprs(nsF), factor(nsF$Group))
rttStat <- rtt$statistic
tA <- as.vector(Am2 %*% rttStat)
tAadj <- tA/sqrt(rowSums(Am2)) # Z score
names(tA) = names(tAadj) = rownames(Am2)

