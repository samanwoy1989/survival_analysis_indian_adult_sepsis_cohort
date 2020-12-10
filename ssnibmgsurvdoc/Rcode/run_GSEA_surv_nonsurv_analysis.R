# This code performs GSEA on each eset without considering common genes
#########################################################

#########################################################################
# Septic Shock data														#
# Gene Set Enrichment Analysis											#
#########################################################################
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

# Reduce the incidence matrix by removing all gene sets that have 10 genes
# or fewer
selectedRows = (rowSums(Am)>10)
Am2 = Am[selectedRows, ]

##########################################################
##############################################
# Permutation-based GSEA for CHILDREN DATASETS
##############################################

gseap.fn <- 
 "Results/gseap_nperm.child.datasets._2018.11.05.rda"
if(file.exists(file=gseap.fn)) {
  cat("Loading gsea_child data from file ...")
  load(file=gseap.fn)
  cat(" done!\n")
} else {
  gseap <- list()
  for(i in 1:length(ss.child)) {
    cat(names(ss.child)[i])
    nsF <- ss.child[[i]]
    set.seed(123)
    NPERM <- 10000
    Am2new <- Am2[, intersect(featureNames(nsF), colnames(Am2))]
    gseap[[i]] <- res<- gseattperm(nsF, factor(nsF$Outcome), Am2new, NPERM)
    names(gseap)[i] <- names(ss.child)[i]
    cat("\n")
  }
  save(gseap, file=gseap.fn)
}
# get all child up-regulation GSEA permutation p-values
gseap.up.child <- sapply(gseap, function(x) x[,2])
# change the column name
colnames(gseap.up.child) <- paste(colnames(gseap.up.child),"GSEA_UP_Child", sep="_")

gseaFp.up.child <- apply(gseap.up.child, 1, Fisher.test)["p.value",]

# get all child down-regulation GSEA permutation p-values
gseap.dn.child <- sapply(gseap, function(x) x[,1])
# change the column name
colnames(gseap.dn.child) <- paste(colnames(gseap.dn.child),"GSEA_DOWN_Child", sep="_")


gseaFp.dn.child <- apply(gseap.dn.child, 1, Fisher.test)["p.value",]

rm(gseap)

###############################################
# create a GSEA dataframe with 10 coloumns
##############################################
#
# study1.dn study2.dn study3.dn study4.dn study1.up study2.up study3.up study4.up  up.f.product  dn.f.study

gsea.res.child <- cbind(gseap.up.child, gseap.dn.child, gseaFp.up.child,  gseaFp.dn.child)


# save the result in a .rda file
################################
if(file.exists(file="Results/gsea.res.child.rda")) {
  cat("File exists ...\n")
  cat("Done!...\n")
  } else {save(gsea.res.child, file="Results/gsea.res.child.rda")
}

rm(gseaFp.dn.child, gseaFp.up.child, gseap.fn, gseap.up.child, gseap.dn.child)
########################################################################################
# Permutation-based GSEA for ADULT DATASETS
############################################################################################
gseap.fn <- 
  "Results/gseap_nperm.adult.datasets._2018.11.05.rda"
if(file.exists(file=gseap.fn)) {
  cat("Loading gsea_adult data from file ...")
  load(file=gseap.fn)
  cat(" done!\n")
} else {
  gseap <- list()
  for(i in 1:length(ss.adult)) {
    cat(names(ss.adult)[i])
    nsF <- ss.adult[[i]]
    set.seed(123)
    NPERM <- 10000
    Am2new <- Am2[, intersect(featureNames(nsF), colnames(Am2))]
    gseap[[i]] <- gseattperm(nsF, factor(nsF$Outcome), Am2new, NPERM)
    names(gseap)[i] <- names(ss.adult)[i]
    cat("\n")
  }
  save(gseap, file=gseap.fn)
}




# get all adult up-regulation GSEA permutation p-values
gseap.up.adult <- sapply(gseap, function(x) x[,2])
# change the column name
colnames(gseap.up.adult) <- paste(colnames(gseap.up.adult),"GSEA_UP_adult", sep="_")

gseaFp.up.adult <- apply(gseap.up.adult, 1, Fisher.test)["p.value",]

# get all adult down-regulation GSEA permutation p-values
gseap.dn.adult <- sapply(gseap, function(x) x[,1])
# change the column name
colnames(gseap.dn.adult) <- paste(colnames(gseap.dn.adult),"GSEA_DOWN_adult", sep="_")

gseaFp.dn.adult <- apply(gseap.dn.adult, 1, Fisher.test)["p.value",]

rm(gseap)

###############################################
# create a GSEA dataframe with 10 coloumns
##############################################
#
# study1.dn study2.dn study3.dn study4.dn study1.up study2.up study3.up study4.up  up.f.product  dn.f.study

gsea.res.adult <- cbind( gseap.up.adult, gseap.dn.adult, gseaFp.up.adult,  gseaFp.dn.adult)

# save the result in a .rda file
################################
if(file.exists(file="Results/gsea.res.adult.rda")){
  cat("File exists ...\n")
  cat("Done!...\n")
  } else {save(gsea.res.adult, file="Results/gsea.res.adult.rda")
}

#Create a GSEA dataframe
gsea.res.dat <-as.data.frame(cbind(gsea.res.child, gsea.res.adult))


##############################################
# Remove unnecessary variables from workspace
##############################################
rm(Am, Am2, gseap.fn,  gseaFp.up.adult,  gseaFp.dn.adult,  i, numGenes, numPathways, pathway.genes, pathway.names, selectedRows)

