######################################################
# ORA: Over-representation Analysis
######################################################
#
# This Part of the code was used to draw VennDiagram of
# overlapping Genes in Adult and Child Studies
#
#

######################################
# ADULT
######################################
ss.adult = ss.surv.list[study.type=="Adult"]

## Get the common genes of Adult datasets
allegs.ad <- featureNames(ss.adult[[1]])
for (i in 2:length(ss.adult)){
  allegs.ad <- intersect(allegs.ad, featureNames(ss.adult[[i]]))
}

# subset the adult eset list for commongenes
ss.adult <- sapply(ss.adult, function(eset) {
  eset <- eset[allegs.ad,]
})

rmat =  diag(length(ss.adult))	 # rmat is identity matrix
rownames(rmat) = names(ss.adult) # because no overlap of samples among the data sets
colnames(rmat) = names(ss.adult)
pmeta.adult <- getpmeta(ss.adult, group="Outcome", rmat=rmat) # fixed-effect meta-analysis
# select the de genes 
fdr.degns.ad <- names(pmeta.adult[which(pmeta.adult<0.05)])

# lfc for adult
lfc.adult = sapply(ss.adult, function(eset)
  rowttests(eset[allegs.ad,], "Outcome")$dm)
rownames(lfc.adult) = allegs.ad
lfc.adult = lfc.adult[fdr.degns.ad,]

# Number of genes in similar direction in :: Adult
up.adult = names(which(rowSums(lfc.adult>0)>=3))
dn.adult = names(which(rowSums(lfc.adult<0)>=3))
others.adult <- setdiff(rownames(lfc.adult), union(up.adult, dn.adult))

de.adult.gns <- c(rep("UP", length(up.adult)), rep("DOWN", length(dn.adult)), rep("OTHERS", length(others.adult)))
names(de.adult.gns) <- c(up.adult, dn.adult, others.adult)

# ORA of UP genes in Adult datasets
up.paths.ad <- getORApvals(names(de.adult.gns[de.adult.gns=="UP"]), allegs.ad)

# ORA of DOWN genes in Adult datasets
dn.paths.ad <- getORApvals(names(de.adult.gns[de.adult.gns=="DOWN"]), allegs.ad)

# ORA of DE genes (UP+DOWN) in Adult datasets
# Added by SKM to get pNDE for modified SPIA method
pNDE.paths.ad <- getORApvals(names(de.adult.gns[de.adult.gns%in%c("UP","DOWN")]),
	featureNames(ss.adult[[1]]))


# clean-up
#rm(i, id, rmat, pmeta.adult, lfc.adult,  up.mat, dn.mat, up.adult, dn.adult, others.adult)


######################################
# CHILD
######################################
ss.child = ss.surv.list[study.type=="Child"]

## Get the common genes of child datasets
allegs.child <- featureNames(ss.child[[1]])
for (i in 2:length(ss.child)){
  allegs.child <- intersect(allegs.child, featureNames(ss.child[[i]]))
}

# subset the child eset list for commongenes
ss.child <- sapply(ss.child, function(eset) {
  eset <- eset[allegs.child,]
})

rmat =  diag(length(ss.child))	 # rmat is identity matrix
rownames(rmat) = names(ss.child) # because no overlap of samples among the data sets
colnames(rmat) = names(ss.child)
pmeta.child <- getpmeta(ss.child, group="Outcome", rmat=rmat)
# select the de genes 
fdr.degns.child <- names(pmeta.child[which(pmeta.child<0.05)])

# lfc for child
lfc.child = sapply(ss.child, function(eset)
  rowttests(eset[allegs.child,], "Outcome")$dm)
rownames(lfc.child) = allegs.child
lfc.child = lfc.child[fdr.degns.child,]

# Number of genes in similar direction in :: Adult
up.child = names(which(rowSums(lfc.child>0)>=3))
dn.child = names(which(rowSums(lfc.child<0)>=3))
others.child <- setdiff(rownames(lfc.child), union(up.child, dn.child))

de.child.gns <- c(rep("UP", length(up.child)), rep("DOWN", length(dn.child)), rep("OTHERS", length(others.child)))
names(de.child.gns) <- c(up.child, dn.child, others.child)

# ORA of UP genes in child datasets
up.paths.child <- getORApvals(names(de.child.gns[de.child.gns=="UP"]), allegs.child)

# ORA of DOWN genes in child datasets
dn.paths.child <- getORApvals(names(de.child.gns[de.child.gns=="DOWN"]), allegs.child)

# ORA of DE genes (UP+DOWN) in child datasets
# Added by SKM to get pNDE for modified SPIA method
pNDE.paths.child <- getORApvals(names(de.child.gns[de.child.gns%in%c("UP","DOWN")]),
	featureNames(ss.child[[1]]))

# Create an ORA dataframe
#####################################
ora.dat <- as.data.frame(cbind(as.character(up.paths.child[,"p"]) , as.character(dn.paths.child[,"p"]), as.character(up.paths.ad[,"p"]), as.character(dn.paths.ad[,"p"])))
colnames(ora.dat) <- c("up.child", "dn.child", "up.adult", "dn.adult")
rownames(ora.dat) <- rownames(up.paths.child)

# save the dataframe into a rda file
#####################################
#save(ora.dat, file="Results/ORA_dat_surv_vs_nsurv.rda")

# remove unnecessary objects from work space
#####################################
rm( "de.child.gns", "dn.paths.ad", "dn.paths.child",  "getORApvals",  "up.paths.ad", "up.paths.child")
