## run_hclust.R
## Performs hierarchical clustering with common genes
## all eight studies
## Get the common genes
allegs <- featureNames(ss.surv.list[[1]])
for (i in 2:length(ss.surv.list)){
  allegs <- intersect(allegs, featureNames(ss.surv.list[[i]]))
}

## calculate LFC for each gene
lfc8 = sapply(study.ids, function(sid) {
  eset = ss.surv.list[[sid]][allegs,]
  dm = rowttests(eset, "Outcome")$dm
  names(dm) = allegs
  dm
})

#####################################################
# make each study LFc as Standard Normal distribution
#####################################################
std.norm <- function(x){
    mn.x <- mean(x)
    std.x<- sd(x)
    nw.x <- (x - mn.x)/ std.x
    return(nw.x)
}
lfc8.std.nrm <- apply(lfc8, 2, std.norm)

# Perform Hierarchical Clustering
d = dist(t(lfc8.std.nrm), method="euclidean")
clstr.lfc <- hclust(d, method = "complete")
clstr.labels = clstr.lfc$labels
clstr.lfc$labels = paste(study.type[clstr.labels], ".", clstr.labels,sep="")
clstr.lfc$order <- c("2", "3", "1", "4", "6", "7", "8", "5")

plot(clstr.lfc, labels = NULL, hang = 0.2, cex= 0.8,   
	 sub ="", main="", 
	xlab = "", ylab = "Height")

# Remove unnecessary variables
rm(d, i, std.norm, clstr.lfc, clstr.labels, allegs,  lfc8, lfc8.std.nrm)

