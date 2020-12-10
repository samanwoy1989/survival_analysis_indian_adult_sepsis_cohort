# Fixed-effect meta-analysis code
# Input:	eset.list - List of expression sets (one each per study)
# 			group - The phenotype variable name for t-test
#			rmat - Sample Correlation matrix across the studies
# Output:	Meta-analyzed p-values (a vector)

getpmeta = function(eset.list=ss.child, group="Group", rmat=rmat) {
  curregs1 = featureNames(eset.list[[1]])
  for(i in 2:length(eset.list)) {
    curregs = featureNames(eset.list[[i]])
    if(!identical(curregs1, curregs)) {
      print(paste("Different genes in esets 1 and ",i, sep=""))
      return(0)
    }
  }  
  rttlist = lapply(eset.list, function(eset) {
      # row-wise t-test for each gene
      rtt = rowttests(eset, group)
      # Number of samples in each group
      nsamples = table(pData(eset)[,group])
      # convert t-statistic to z-statistic
      degf = sum(nsamples)-2 
      ptres = pt(rtt$statistic, df=degf)
      z = qnorm(ptres)
      # Product of group-size / sum of group-size
      rneff = sqrt(prod(nsamples)/sum(nsamples))
      # return z-statistics (vector of size N) and 
      list("z"=z, "rneff"=rneff)
  })

  # Effective sample size
  rneff = sapply(rttlist, function(x) x$rneff) # k studies vector
  
  # Zstatistic
  z = sapply(rttlist, function(x) x[["z"]]) # N genes-by-k studies matrix

  # Meta-analyze
  numr = z%*%as.matrix(rneff, ncol=1)
  denr = sqrt((matrix(rneff, nrow=1) %*% rmat %*% matrix(rneff, ncol=1)))
  z.meta <- numr/as.numeric(denr)
  pvals <- 2*pnorm(abs(z.meta), lower.tail=F)
  padj <- p.adjust(pvals, method="BH")
  names(padj) = curregs1
  return(padj)
}
