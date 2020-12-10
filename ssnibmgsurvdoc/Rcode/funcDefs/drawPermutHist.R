drawPermutHist <- function(keggid=keggid, geneids=NULL, titlestr = NULL, nperm=10000, 
                           eset=eset, fac, drawHisto=TRUE, lwd=3, cex=1.5) {
  if(is.null(geneids)) {
    pathegs = genes.by.pathway[[keggid]]
    titlestr = KEGGPATHID2NAME[[strsplit(keggid, "hsa")[[1]][2]]]
  } else {
    pathegs = geneids
    titlestr = titlestr
  }
  pathegs = intersect(pathegs, featureNames(eset))
  eset.sub = eset[pathegs,]
  rttStat = rowttests(eset.sub, fac, tstatOnly = TRUE)[["statistic"]]  
  zobs = sum(rttStat)/sqrt(length(pathegs))
  
  permVec <- rep(0, nperm)
  set.seed(123)
  j <- 1L
  while (j < (nperm + 1)) {
    p1 <- sample(fac)
    rstatSim <- rowttests(eset.sub, p1, tstatOnly = TRUE)[["statistic"]]
    permVec[j] = sum(rstatSim)/sqrt(length(pathegs))
    j <- j + 1L
  }
  permut_pval = ifelse(zobs<0, sum(permVec<zobs), sum(permVec>zobs))
  permut_pval = permut_pval/nperm
  res <- c("p.val"=permut_pval,"zobs"=zobs)
  if(drawHisto==TRUE) {
    hist(permVec, xlab="Simulated Pathway Score",col="gray",
         xlim=range(c(permVec, zobs)), 
         main= titlestr)
    abline(v=zobs, lwd=lwd, col="red")
    text(x=zobs+sign(zobs), y=500, srt=90,
         pos= 4, labels=c(paste0("p.val = ", formatC(permut_pval, digit=1))),
         #labels=c("Observed Score", "Null distribution"), 
         col=c("red"), cex=cex )
  }
  return(res)
}

