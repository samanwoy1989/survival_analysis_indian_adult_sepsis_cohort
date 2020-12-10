drawPermutHist_ggplot <- function(keggid=keggid, geneids=NULL, titlestr = NULL, nperm=10000, 
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
        # library(ggplot2)
    # Basic histogram
    # ggplot(df, aes(x=perm)) + geom_histogram()
    # Change the width of bins
    # ggplot(df, aes(x=perm)) + 
    #   geom_histogram(binwidth=1)
    # Change colors
    
    theme_update(plot.title = element_text(hjust = 0.5))
    p <- ggplot(df, aes(x=perm)) + 
      geom_histogram(color="black", fill="gray", binwidth = 0.75)
    
    p+ ylab("Frequency")+ xlab("Simulated Pathway Score") + ggtitle(titlestr) + 
      theme(plot.title = element_text(color="red", size=18, face="bold"), 
            axis.title.x = element_text(color="blue", size=12, face="bold"), 
            axis.title.y = element_text(color="blue", size=12, face="bold")) + 
      geom_vline(xintercept=zobs, linetype="solid", color = "red", size=1.2) +
      geom_text(x=zobs+sign(zobs), angle= 90, y=800, label=c(paste0("p.val = ", formatC(permut_pval, digit=1))),
                #labels=c("Observed Score", "Null distribution"), 
                col=c("blue"), size= 6,  fontface="italic")
  }
  return(res)
}

