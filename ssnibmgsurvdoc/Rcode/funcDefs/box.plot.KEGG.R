###############################################################
# This is the code to get Box-plots for a DE pathway
# Input: pathid, direction
# Output: A box plot of control, SS-S1, SS-D2
# Sepsis Genomics Project
# Author: Saroj Kant Mohapatra, NIBMG
# Date: Novemember, 2017
###############################################################
box.plot.KEGG<-function(id="hsa04380", direction="up"){
  is.up <- direction=="up"
  
  genes <- intersect(featureNames(eset),genes.by.pathway[[id]])
  if(is.up) {
    genes <- intersect(genes, upg)
  } else {
    genes <- intersect(genes, downg)
  }

  pathway.name <- strsplit(pathways.list[paste("path:",id,sep="")], 
  	split=" - Homo sapiens")[[1]][1] 
  
  xcon <- rowMeans(exprs(eset[genes, which.ctrl]))
  xd1 <- rowMeans(exprs(eset[genes, which.d1]))
  xd2 <- rowMeans(exprs(eset[genes, which.d2]))

  par(pty="s", mar=c(5,4,4,2))
  if(is.up) {
  	boxcols <- c("lightblue","red","pink")
  } else {
  	boxcols <- c("lightblue","darkgreen","lightgreen")
  }
  	
  boxplot(list(xcon, xd1, xd2), col=boxcols,
		outline=F, at=c(1,4,6), 
		ylim=range(xd1)+c(-2,2), axes=F, 
		ylab="Gene expression (Log 2 scale)",
		main="\nTemporal change in \nwhole blood gene expression\n\n ")
  axis(1, at=c(1,4,6), 
		labels=c("Control","SS-Day1","SS-Day2"),las=2); 
  axis(2); box()
  lines(x=c(1,4,6), 
		y=c(median(xcon), median(xd1), median(xd2)), 
		col="gray50", lty=2, lend=2, lwd=3)
  legend(x="bottomright", legend=pathway.name, 
    text.col=ifelse(is.up, "red2", "darkgreen"), cex=0.7)
  pval <- t.test(xd1, xd2, paired=TRUE, 
    alternative=ifelse(is.up,"greater","less"))$p.value
  text.str <- paste("p = ", formatC(c(pval), digits=3), sep="")
  text(x=c(5,5), y=c(median(xd1)+2.95), 
  labels=text.str, font=4, cex=.8)
}
