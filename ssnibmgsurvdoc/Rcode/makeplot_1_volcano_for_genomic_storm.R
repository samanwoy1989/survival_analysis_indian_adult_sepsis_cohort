###################################################
# Figure 1
# Draw a volcano plot to show the differentially
# expressed genes on day 1
###################################################
upg <- egs.all[which(pdm<0.05 & lfcdm>0)]
downg <- egs.all[which(pdm<0.05 & lfcdm<0)]
deg.d1 <- c(upg, downg)

par(mar=c(4,4,3,3))
neglogp <- (-1)*log(pdm, 10)
percde <- formatC(100*length(deg.d1)/as.numeric(nrow(eset)), digits=2)
#title.str <- paste(percde, "% of the genome (", nrow(eset),	" genes) are\n significantly differentially expressed.", sep="")
plot(x=lfcdm, y=neglogp, pch=20,
	cex=0.5, xlab="Log-fold change",
	#main=title.str,
	col.main="blue",
	ylab="- log10 (FDR p)",
	col="gray", xlim=c(-8,8), font.lab=2)
abline(h=(-1)*log(0.05,10), lty=2)
is.lowp.d1 <- pdm<0.05
colvec <- c(rep("red", length(upg)),
			rep("darkgreen", length(downg)))
sellfcdm <- lfcdm[deg.d1]
selneglogp <- neglogp[deg.d1]
points(x=sellfcdm, y=selneglogp, col=colvec,
	pch=20, cex=0.5)
text(x=-5, y=(-1)*log(0.01,10)-0.5, labels="p = 0.05", font=4)
box()
