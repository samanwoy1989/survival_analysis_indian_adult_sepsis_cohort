###################################################
# Figure 2
# Draw box plot of group-mean expression across days
###################################################
xcon <- rowMeans(exprs(eset[upg, which.ctrl]))
xupd1 <- rowMeans(exprs(eset[upg, which.d1]))-xcon
xupd2 <- rowMeans(exprs(eset[upg, which.d2]))-xcon

xcon <- rowMeans(exprs(eset[downg, which.ctrl]))
xdnd1 <- rowMeans(exprs(eset[downg, which.d1]))-xcon
xdnd2 <- rowMeans(exprs(eset[downg, which.d2]))-xcon

par(mar=c(4,5,2,2))
boxplot(list(NA, xupd1, xupd2), col=c("lightblue","red","pink"), outline=F,
	at=c(1,4,6), ylim=c(-3,3), axes=F, 
	ylab="Log-fold change \n (normalized to control value)"#,
	#main="Temporal change in whole blood gene expression \nin adult sepsis"
	)
lines(x=c(1.7,7), y=c(0,0))
boxplot(list(NA, xdnd1, xdnd2), col=c("lightblue", "darkgreen","lightgreen"),
	add=TRUE, outline=F, at=c(1,4,6), axes=F)
axis(1, at=c(4,6), labels=c("SS-Day1", "SS-Day2"), font=2, tick=F)
axis(2); box()
lines(x=c(1,4,6), y=c(0,median(xupd1), median(xupd2)), col="gray50", lty=2,
	lend=2, lwd=3)
lines(x=c(1,4,6), y=c(0,median(xdnd1), median(xdnd2)), col="gray50", lty=2,
	lend=2, lwd=3)
text(x=1,y=0, labels=c("CONTROL"), font=2, col="blue")
legend.up <- paste(length(upg), " genes up-regulated \nin sepsis", sep="")
legend.down <- paste(length(downg), " genes down-regulated \nin sepsis", sep="")
legend(x="topleft", legend=legend.up, text.col="red", cex=1.2)
legend(x="bottomleft", legend=legend.down, text.col="darkgreen", cex=1.2)

pval.upg <- t.test(xupd1, xupd2, paired=TRUE, alternative="greater")$p.value
pval.downg <- t.test(xdnd1, xdnd2, paired=TRUE, alternative="less")$p.value
text.str <- paste("p = ", formatC(c(pval.upg, pval.downg), digits=3), sep="")
text(x=c(5,5), y=c(median(xupd1)+0.5, median(xdnd1)-0.5), 
	labels=text.str, font=4)
