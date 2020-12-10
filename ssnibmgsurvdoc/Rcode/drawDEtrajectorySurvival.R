# drawDEtrajectorySurvival.R
###########################################
############################################
# > length(deg.d1)
# [1] 1109
# Use the 1109 DE genes to show difference in trajectories
# between survivor and non-survivor
############################################
#graphics.off()
############################################
which.surv <- which(eset$Outcome=="Surv")
which.nonsurv <- which(eset$Outcome=="Nonsurv")
which.d1.surv <- intersect(which.d1, which.surv)
which.d1.nonsurv <- intersect(which.d1, which.nonsurv)
which.d2.surv <- intersect(which.d2, which.surv)
which.d2.nonsurv <- intersect(which.d2, which.nonsurv)

up.con <- exprs(eset[upg, which.ctrl])
up.d1.surv <- exprs(eset[upg, which.d1.surv])
up.d1.nonsurv <- exprs(eset[upg, which.d1.nonsurv])
up.d2.surv <- exprs(eset[upg, which.d2.surv])
up.d2.nonsurv <- exprs(eset[upg, which.d2.nonsurv])

ucon <- rowMeans(up.con)
ud1s <- rowMeans(up.d1.surv)
ud2s <- rowMeans(up.d2.surv)
ud1n <- rowMeans(up.d1.nonsurv)
ud2n <- rowMeans(up.d2.nonsurv)

xsurv <- c(median(ucon), median(ud1s),
           median(ud2s))-median(ucon)
xnonsurv <- c(median(ucon), median(ud1n),
              median(ud2n))-median(ucon)

down.con <- exprs(eset[downg, which.ctrl])
down.d1.surv <- exprs(eset[downg, which.d1.surv])
down.d1.nonsurv <- exprs(eset[downg, which.d1.nonsurv])
down.d2.surv <- exprs(eset[downg, which.d2.surv])
down.d2.nonsurv <- exprs(eset[downg, which.d2.nonsurv])

dcon <- rowMeans(down.con)
dd1s <- rowMeans(down.d1.surv)
dd2s <- rowMeans(down.d2.surv)
dd1n <- rowMeans(down.d1.nonsurv)
dd2n <- rowMeans(down.d2.nonsurv)

ysurv <- c(median(dcon), median(dd1s),
           median(dd2s))-median(dcon)
ynonsurv <- c(median(dcon), median(dd1n),
              median(dd2n))-median(dcon)

myrange = range(c(xsurv, xnonsurv, ysurv, ynonsurv))+c(-0.3,0.3)

which.d1.surv <- intersect(which.d1, which.surv)
which.d1.nonsurv <- intersect(which.d1, which.nonsurv)
which.d2.surv <- intersect(which.d2, which.surv)
which.d2.nonsurv <- intersect(which.d2, which.nonsurv)
titlestr <- paste("Temporal progression of ", length(deg.d1),
                  " \ndifferentially regulated genes", sep="")
par(mar=c(4,5,4,2))
plot(x=1:3, y=xnonsurv, col="red", type="l",
     lty=2, lwd=3, ylim=myrange,
     xlab="",
     ylab="Log-fold change \n (normalized to control value)",
     axes=F, cex.lab=1.1, main=titlestr)
axis(1, at=1:3, labels=c("Control", "SS-Day1", "SS-Day2"), font=2, tick=F)
axis(2, font=2)
box()
lines(x=1:3, y=xsurv, col="pink", lwd=3)
lines(x=1:3, y=ynonsurv, col="darkgreen",
      lty=2, lwd=3)
lines(x=1:3, y=ysurv, col="green", lwd=3)
legend("topleft",
       legend=c("Non-survivor","Survivor"),
       lty=c(2,1), lwd=3, text.font=2)
