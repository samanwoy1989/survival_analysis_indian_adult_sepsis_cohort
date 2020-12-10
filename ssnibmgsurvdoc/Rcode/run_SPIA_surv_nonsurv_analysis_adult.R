#####################################################################
# Septic Shock data													#
# Signaling Pathway Impact Analysis									#
#####################################################################
spia.file <- "Results/spia.surv.non.surv.adult.res_2018.11.05.rda"
if(file.exists(spia.file)) {
  load(file=spia.file)
} else {
  spia.res <- list()
  print.noquote("Running SPIA... Process will take Time...")
  for(i in 1:length(ss.adult)) {
    study.id <-     pb <- x[, "pPERT"]
    ph <- x[, "pNDE"]
    combinemethod = ifelse(sum(combfunc(pb, ph, "fisher") == 
        x$pG) > sum(combfunc(pb, ph, "norminv") == x$pG), "fisher", 
        "norminv")
    okx <- (ph < 1e-06)
    oky <- (pb < 1e-06)
    ph[ph < 1e-06] <- 1e-06
    pb[pb < 1e-06] <- 1e-06
    plot(-log(ph), -log(pb), xlim = c(0, max(c(-log(ph), -log(pb)) + 
        1, na.rm = TRUE)), ylim = c(0, max(c(-log(ph), -log(pb) + 
        1), na.rm = TRUE)), pch = 19, main = "SPIA two-way evidence plot", 
        cex = 1.5, xlab = "-log(P NDE)", ylab = "-log(P PERT)")
    tr <- threshold/dim(na.omit(x))[1]
    abline(v = -log(tr), lwd = 1, col = "red", lty = 2)
    abline(h = -log(tr), lwd = 1, col = "red", lty = 2)
    if (combinemethod == "fisher") {
        points(c(0, -log(getP2(tr, "fisher")^2)), c(-log(getP2(tr, 
            "fisher")^2), 0), col = "red", lwd = 2, cex = 0.7, 
            type = "l")
    }
    else {
        somep1 = exp(seq(from = min(log(ph)), to = max(log(ph)), 
            length = 200))
        somep2 = pnorm(qnorm(tr) * sqrt(2) - qnorm(somep1))
        points(-log(somep1), -log(somep2), col = "red", lwd = 2, 
            cex = 0.7, type = "l")
    }
    oks <- x[, "pGFWER"] <= threshold
    trold = tr
    tr <- max(x[, "pG"][x[, "pGFdr"] <= threshold])
    if (tr <= trold) {
        tr = trold * 1.03
    }
    if (combinemethod == "fisher") {
        points(c(0, -log(getP2(tr, "fisher")^2)), c(-log(getP2(tr, 
            "fisher")^2), 0), col = "blue", lwd = 2, cex = 0.7, 
            type = "l")
    }
names(ss.adult)[i]
    print.noquote(study.id)
    eset <- ss.adult[[study.id]]
    rtt <- rowttests(eset, "Outcome")
    #rtt <- all.con.rtt[[study.id]]
    sel <- which(rtt$p.value<0.05)
    egs.sel <- rownames(rtt)[sel]
    lfc.sel <- rtt[egs.sel, "dm"]
    names(lfc.sel) <- egs.sel
    universeGeneIds <- rownames(rtt)

    # SPIA
    res<-spia(de=lfc.sel, all=universeGeneIds, organism="hsa", nB=2000,
	  plots=F, beta=NULL, combine="fisher", verbose=FALSE)
    spia.res[[i]] <- res
  }
  spia.date <- date()
  save(spia.res, spia.date, file=spia.file)
}

