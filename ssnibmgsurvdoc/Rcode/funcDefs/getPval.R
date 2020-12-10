# Define a function to calculate paired t-test p-value for
# survivor versus non-survivor for a given pathway
# smanwoy 
##############################################################
getPval <- function(keggid="hsa04612", drawPlot=F, getSigGenes=F, pval=FALSE) {
    #pathegs<-as.character(unlist(get(keggid, org.Hs.egPATH2EG)))
    pathegs <- genes.by.pathway[[keggid]]
    pathegs <- intersect(featureNames(eset.s), pathegs)
    path.eset<- eset.s[pathegs, ]
    #is.day1 <- path.eset$Group=="0"
    #path.eset <- path.eset[, is.day1]
    
    # Survivor versus non-survivors
    xCon <- rowMeans(exprs(eset[pathegs, eset$Group=="Zcontrol"]))
    xSurv <- rowMeans(exprs(path.eset[, path.eset$Outcome=="Surv"]))
    xNonsurv <- rowMeans(exprs(path.eset[, path.eset$Outcome=="Nonsurv"]))
    
    # Get the genes for which difference is more than 10% of survivor
    sigGenes <- names(which(abs(xSurv-xNonsurv)>0.1*xSurv))
    
    sigGenes.sym <- as.character(unlist(mget(sigGenes, org.Hs.egSYMBOL)))
    
    if(getSigGenes) {
        print(sigGenes.sym)
    }
    
    if(pval== "TRUE"){
        pval = t.test(xSurv, xNonsurv, paired=T)$p.value
        ttlstr = paste0("p = ", formatC(pval, digits=2))
    }else{
        ttlstr = " "
    }
    # Draw scatter-plot if necessary
    require("KEGG.db")
    pathname <- KEGGPATHID2NAME[[strsplit(keggid, "hsa")[[1]][2]]]
    if(drawPlot) {
        par(mfrow=c(1,2),  pty="s")
        mycols <- rep("black", length(pathegs))
        #names(mycols) <- pathegs
        #mycols[sigGenes] <- "blue"
        #mycols <- as.character(mycols)
        
        # draw boxplot First
        #####################
        boxplot(list(xCon, xSurv, xNonsurv), outline=F, names=c("Control", "Survivor", "Nonsurvivor"), 
                main=pathname, boxwex=0.4,
                col= c("chartreuse", "cyan2", "brown1"), notch=F, las=2, ylab="Gene expression at the time of diagnosis")
        
        # draw scatterplot then
        #####################
        plot(x=xSurv, y=xNonsurv,  xlab="Mean expression in Survivors",
             ylab="Mean expression in Non-survivors", 
             main= ttlstr,
             pch=16, col=mycols, cex=0.8)
        abline(0,1, col="blue", lwd=2, lty=2)
    }
    
    # calculate p-value
    #pval
}
