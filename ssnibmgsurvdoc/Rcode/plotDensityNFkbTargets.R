################################################################################
# drawing a histogram for NFkB targets
################################################################################
# all genes-NFkB targets (control normalised, survivor normalised)
# density plot of background genes 
# density plot of background genes
# Get the control mean expression
xCon <- rowMeans(exprs(eset[,eset$Group=="Zcontrol"]))
xSurv <- rowMeans(exprs(eset[,which(eset$Group=="D1" & eset$Outcome=="Surv")]))
xNsurv <- rowMeans(exprs(eset[,which(eset$Group=="D1" & eset$Outcome=="Nonsurv")]))
gexp <- data.frame("Con"=xCon, "S"=xSurv, "NS"=xNsurv)

# read the targets list
tg.gns <- read.delim("metadata/NFkBtargetsall.txt")
syms <- as.character(tg.gns$Gene)
syms[syms==""] <- "XXX"
egids <- as.character(unlist(mget(syms, org.Hs.egALIAS2EG, ifnotfound=NA), "[", 1))
egids <- egids[!is.na(egids)]
xs <- gexp[,"S"]-gexp[,"Con"]
xns <- gexp[,"NS"]-gexp[,"Con"]
names(xs)=names(xns)=rownames(gexp)

# check AgPresent genes
apGenes <- as.character(tg.gns$Gene[tg.gns$Type=="AntigenPresent"])
apEgids <- as.character(unlist(mget(apGenes, org.Hs.egSYMBOL2EG, ifnotfound=NA)))

# check Immunoreceptor genes
imGenes <- as.character(tg.gns$Gene[tg.gns$Type=="Immunoreceptor"])
imGenes <- imGenes[imGenes!=""]
imEgids <- as.character(unlist(mget(imGenes, org.Hs.egSYMBOL2EG, ifnotfound=NA)))
imEgids <- imEgids[!is.na(imEgids)]
imEgids <- intersect(imEgids, names(xs))


# Find the highest overlap with pathway APC genes / IM receptors
pathways.list[paste0("path:", 
                     names(which.max(
                       sapply(sapply(genes.by.pathway, 
                                     function(x){intersect(x, apEgids)}), length))))]
pathways.list[paste0("path:", 
                     names(which.max(
                       sapply(sapply(genes.by.pathway, 
                                     function(x){intersect(x, imEgids)}), length))))]

# take unoun of two gene sets so that these can be 
# substracted from genome for background
sel.tg.gns <- union(apEgids, imEgids)

# set the background genes
bck.gns <- xns-xs
bck.gns.n <- bck.gns[which(names(bck.gns)%in%sel.tg.gns=="FALSE")]

##################################
d.plot <- density(bck.gns.n)
apc.d <- density(bck.gns[apEgids])
imr.d <- density(bck.gns[imEgids])
sel.tg.gns.d <- density(bck.gns[sel.tg.gns])
##################################
#hist(bck.gns.n, col="gray10", 
#     main="Relative change in gene expression \ncompared to survivors (control normalised)", xlab="")
plot(d.plot, col="black", lwd=1.5, 
     main="", xlab="")
#lines(apc.d,  col="red", 
#     main="Relative change in gene expression \ncompared to survivors (control normalised)", xlab="")
#lines(imr.d, col="blue", 
#     main="Relative change in gene expression \ncompared to survivors (control normalised)", xlab="")
lines(sel.tg.gns.d, col="blue", lwd=1.5,
      #main="Relative change in gene expression \ncompared to survivors (control normalised)", 
      xlab="")
legend("topleft",
       legend=c("NF-kB targets", "Genome"), col=c("blue", "black", xlab=""), lty=1,
       text.font=2)

# test of significance
pval <- as.numeric(t.test(bck.gns[sel.tg.gns], bck.gns.n)$p.value)
#t.test(bck.gns[imEgids], bck.gns.n)$p.value
legend.str <- paste("p = ", formatC(pval, digits=1), sep="")
legend("topright", legend.str, bty="n", text.font=4)

