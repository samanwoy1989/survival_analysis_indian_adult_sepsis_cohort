## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(fig.show = "hold")
knitr::opts_chunk$set(eval = TRUE)
knitr::opts_chunk$set(cache = TRUE)
knitr::opts_chunk$set(message = FALSE)
knitr::opts_chunk$set(warning = FALSE)
options(width=60)


## ---- eval=FALSE---------------------------------------------------------
## install.packages(pkgs="ssnibmgsurv_1.0.tar.gz", repos=NULL)
## # Now the data package ssnibmgsurv is installed on your computer.
## # Check with the following command:
## library("ssnibmgsurv")


## ---- eval=FALSE---------------------------------------------------------
## install.packages(pkgs="ssgeosurv_1.0.tar.gz", repos=NULL)
## # Now the data package ssgeosurv is installed on your computer.
## # Check with the following command:
## library("ssgeosurv")


## ---- echo=FALSE, fig.align='center', fig.cap="Analysis flow with results: The left arm describes the steps for case-vs-control analysis whereas the right arm describes the steps for identification of pathways associated with survival", out.width = '50%'----
knitr::include_graphics('images/analysis_layout.jpg')


## ----prilim_getdata, echo=TRUE, message=FALSE, warning=FALSE-------------
# Clearing the workspace and close any graphics window if open
rm(list=ls()) 
graphics.off()
# Loading the preliminary libraries including data packages
source("Rcode/prelim.R")
source("Rcode/getData.R")
# Using age- and gender-matched controls
matched.12 <- c("C11","C8","C1","C7","C17","C10","C21","C18","C20","C4","C9",
                "C12","42D1","1D1","8D1","50D1","60D1", "90D1","62D1","70D1",
                "19D1","32D1","14D1","61D1")
esetm <- eset[, matched.12]
rttm <- rowttests(esetm, "Group")
lfcdm<-rttm$dm
pdm <- p.adjust(rttm$p.value, method="BH")
names(pdm) <- rownames(rttm)
names(lfcdm) <- rownames(rttm)
egs.all <- featureNames(esetm)
# Removing some variables that are not to be used for further analysis
rm(snames, ptids)


## ----volcano_plot_case_control, echo=TRUE, fig.align='center', fig.pos= "h",  fig.cap=c("Volcano plot showing 24 percent of the genome perturbed in sepsis compared to healthy control (FDR p < 0.05). This establishes large scale change in gene expression in sepsis, and possible multiple pathways being perturbed.\\label{fig:volcano_plot_case_control}"), fig.height=5,  fig.width=5, message=FALSE, warning=FALSE----
# Draw a volcano plot to show that 24% of the genome
# are perturbed in sepsis (FDR p < 0.05)
# sepsis and matched control samples are used here
source("Rcode/makeplot_1_volcano_for_genomic_storm.R")



## ----De_genes_temporal, echo=TRUE, fig.align='center', fig.cap=c("Temporal change of DE genes (FDR p < 0.05, 2 fold-change or more), there is a non-random trend toward the baseline with time (p-values from paired t-tests are provided in the legend). This is consistent with earlier findings from patients with trauma."), fig.height=5,  fig.width=5, message=FALSE, warning=FALSE----
# Detect the highly significant DE genes of sepsis - FDR p < 0.01; 2-fold
# Draw box-plot to show temporal changes in control vs cases
upg <- egs.all[which(pdm<0.01 & lfcdm>1)]
downg <- egs.all[which(pdm<0.01 & lfcdm<(-1))]
deg.d1<-union(upg,downg)

# Line plot showing slow return or non-survivors to baseline gene expression
source("Rcode/makeplot_2_temporal_plot_for_DEgns.R")



## ----DE_survival_trajectory, echo=TRUE, fig.align='center', fig.cap=c("Temporal change of DE genes (FDR p < 0.05, 2 fold-change or more), there is a non-random trend toward the baseline with time (p-values from paired t-tests are provided in the legend). The delayed return to baseline is associated with non-recovery from sepsis."), fig.height=5,  fig.width=5, message=FALSE, warning=FALSE----
# draw trajectory of DE genes survivor versus non-survivor
source("Rcode/drawDEtrajectorySurvival.R")


## ------------------------------------------------------------------------
# Pathway Analysis: ORA, GSEA, SPIA
# ORA - Over representation analysis
pORAup <- getORApvals(upg, egs.all)[,"p"]
pORAup <- p.adjust(pORAup, method="BH")
pORAdown <- getORApvals(downg, egs.all)[,"p"]
pORAdown <- p.adjust(pORAdown, method="BH")

NIBMG.ORA.disease <- cbind(pORAup, pORAdown)
colnames(NIBMG.ORA.disease) <- paste(colnames(NIBMG.ORA.disease),
                                     "NIBMG.disease", sep="_")

# GSEA - Gene Set Enrichment Analysis
# SPIA - Signaling Pathway Impact Analysis
# Warning running this code will take long time (approx 10 minutes each)
source("Rcode/run_GSEA_sepsis_vs_control.R") # GSEA
source("Rcode/run_SPIA_sepsis_vs_control.R") # SPIA

# Combine result from 3 pathway analyses and print the Down and Up pathways#
pathsDown = intersect(intersect(names(which(pGSEAdown < 0.01)), 
                                names(which(pORAdown < 0.01))), names(which(pG< 0.01)))
pathways.list[ paste0("path:", pathsDown)]
# Removing some variables that are not to be used for further analysis
rm(pathsDown)

pathsUp = intersect(intersect(names(which(pGSEAup < 0.01)), 
                              names(which(pORAup < 0.01))), names(which(pG< 0.01)))
pathways.list[ paste0("path:", pathsUp)]
# Removing some variables that are not to be used for further analysis
rm(pathsUp)


## ----fig-margin,  fig.align='center', fig.cap=c("Temporal boxplot of 3 key pathwayes altered in Sepsis by combined pathway analysis."), echo=TRUE----
par(mfrow=c(1,3))
# Box plot of two pathways down-regulated in sepsis
box.plot.KEGG(id="hsa04612", direction="down") # Antigen processing and presentation
box.plot.KEGG(id="hsa04660", direction="down") # T cell receptor signaling
# Box plot up-regulated pathway: Osteoclast Differentiation
box.plot.KEGG(id="hsa04380", direction="up")



## ----loading_data--------------------------------------------------------
# Section B: Survival analysis 
# Getting genes and pathways associated with survival
# uses NIBMG and published data sets
###################################################
rm(list=ls())
graphics.off()

##########################################################
# Preliminries
source("Rcode/prelim.R")



## ----hierarchical_clustering_in_GEO_data, fig.align='center', fig.cap=c("Survival analysis with eight published data sets: human adults and children with sepsis; hierarchical clustering with log-fold change in gene expression led to evidence of developmental age-specific differential perturbation; i.e., separate clusters for adult and child data sets. In view of this difference, further analysis was confined to to Adult data sets when combined with NIBMG data."), fig.height=5, fig.margin=TRUE, fig.width=6----
#############################
# Get the expression set
source("Rcode/getData.R")
library(ssgeosurv)
data(ss.list) # eight data sets = 4 adult + 4 child
data(ss.surv.list) # 8 datasets with survivors and non-survivors
studies = read.table(file="metadata/studies.txt", header=TRUE, sep="\t")
study.ids = as.character(studies$study.id)
study.type = as.character(studies$age)
names(study.type) = study.ids

##########################################
# Day 1 non-survivor vs survivor with FDR p cutoff 0.01
rtt1 <-rowttests(eset.s, factor(eset.s$Outcome))
sel1 <- rownames(rtt1)[which(rtt1$p.value < 0.01)]
sel1.nibmg.lfc <- rtt1[sel1,]

# Hierarchical clustering of log-fold change of 8 studies
source("Rcode/run_hclust.R")
box()


## ------------------------------------------------------------------------
# Pathway analysis in GEO data
# ORA: over-representation Analysis
source("Rcode/run_ORA_surv_nonsurv_analysis.R")
# Perform permutation-based GSEA analysis
source("Rcode/run_GSEA_surv_nonsurv_analysis.R")
#############################
# Two-way evidence plot adult
source("Rcode/run_SPIA_surv_nonsurv_analysis_adult.R")


## ----two_way_evidence_plot, eval=TRUE, echo=TRUE, fig.align='center', fig.cap=c("KEGG pathways associated with survival in GEO data; NF-kappaB signalling pathway, Osteoclast differentiation and NOD-like receptor signalling pathway."), fig.height=5, fig.width=5----
res <- spia.res[[1]]
resall.adult <-  data.frame(rep(names(ss.adult)[1],  nrow(res)), 
                            res$Name, res$ID, res$NDE, 
                            res$pNDE, res$tA, res$pPERT, 
                            res$pG, res$pGFdr, res$Status) 

col.nm <- as.character(sapply(strsplit(colnames(resall.adult), "res."), "[[", 2))
col.nm[1] <- "Study" 

colnames(resall.adult)<- col.nm

for(i in 2:length(ss.adult)) {
  res <- spia.res[[i]]
  resedited <- data.frame(rep(names(ss.adult)[i], 
                          nrow(res)), res$Name, res$ID, res$NDE, res$pNDE, 
                          res$tA, res$pPERT, res$pG, res$pGFdr, res$Status) 
  
  colnames(resedited)<- col.nm
  resall.adult <- rbind(resall.adult, resedited, deparse.level=0)
}

########################################################
# Calculate Fisher's product of p-values of pertabation for all pathways
keggs <- as.character(unique(resall.adult$ID))

# Create empty vector for capturing fisher product of pG
pPERT.Fp.adult <- vector(mode="numeric", length=length(keggs))
names(pPERT.Fp.adult) <- keggs

for(id in keggs) {
  pvec <- resall.adult[resall.adult$ID==id, "pPERT"]
  pPERT.Fp.adult[id] <- Fisher.test(pvec)["p.value"]
  
}
rm(pvec)

# Combine by Fisher product the two p values 
# for perturbation (pb) and hypergeometric test (ph)
pb <- pPERT.Fp.adult
pb <- p.adjust(pb, "fdr")
ph <- as.numeric(pNDE.paths.ad[paste("hsa",names(pb),sep=""),"p"])
names(ph) <- names(pb)

# Use a floor value for p
ph[ph < 1e-07] <- 1e-07
pb[pb < 1e-07] <- 1e-07

pGmeta.adult <- combfunc(pb,ph, "fisher")

# Capture the fisher product of pPERT pG Meta p values into a dataframe
fisher.prod.spia.adult <- as.data.frame(cbind(paste("hsa", 
                                                    names(pPERT.Fp.adult), sep=""), 
                                              as.numeric(ph), as.numeric(pPERT.Fp.adult), 
                                              as.numeric(pGmeta.adult)))
colnames(fisher.prod.spia.adult) <- c("paths", "pNDE.adult", 
                                      "pPERT.adult", "pG.meta.adult")

#####################
# Following code is derived from SPIA::plotP
# The pG threshold is the p-value 0.05 corrected for the number of 
# pathways being considered
tr= 0.05 
#tr<-  0.05/length(pb)

# plot neg.log.p_PERT against neg.log.p_NDE
plot(-log(ph), -log(pb), col="gray80",
     xlim = c(0, max(c(-log(ph), -log(pb)) +1, na.rm = TRUE)), 
     ylim = c(0, max(c(-log(ph), -log(pb) +1), na.rm = TRUE)), 
     pch = 19, main = "Two-way evidence plot : Adult Sepsis", cex = 1.5, 
     xlab = "Evidence of Over-representation, -log(p_ORA)", 
     ylab = "Evidence of Perturbation, -log(p_PERT)")

# For selected  pathways for visualisation: NLR, NFkB, Osteoclast
##################################
sel.paths.ad <- c("04621", "04064", "04380")
col.vec <- c("red2", "purple2", "darkblue")
points(-log(ph)[sel.paths.ad ], -log(pb)[sel.paths.ad ], pch = 19, col = col.vec, 
       cex = 1.5)
abline(v = -log(tr), lwd = 1, col = "red", lty = 2)
abline(h = -log(tr), lwd = 1, col = "red", lty = 2)
path.nms <- as.character(sapply(strsplit(
  pathways.list[paste("path:hsa", sel.paths.ad, sep="")],
  " -"), "[[", 1))
# Add a legend to the plot         
legend("topright",  title="Pathway Names", 
       legend= paste(sel.paths.ad , path.nms, sep="  "),  
       text.font=2, text.col = col.vec, horiz=F, 
       cex=0.8, pch=20, col= col.vec)



## ----NFKB_box_scatter, fig.align='center', fig.cap=c("NF-kappa B signalling pathway boxplot and scatterplot")----
rm(list=ls())
##########################################################
# Preliminries
source("Rcode/prelim.R")
# Get the expression set
source("Rcode/getData.R")

# show gene expression trend for NF-kB signaling pathway
getPval(keggid="hsa04064", drawPlot=TRUE, getSigGenes=TRUE)



## ----NFKB_permut_histogram_NFKB_pathway, echo=TRUE, fig.align='center', fig.cap=c("NF-kappa B signalling pathway histogram. Permutation bases Gene set enrichment analysis creats a histogram of simulated pathway scores (in gray bars). The red line shows the observed pathway score in nonsurvivors when compared to survivor."), fig.height=4, fig.width=4----
# Draw the permutation histogram of GSEA for NF-kB signaling pathway
drawPermutHist(keggid="hsa04064", eset=eset.s, fac=factor(as.character(eset.s$Outcome)))



## ---- fig.align='center', fig.cap=c("Relative gene expression of the targets of NFkB (i.e., Antigen processing and presentation genes and Immune receptor Genes). The gray peak in the background representsdistribution of all genes in the genome. There is significant down-regulation of the targets in the non-survivors (blue line).")----
################################################################################
# drawing a histogram for NFkB targets
################################################################################
source("Rcode/plotDensityNFkbTargets.R")



## ----M2_func_def---------------------------------------------------------
# check the expression of M1 vs M2 markers in data from SCB cohort
# read the file containing gene IDs
##############
sel.gns.dat <- read.table("metadata/M1_M2_markers.txt", sep="\t", header=T)
m1.gns <- sel.gns.dat[,1]
m1.gns <- intersect(m1.gns, featureNames(eset))
m2.gns <- sel.gns.dat[,2]
m2.gns <- intersect(m2.gns, featureNames(eset))

# function for plotting macrophage-specific gene expression
plotMgexp = function(type="M1", normalizeByControl=FALSE) {
  if(type=="M1") {
    egs = m1.gns 
  } else {
    egs = m2.gns
  }
  
  # M1 gene expression
  gexp.s = rowMeans(exprs(eset[egs, eset$Outcome=="Surv" & eset$Group=="D1"]))
  gexp.ns = rowMeans(exprs(eset[egs, eset$Outcome=="Nonsurv" & eset$Group=="D1"]))
  
  if(normalizeByControl==TRUE) {
    gexp.c = rowMeans(exprs(eset[egs, which.ctrl]))
    gexp.s  = gexp.s-gexp.c
    gexp.ns = gexp.ns-gexp.c
  }
  plot(x = gexp.s, y=gexp.ns, las=1, cex=0.62, col= "blue", pch=16, 
       ylab = "Mean expression in non-survivors", xlab = "Mean expression in survivors",  
       main=paste0(type, "-specific gene expression"))
  abline(0,1, lty=2, xlim= c(-3, 10), ylim=c(-3, 10))
  gsyms <- as.character(unlist(mget(m1.gns, org.Hs.egSYMBOL)))
  #gsyms[intersect(which(gexp.1[,1] < -1), which(gexp.1[,2] < -1) )] <- ""
  text(x = gexp.s, y=gexp.ns, labels=gsyms, cex= 0.58, pos=2, offset = 0.75, font=2)
  
  pval <- t.test(gexp.s, gexp.ns, paired = T)$p.value
  legend.str <- paste("p = ", formatC(pval, digits=1), sep="")
  legend("topleft", legend.str, bty="n", text.font=4)
}



## ----M1_gexp_scatterplot, fig.align='center', out.width= '70%', fig.cap=c("M1 macrophages (classically activated macrophages) are pro-inflammatory, important in host defence against the pathogens, phagocytosis, secretion of pro-inflammatory cytokines and microbicidal molecules. M2 macrophages (alternatively activated macrophages) participate in regulation of resolution of inflammation and repair of damaged tissues. M2-specific under-expression is observed in non-survivors (p = 0.02)."), fig.height=6, fig.width=8----
plotMgexp("M2")



## ----AgPP_box_scplot, fig.align="center", fig.pos= 'h', out.width= '70%', fig.cap=c("Status of AgPP pathway in NIBMG data with Box/Scatterplot "), fig.height=5, fig.width=8----
###################################################
# Survivor versus non-survivor
###################################################
getPval("hsa04612", drawPlot=TRUE, getSigGenes=TRUE) # AgPP


## ----AgPP_permute_hist, fig.align="center", fig.pos= 'h', out.width= '50%', fig.cap=c("Status of AgPP pathway in NIBMG data, permutation based Gene set enrichment histogram.", out.width= '60%')----
drawPermutHist("hsa04612", eset=eset.s, fac=factor(as.character(eset.s$Outcome)))# AgPP



## ----TCR_box_sc_plot, fig.align="center", fig.pos= 'h', out.width= '70%', fig.cap=c("Box plot and scatter plot showing down-regulation of TCR pathway in non-survivors")----
###################################################
# Survivor versus non-survivor
###################################################
getPval("hsa04660", drawPlot=TRUE, getSigGenes=TRUE) # TCR



## ----TCR_permute_hist, fig.align="center",fig.pos= 'h',  out.width= '50%', fig.cap=c("Down-regulation of TCR pathway in SCB cohort, with histogram from permutation based Gene set enrichment test.")----

drawPermutHist("hsa04660", eset=eset.s, fac=factor(as.character(eset.s$Outcome)))# TCR



## ----3_biological_network_construction, width= 0.15, fig.pos= 'b', fig.align="center", fig.cap=c("Three Biological processes found to be differentially enriched in Nonsurvivors. Refer to the main manuscript for more details.")----
net.dat <- read.table("metadata//network_input.csv", sep="\t", header=T)
net.list = with(net.dat, split(x=Egid, f=Process))

par(mfrow=c(3, 1), mar=c(2, 16, 2, 16))
for(i in 1:3) {
  drawPermutHist(geneids=net.list[[i]], 
                 titlestr=names(net.list)[i], 
                 fac=eset.s$Outcome, 
                 eset=eset.s, cex= 0.80)
  box()
}



## ----NIBMg_sepsis_heterogeneity, fig.align="center", fig.cap=c("Barplot showing the score of 3 key modules in each sepsis patient"), fig.height=4, fig.width=12----
# get scores for the three modules from NIBMG data
fn = "Results/nibmgModuleScores.rda"
if(file.exists(fn)) {
  cat("Reading NIBMG module scores from file ...")
  load(fn)
  cat(" done!\n")
} else {
  sids = sampleNames(eset)[-c(which.ctrl)]
  ids.con = sampleNames(eset)[which.ctrl]
  modScores = sapply(names(net.list), function(id) {
    egs = as.character(net.list[[id]])
    sapply(sids, function(sid) {
      eset.curr = eset[,c(sid, ids.con)]
      eset.curr$Group = factor(as.character(eset.curr$Group))
      rtt = rowttests(eset.curr,"Group")
      rttstats = rtt[egs, "statistic"]
      z = sum(rttstats)/sqrt(length(egs))
      return(z)  
    })
  })
  save(modScores, file=fn)
}

# format the data
sids = rownames(modScores)
sids.df = do.call(rbind,strsplit(sids, split="D"))
colnames(sids.df) = c("Pt","Day")
rm(sids)
dat = data.frame(sids.df, sids=rownames(modScores), modScores)

dat.split = with(dat, split(sids, Pt))
dat.split = sapply(dat.split, as.character)

par(mar=c(3,5,2,1))
dat = sapply(dat.split, function(x)
      as.matrix(dat[unlist(x), 
      c("coagulation","immunosuppression","inflammation")]))
dat1 = dat
for(i in 1:length(dat)) {dat1[[i]] = rbind(dat[[i]], c(0,0,0))}
o = order(sapply(dat, function(x) min(x)))
dat1 = dat1[o]
plotdat = t(do.call(rbind, dat1))
mycols = c("blue", "darkgreen", "red")
b = barplot(plotdat, width=1, ylim=c(-15, 15), beside = T, 
        col=mycols, , ylab="Module score",
        border=mycols, las=2, cex.names=0.9)
legend("top", c("Coagulation", "Immunosuppression", "Inflammation"), horiz=T,  
       col=mycols, 
       border=mycols,
       pch=15, pt.cex= 1.7, bty="n", inset=c(-0.05, 0))
box()


## ----NIBMg_sepsis_modscore, fig.align="center", fig.pos= "b", fig.cap=c("Barplot showing the magnitude of 3 key modules in each sepsis patients. Each bar represents mean of patient-level module score with SEM as the error bar. For all three modules, absolute z-scores have been used. It is clear that there is much greater immunnosuppression compared to inflammtion."), fig.height=4, fig.width=6----
require(ggplot2)
dat <- t(abs(plotdat))
apply(t(abs(dat)), 1, median)
df <- dat
sem <- function(x){
  sd(x)/sqrt(length(x))
}
my_mean <- apply(df, 2, mean)
my_sem <- apply(df, 2, sem)
# new data frame for storing the mean and sem
mean_sem_old <- data.frame(means=my_mean, sems=my_sem, group=colnames(df))
mean_sem <- rbind(mean_sem_old[2,], mean_sem_old[1,], mean_sem_old[3,])
rownames(mean_sem) < c("A", "B", "C")
# larger font
theme_set(theme_gray(base_size = 2))
# factorize the variable for legend  
mean_sem$group <- factor(mean_sem$group, levels = mean_sem$group)
ggplot(mean_sem, aes(x=group, y=means, fill=group)) + theme(legend.position="right") +
   geom_bar(stat='identity', width=0.4) + 
  geom_errorbar(aes(ymin=means-sems, ymax=means+sems), width=.12) + 
    scale_fill_brewer(palette="Dark2")+
    xlab('')  + theme(legend.text=element_text(size=12, face="bold")) + 
    theme(legend.title = element_blank())+
    ylab('Absolute module score (z)') + 
    theme(plot.title = element_text(color="red", size=18, face="bold"), 
          axis.title.x = element_text(color="blue", size=12, face="bold"), 
          axis.title.y= element_text(color="blue", size=12, face="bold"))+
  theme(axis.text.x = element_text(size=0.0012)) + 
  theme(axis.text.y = element_text(size=11, face="bold")) +
  geom_hline(yintercept=0, color="black", size= 0.46) + 
  geom_vline(xintercept=0.61, color="black", size= 0.46)


## ---- echo= FALSE--------------------------------------------------------
sessionInfo()

