# runEndo.R
# Endotype-specific analysis
# CIA     Sweeny2018    Check for Gene group expression changing in the predicted direction
# MARS2   Scicluna2017  Check gene expression ratio for each MARS category
# MARS140 Scicluna2017  Check gene expression correlation with given MARS category-mean expression; PCA
# PED     Wong2009      PCA
# SRSDE   Davenport2016 PCA
# SRS7    Davenport2016 PCA
# SRSF6   Burnham2017   PCA

# if endogenes object is absent, get it
if(!exists("endogenes")) {
  source("Rcode/getEndogenes.R")
}

# Should the gene expression data be normalized to mean control group expression
normalize2con = TRUE
# get the gene expression data for the patients
which.case = c(which.d1, which.d2)
sids = sampleNames(eset[,which.case])

meanCon = rowMeans(exprs(eset[,which.ctrl]))
gexp = exprs(eset[,sids])
if(normalize2con==TRUE) {
  gexp = gexp-meanCon
}

###################################################################
# CIA   Check with specific genes direction; derived from Sweeny2018
###################################################################
# ciaGeneScore: a cia internal function for getting gene score (based on t-test and log-fold change) 
# input: gene id (eg), patient id (id)
# output: gene socre scalar (+1, -1 or 0)
ciaGeneScore = function(eg, id) {
  which.id = which(sampleNames(eset)==id)
  eset.sel = eset[, c(which.ctrl, which.id)]
  rtt = rowttests(eset.sel[eg, ], "Group")
  pval = rtt[eg, "p.value"]
  gscore = 0
  lfc = rtt[eg, "dm"]
  gscore = sign(lfc)
  return(gscore)
}

# ciaDirection2integer: a cia internal function for converting gene direction to signed integer
# input: cia category (group)
# output: vector of integers (+1 or -1)
ciaDirection2integer = function(group) {
  attach(endogenes)
  sel = which(cia$Category==group) 
  egs = cia$EntrezID[sel]
  result = c(-1,1)[factor(cia$Direction[sel])]
  names(result) = egs
  detach(endogenes)
  return(result)
}

# ciaGroupScore: a cia internal function
# input: cia category id (ciagroup), patient id (id)
ciaGroupScore = function(ciagroup, id) {
  delta = ciaDirection2integer(ciagroup)
  egs = names(delta)
  scores = sapply(egs, ciaGeneScore, id)
  grpscore = sum(scores==delta)*100/length(delta)
  return(grpscore)
}

# generate the patient-specific category scores
nibmg.cia = t(sapply(sort(sids), function(id) {
  sapply(as.character(unique(endogenes$cia$Category)), function(ciagroup) {
    ciaGroupScore(ciagroup, id)
  }) 
}))

# make binary result for each sample
nibmg.cia.binary = t(apply(nibmg.cia, 1, function(x) {
  res = rep(0, length(x))
  sel = which.max(x)
  res[sel] = 1
  res
    
}))
colnames(nibmg.cia.binary) = colnames(nibmg.cia)
nibmg.cia = nibmg.cia.binary


###################################################################
# MARS2   Check ratio
###################################################################
# mars2scores: a mars2 internal function for getting mars2 scores for a given patient id
# input: patient id (id)
# output: a vector of 4 values for 4 MARS categoris
mars2scores = function(id) {
  attach(endogenes)
  mnames = paste0("Mars", 1:4)
  mscores = rep(NA, 4)
  names(mscores) = mnames
  for(i in 1:length(mscores)) {
    numr.eg = as.character(mars2$EntrezID[mars2$Category== mnames[i] & mars2$Math=="Numerator"])
    denr.eg = as.character(mars2$EntrezID[mars2$Category== mnames[i] & mars2$Math=="Denominator"])
    mscores[i] = 2^(gexp[numr.eg,id]-gexp[denr.eg,id])
  }
  detach(endogenes)
  return(mscores)
}

ratio.mars2 = t(sapply(sort(sids), mars2scores))

nibmg.mars2 = t(apply(ratio.mars2, 1, function(x){
  sel = which.max(x)
  y = rep(0, length(x))
  y[sel] = 1
  return(y)
}))
colnames(nibmg.mars2) = colnames(ratio.mars2)
colnames(nibmg.mars2) = paste0(colnames(nibmg.mars2), ".2gr")

#scan()
###################################################################
# MARS140   check for correlation with our data
###################################################################
# get data for endotype-associated genes; intersect with NIBMG expression set
dat.mars140 = endogenes[["mars140"]]
egs.mars140 = dat.mars140$EntrezID
rownames(dat.mars140) = egs.mars140

egs.mars140 = intersect(egs.mars140, featureNames(eset))
dat.mars140 = dat.mars140[egs.mars140,]

gexp.mars140 = gexp[egs.mars140,]
r.mars140 = t(apply(gexp.mars140, 2, function(x) {
  apply(dat.mars140[, paste0("Mars",1:4)], 2, function(y) {
    rx=cor.test(x, y)
    res=ifelse(rx$p.value<0.01, rx$estimate, 0)
    return(res)
  })
}))
require(corrplot)
o = order(rownames(r.mars140))
r.mars140 = r.mars140[o,]

# find the ids without any category result; correlation was insignificant
# for each of the 4 MARS categories
noCategory = names(which(apply(r.mars140, 1, max)==0))

nibmg.mars140 = t(apply(r.mars140, 1, function(x){
  sel = which.max(x)
  y = rep(0, length(x))
  y[sel] = 1
  return(y)
}))
nibmg.mars140[noCategory,] = 0
colnames(nibmg.mars140) = colnames(r.mars140)

cols3 = rainbow(3)
plotdat = as.matrix(data.frame(nibmg.cia, nibmg.mars2, nibmg.mars140))
heatmap.2(plotdat, trace="none", mar=c(8,22),  dendrogram = "none",
          ColSideColors = c(rep(cols3[1],3),rep(cols3[2],4),rep(cols3[3],4)),
          Rowv = F, Colv = F, key=FALSE, keysize=0.7,
          col=c("white","darkblue"), 
          rowsep=0:nrow(plotdat), colsep=0:ncol(plotdat), cexRow = 0.7,
          sepcolor="black", sepwidth=c(.01,.01))
legend(x=0.68, y=1, legend=c("CIA","MARS with 2-gene ratio","MARS with 140 genes"),
       col=cols3, bty="n", pch=15, horiz=F)

# write the result from 3 algorithms to a file
#write.table(data.frame(nibmg.cia, nibmg.mars2, nibmg.mars140), 
#          file="Results/nibmgEndotypePrediction.txt",
#          sep="\t", quote=F, row.names=TRUE, col.names=TRUE)

nibmg.cia.vector = as.integer(apply(nibmg.cia, 1, function(x) which(x==1)))
names(nibmg.cia.vector) = rownames(nibmg.cia)

nibmg.mars2.vector = as.integer(apply(nibmg.mars2, 1, function(x) which(x==1)))
names(nibmg.mars2.vector) = rownames(nibmg.mars2)

nibmg.mars140.vector = as.integer(apply(nibmg.mars140, 1, function(x) which(x==1)))
names(nibmg.mars140.vector) = rownames(nibmg.mars140)

#write.table(data.frame(nibmg.cia.vector, nibmg.mars2.vector, nibmg.mars140.vector), 
#          file="Results/nibmgEndotypePredictionStatus.csv",
#          sep="\t", quote=F, row.names=TRUE, col.names=TRUE)


#scan()
###################################################################
# PCA for MARS140, PED, SRSDE, SRS7, SRSF6
###################################################################
# endotypes for PCA only
for(endo in c("mars140","ped","srsde","srs7", "srsf6")) {
  cat(paste0("Now doing for the endotype ", endo, "\n"))
  
  # get data for endotype-associated genes; intersect with NIBMG expression set
  endoegs = intersect(as.character(endogenes[[endo]]$EntrezID), featureNames(eset))
  gexp.endo = gexp[endoegs,]
  
  # PCA
  pca = prcomp(t(gexp.endo))
  for(i in 1:3) {
    ids = pData(eset)[which.case,c("Outcome","Group","PTID")[i]]
    mycols=ids2colors(ids)
    legendcols = unique(mycols)
    names(legendcols) = names(mycols)[match(legendcols, mycols)]
    pca3d(pca, col=mycols)
    legend3d("topright", names(legendcols), col=legendcols, pch=19, pt.cex=2, bty="n", horiz=FALSE)
    #fn = paste0("Results/endogenes/", endo, "_", c("survival","day","subject")[i],".png")
    #snapshotPCA3d(file=fn)
    #scan()
    }
}

rm("ciaDirection2integer", "ciaGeneScore", "ciaGroupScore", "dat.mars140", "egs.mars140", "endo", "endoegs", "endogenes", "gexp", "gexp.endo", "gexp.mars140", "i", "ids", "legendcols", 
"mars2scores", "meanCon", "mycols", "noCategory", "normalize2con", "o", "pca", "r.mars140", "sids", "which.case")

