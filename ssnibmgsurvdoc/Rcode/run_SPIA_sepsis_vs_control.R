###################################################
# SPIA: Signaling Pathway Impact Analysis
#				permutation t-test for
#				up- and down-regulation
###################################################
spia.file <- "Results/spia.res_2019.05.13.rda"
if(file.exists(spia.file)) {
  cat(" Loading SPIA \n result from file ... ")
  load(file=spia.file)
  cat(" done!\n")
} else {
  spia.res <- list()
  print.noquote(" Running SPIA... \n Process will take Time...")
  egs.sel <- deg.d1
  lfc.sel <- rttm[egs.sel, "dm"]
  names(lfc.sel)  <- egs.sel
  universeGeneIds <- egs.all
  # SPIA
  res<-spia(de=lfc.sel, all=universeGeneIds, organism="hsa", nB=10000)
  spia.res <- res
  spia.date <- date()
  save(spia.res, spia.date, file=spia.file)
}
res <- spia.res
resall.spia <-  data.frame(paste("hsa",res$ID,sep=""), res$NDE, res$pNDE, res$tA, res$pPERT, res$pG,res$Name)
colnames(resall.spia)<-c( "KEGG_id", "NDE", "pNDE", "tA", "pPERT", "pG_score","KEGG_name")
rownames(resall.spia)<-resall.spia$KEGG_id


pPERT = resall.spia$pPERT
names(pPERT) = resall.spia$KEGG_id



pG = resall.spia$pG_score
names(pG) = resall.spia$KEGG_id

###################################################
###
### spia.file <- "Results/spia.res_2017.04.02.rda"
### if(file.exists(spia.file)) {
###  cat(" Loading SPIA \n result from file ... ")
###  load(file=spia.file)
###   cat(" done!\n")
### } else {
###   spia.res <- list()
###   print.noquote(" Running SPIA... \n Process will take ### Time...")
###   egs.sel <- deg.d1
###   lfc.sel <- rttm[egs.sel, "dm"]
###   names(lfc.sel)  <- egs.sel
###   universeGeneIds <- egs.all
###   # SPIA
###   res<-spia(de=lfc.sel, all=universeGeneIds, ### ### organism="hsa", nB=20000)
###   spia.res <- res
###   spia.date <- date()
###   save(spia.res, spia.date, file=spia.file)
### }
### res <- spia.res
### resall.spia <-  data.frame(paste("hsa",res$ID,sep=""), res$NDE, res$pNDE, res$tA, res$pPERT, res$pG, res$Name, res$Status )
### colnames(resall.spia)<-c( "KEGG_id", "NDE", "pNDE", "tA", "pPERT", "pG_score","KEGG_name", "Status")
### rownames(resall.spia)<-resall.spia$KEGG_id









