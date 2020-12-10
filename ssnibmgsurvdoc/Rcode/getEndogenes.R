# getEndogenes.R
# this code creates a list of gene groups for endotypes
# corresponding to endotypes of Sepsis reported in
# Davenport 2016 Lancet Resp Med, Scicluna 2017 Lancet Resp Med, Sweeny 2018 Crit Care Med, Wong 2009 BMC Medicine
# Code written on 06 June 2019

fn = "metadata/endogenes.rda"
if(file.exists(file=fn)) {
  cat(paste0("Loading endotype-specific gene lists from file ... "))
  load(fn)
  cat(" done!\n")
} else {
  require("org.Hs.eg.db")
  datapath = "metadata/genes4endotypes/"
  endogenes = list()

  # cia: coagulopathy-inflammopathy-adaptive Sweeny2018
  cia = read.table(file=paste0(datapath, "cia.txt"), header=TRUE, sep="\t")
  
  # mars2: Scicluna2017, four MARS categories using 2 gene ratio
  mars2 = read.table(file=paste0(datapath, "mars2.txt"), header=TRUE, sep="\t")
  
  # mars140: Scicluna2017, four MARS categories using 140 genes in original data
  #          one is removed from analysis because
  #          LOC100131541 is not annotated in NCBI <https://www.ncbi.nlm.nih.gov/gene/100131541>
  #          On checking three entrez ids are represented by aliases not symbols
  #          > mars140[which(!as.character(unlist(mget(as.character(mars140$EntrezID), org.Hs.egSYMBOL, ifnotfound=NA)))%in%as.character(mars140$Symbol)),1:2]
  #             EntrezID  Symbol
  #          12    80127 CCDC176
  #          16    51816   CECR1
  #          50    29997 GLTSCR2
  #         > mget(c("80127","51816","29997"), org.Hs.egSYMBOL)
  #               $`80127`  [1] "BBOF1"
  #               $`51816`  [1] "ADA2"
  #               $`29997`  [1] "NOP53"
  mars140 = read.table(file=paste0(datapath, "mars140.txt"), header=TRUE, sep="\t")
  
  
  # ped: Wong 2009, 3 subclasses of pediatric septic shock, A, B, C
  #      for one gene, authors' list contains alias, not symbol
  #      > ped[which(!as.character(unlist(mget(as.character(ped$EntrezID), org.Hs.egSYMBOL, ifnotfound=NA)))%in%as.character(ped$Symbol)), 1:2]
  #      EntrezID Symbol
  #      26     5579 PRKCB1
  #      get("5579", org.Hs.egSYMBOL)
  #      [1] "PRKCB"
  ped = read.table(file=paste0(datapath, "ped.txt"), header=TRUE, sep="\t")
  
  # srsde: Davenport2016, two SRS categories of sepsis due to CAP, using a large number of DE genes
  #      The original file does not contain Entrez IDs, but only gene symbols (or aliases)
  #      Additionally, there are duplicated genes
  #
  #      Step 1: remove duplicated genes by keeping gene with highest inter-SRS difference
  srsde = read.table(file=paste0(datapath, "srsde.txt"), header=TRUE, sep="\t")
  gsyms = as.character(srsde$Symbol)
  dupgenes = unique(gsyms[duplicated(gsyms)])
  unigenes = setdiff(unique(gsyms), dupgenes)
  unidat = srsde[srsde$Symbol%in%unigenes,]
  dupinds = as.integer(sapply(dupgenes, function(gene) {
    which.gene = which(srsde$Symbol==gene)
    currdat = srsde[which.gene,c("SRS1","SRS2")]
    sel = names(which.max(apply(currdat, 1, function(x) abs(diff(x)))))
    as.integer(sel)
  }))
  dupdat = srsde[dupinds,]
  srsde = rbind(unidat, dupdat)
  rm(unigenes, dupgenes, unidat, dupdat, gsyms, dupinds)
  
  #      Step 2: map gene symbol/alias to Entrez ID
  gsyms = as.character(srsde$Symbol)
  #     first we tried the following code
  #     > egids = as.character(unlist(mget(gsyms, org.Hs.egSYMBOL2EG, ifnotfound=NA)))
  #     due to the presence of a gene symbol "HBD" mapping to two entrez ids
  #     > which(sapply(mget(gsyms, org.Hs.egSYMBOL2EG, ifnotfound=NA), length)>1)
  #      HBD 
  #      866 
  #     > get("HBD", org.Hs.egSYMBOL2EG)
  #     [1] "3045"      "100187828"
  #     This is the only such entry in the list. The code was modified to include the first mapped Entrez id for a gene
  #     > egids = as.character(unlist(sapply(mget(gsyms, org.Hs.egSYMBOL2EG, ifnotfound=NA), "[", 1)))
  egids = as.character(unlist(sapply(mget(gsyms, org.Hs.egSYMBOL2EG, ifnotfound=NA), "[", 1)))
  #     Some are aliases and not symbols, of these some also map to multiple Entrez ids
  #     The first problem is resolved using org.Hs.egALIASES2EG
  #     The second problme is resolved by taking the first entrez id
  gsyms.unannot = gsyms[is.na(egids)]
  egids[is.na(egids)] = as.character(unlist(sapply(mget(gsyms.unannot, org.Hs.egALIAS2EG, ifnotfound=NA), "[", 1)))
  #     After mapping the gene names with org.Hs.egSYMBOL2EG and org.Hs.egALIASES2EG
  #     there are missing entrez ids
  #     > sum(is.na(egids))
  #     [1] 439
  #     These 439 are not considered.
  which.not.missing = which(!is.na(egids))
  srsde = data.frame("EntrezID"=egids, srsde)
  srsde = srsde[which.not.missing,]
  #     There are 2641 Entrez ID selected after mapping
  #     > dim(srsde)
  #     [1] 2641    5
  rm(which.not.missing, egids, gsyms, gsyms.unannot)
  
  # srs7: Davenport2016, two SRS categories of sepsis due to CAP, using 7 genes
  srs7 = read.table(file=paste0(datapath, "srs7.txt"), header=TRUE, sep="\t")
  
  # srsf6: Burnham2017, two categories of sepsis due to FP, using 6 genes
  srsf6 = read.table(file=paste0(datapath, "srsf6.txt"), header=TRUE, sep="\t")
  
  # save to file
  endogenes = list(cia=cia, mars2=mars2, mars140=mars140, ped=ped, srsde=srsde, srs7=srs7, srsf6=srsf6)
  save(endogenes, file=fn)
}
rm(fn)
