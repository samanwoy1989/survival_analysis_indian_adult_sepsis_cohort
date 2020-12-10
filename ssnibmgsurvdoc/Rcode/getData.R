# Loading data from datapackage ssnibmgsurv
# samanwoy "2019-05-29"
################################
library(ssnibmgsurv)
data(eset)

###################################################
### There are three groups: 
### Control (Zcontrol), at the time of admission or (D1)
### 24 hours later (D2)
### Control has been so named to help with the
### seamless operation of rowttests function (used later)
###################################################
which.ctrl <- which(eset$Group=="Ctrl")
which.d1 <- which(eset$Group=="0")
which.d2 <- which(eset$Group=="24")

grpnames <- as.character(eset$Group)
grpnames <- gsub("Ctrl","Zcontrol", grpnames)
grpnames <- gsub("0","D1", grpnames)
grpnames <- gsub("24","D2", grpnames)
eset$Group <- as.factor(grpnames)


###################################################
### Edit a sample name to make anonymous
### Cases with missing outcome information are survivors
### 90 and 70
###################################################
snames <- sampleNames(eset)
snames[snames=="ctrl_sm"] <- "C1"
sampleNames(eset) <- snames
ptids = as.character(eset$PTID)
ptids[ptids=="ctrl_sm"] = "C1"
eset$PTID = factor(ptids)
eset$Outcome[sampleNames(eset)%in%c("70D1","70D2","90D1","90D2")] = "Surv"

# Day 1 sepsis samples of NIBMG; Do not include samples without outcome information
eset.s <- eset[,eset$Group=="D1" & !is.na(eset$Outcome)]
eset.s$Outcome = factor(as.character(eset.s$Outcome))

