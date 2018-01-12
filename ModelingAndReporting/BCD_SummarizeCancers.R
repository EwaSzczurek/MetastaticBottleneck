#######################
## This code re-generates all plots in the manuscript
######################
source("A_Handy.R")

### overwrite the set of functions (model equations) to such, which will be used for visualisations of both fit and validation

toinclude=c("medSurvExpV","obsMetProbExpV","PmetExpVCond","PmetExpVConda","PmetExpVCondc","medSurvExpVObsMet")
distribution = "lnorm"
currf = getwd()
setwd(fitFolder)
res.files = dir(pattern = "^[R]")
setwd(currf)
res.files = paste(fitFolder,res.files, sep="")

########################################################
# Gather information contained in Figure 2, 3 and 4, 
# with the fits to the data, parameters and model predictions
########################################################

#### This gathers all parameters per cancer into one df
SIM.RES <- getSimRes()
SIM.RES <- cbind(SIM.RES, medianB = exp(SIM.RES[,"c1"]))
SIM.RES <- cbind(SIM.RES, stdB = sqrt( ( exp( SIM.RES[,"c2"]^2 ) - 1)*exp(2*SIM.RES[,"c1"] + (SIM.RES[,"c2"])^2)))
SIM.RES <- cbind(SIM.RES, a=1/SIM.RES[,"f"])

parameter = c(SIM.RES[,"T0"],SIM.RES[,"T1"], SIM.RES[,"h"], SIM.RES[,"a"] , SIM.RES[,"medianB"],( SIM.RES[,"stdB"]), SIM.RES[,"c1"], SIM.RES[,"c2"])
parameters = data.frame(Time= parameter,  Parameter=c( rep("t0", nrow(SIM.RES)), rep("t1", nrow(SIM.RES)), rep("h", nrow(SIM.RES)), rep("a", nrow(SIM.RES)), rep("median", nrow(SIM.RES)), rep("std", nrow(SIM.RES)), rep("location", nrow(SIM.RES)), rep("scale", nrow(SIM.RES)) ), Cancer = rep(rownames(SIM.RES), 8))
parameters$Parameter = factor( parameters$Parameter , levels = c("t0","t1","h","a","median","std","scale","location"))
par.h = parameters[ parameters$Parameter == "h",]


preds = GetPreds()
save(preds, file = "preds.RData")
load("preds.RData")


funcnams = c("medSurvExpV", "obsMetProbExpV")
r.fits = MeasureFitsTog(funcnams)
r.fits$Type = "Fit"
r.fits<- r.fits[r.fits$Data!="0.5-th time to death quantile w. mets",] #not a fit
r.fits<- r.fits[r.fits$Data!="Metastasis probability",] #not a fit
funcnams.val = c("PmetExpVCond","medSurvExpVObsMet")
r.val1 = MeasureFitsTog(funcnams.val)
r.val1$Type="Validation"


funcnams = c("medSurvExpV", "obsMetProbExpV","PmetExpVCond","PmetExpVConda","PmetExpVCondc","medSurvExpVObsMet")
r.valfitsAll = GetFits(cancers, funcnams, c(0.5))

## Here we generate data needed to compute NRSME with and average of "PmetExpVCond","PmetExpVConda","PmetExpVCondc"
r.valfits.MP = GetFits(cancers, c("PmetExpVCond","PmetExpVConda","PmetExpVCondc"), c(0.5))
r.valfits.U = r.valfits.MP[r.valfits.MP$axtit %in% c("Met probability up"),]
r.valfits.D = r.valfits.MP[r.valfits.MP$axtit %in% c("Met probability down"),]
r.valfits.MP = r.valfits.MP[r.valfits.MP$axtit == "Met probability",]

rmses.KM = calculateRMSE(r.valfits.MP)
rmses.U = calculateRMSE(r.valfits.U)
rmses.D = calculateRMSE(r.valfits.D)

rmses = cbind(rmses.KM,rmses.U, rmses.D)
rmses.min = apply(rmses, 1, min)
rmses.wmin  = apply(rmses, 1, which.min)
rmses.Type = rep( "Upper bound", length(rmses.wmin))
rmses.Type[rmses.wmin == 1] = "KM-derived"
rmses.Type[rmses.wmin == 3] = "Lower bound"

rv = data.frame(NRMSE = rmses.min, Data = rep("Min metastasis probability", length(cancers)), Cancer = cancers, Type = rmses.Type)
r.val = rbind(r.val1, rv)


### run separate code evaluating independent validation of the model on pancreatic cancer survival
source("BCD_ValidateHaeno.R")

#### Finally we can call functions drawing the figures we want in the paper.

#### Generate Figure2.pdf, ExtendedData_Figure1.pdf, and ExtendedData_Figure2.pdf
print("Drawing Figure 2, and ExtendedData_Figure1.pdf- fits and predictions, independent validation, as well as ExtendedData_Figure2.pdf with RMSE for fits to quantile survival")
Figure2(r.fits, r.val, r.valfitsAll, km.data, appr, par.h)

### Generate Figure3.pdf
print("Drawing Figure 3- fitted parameters")
Figure3(parameters, preds, r.valfitsAll ,par.h)

### Generate Figure4.pdf
print("Drawing Figure 4- model predictions, main txt")
Figure4(preds, ver="main")

### Generate "ExtendedData_Figure3.pdf"
print("Drawing Figure 4- model predictions, extended data")
Figure4(preds, ver="supp")

### This will generate the trend analysis pvalues table, and the plot of clinical variables dependence on tumor diameter used in Figure 1.
source("BCD_TrendAnalysis.R")
doTrendSigniff(r.valfitsAll)
drawTrends(r.valfitsAll)

### Subfigure for Figure 1 with survival probability for metastasis
source("BCD_MetInitFigure.R")

