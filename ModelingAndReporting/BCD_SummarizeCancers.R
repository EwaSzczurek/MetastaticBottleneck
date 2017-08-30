source("BCD_FittingLib.R")
source("BCD_DataAccess.R")

calculateRMSE<- function(r.valfits){
	rmses = NULL
	for( cancer in cancers){
		y1 = r.valfits[r.valfits$cancer==cancer,"pred"]
		y = r.valfits[r.valfits$cancer==cancer,"y"]

		is.na.y1=is.na(y1)
		is.na.y = is.na(y)
		is.nas = is.na.y1 | is.na.y
		y1 = y1[!is.nas]
		y = y[!is.nas]
		
		n = length(y)
		rmse = sqrt( sum( (y1 - y )^2 )/n )
		rmses = c(rmses, rmse)
	}
	rmses
}


toinclude=c("medSurvExpV","obsMetProbExpV","PmetExpVCond","PmetExpVConda","PmetExpVCondc","medSurvExpVObsMet")
distribution = "lnorm"
currf = getwd()
setwd(fitFolder)
res.files = dir(pattern = "^[R]")
setwd(currf)
res.files = paste(fitFolder,res.files, sep="")

# These commented out calles generate old plots for own use of the authors, disclaimer: not double checked for errors.
# DrawCancersFuncnam("PmetExpVCond"	,"BCD")
# DrawCancersFuncnam("medSurvExpVObsMet"	,"BCD")
# DrawCancersFuncnam("Pcure"	,"BCD")
# DrawCancersFuncnam("Prem"	,"BCD")
# DrawCancersFuncnam("DBottleneck"	,"BCD")
# DrawCancersFuncnam("ENMets"	,"BCD")
# DrawCancersFuncnam("MedNMets"	,"BCD")
# DrawCancersFuncnam("DistrNMets.curing","BCD")
# DrawCancersDelay("medSurvExpVDiff","BCD",c(100, 50, 25))
# DrawCancersDelay("PmetExpVCondDiff","BCD",c(100, 50, 25))

#### This produces summary plots to the result folder:
#### Achtung! DrawSummaries() takes ages and generates old plots for own use of the authors, disclaimer: not double checked for errors.
#### DrawSummaries()

#### This gathers all parameters per cancer into one df
SIM.RES <- getSimRes()
SIM.RES <- cbind(SIM.RES, medianB = exp(SIM.RES[,"c1"]))
SIM.RES <- cbind(SIM.RES, stdB = sqrt( ( exp( SIM.RES[,"c2"]^2 ) - 1)*exp(2*SIM.RES[,"c1"] + (SIM.RES[,"c2"])^2)))
SIM.RES <- cbind(SIM.RES, avgRemovableAge = 1/SIM.RES[,"f"])

#sapply(colnames(SIM.RES), drawPars, SIM.RES, col = "lightblue",sor=TRUE)
#drawPars("rmse", SIM.RES,  "All datasets", sor = TRUE)

parameter = c(SIM.RES[,"T0"],SIM.RES[,"T1"], SIM.RES[,"h"], SIM.RES[,"avgRemovableAge"] , SIM.RES[,"medianB"],( SIM.RES[,"stdB"]), SIM.RES[,"c1"], SIM.RES[,"c2"])
parameters = data.frame(Time= parameter,  Parameter=c( rep("t0", nrow(SIM.RES)), rep("t1", nrow(SIM.RES)), rep("h", nrow(SIM.RES)), rep("a", nrow(SIM.RES)), rep("median", nrow(SIM.RES)), rep("std", nrow(SIM.RES)), rep("location", nrow(SIM.RES)), rep("scale", nrow(SIM.RES)) ), Cancer = rep(rownames(SIM.RES), 8))
parameters$Parameter = factor( parameters$Parameter , levels = c("t0","t1","h","a","median","std","scale","location"))
par.h = parameters[ parameters$Parameter == "h",]

preds = GetPreds()
save(preds, file = "preds.RData")
load("preds.RData")




########################################################
# This is to generate Figure 2 with the fitted data
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
#r.valfits.MP= r.valfits.MP[!is.na( r.valfits.MP$y),]

rmses.KM = calculateRMSE(r.valfits.MP)
rmses.U = calculateRMSE(r.valfits.U)
rmses.D = calculateRMSE(r.valfits.D)

rmses = cbind(rmses.KM,rmses.U, rmses.D)
rmses.min = apply(rmses, 1, min)
rmses.wmin  = apply(rmses, 1, which.min)
rmses.Type = rep( "Upper bound", length(rmses.wmin))
rmses.Type[rmses.wmin == 1] = "KM-derived"
rmses.Type[rmses.wmin == 3] = "Lower bound"


#rv.km = data.frame(NRMSE = rmses.KM, Data = rep("Metastasis probablity", length(cancers)), Cancer = cancers)
rv = data.frame(NRMSE = rmses.min, Data = rep("Min metastasis probability", length(cancers)), Cancer = cancers, Type = rmses.Type)

r.val = rbind(r.val1, rv)


source("BCD_ValidateHaeno.R")

#### Finally we can call functions drawing the figures we want in the paper.

####This will generate Figure2.pdf, ExtendedData_Figure1.pdf, and ExtendedData_Figure2.pdf
Figure2(r.fits, r.val, r.valfitsAll, km.data, appr, par.h)

###This generates Figure3.pdf
Figure3(parameters, preds, r.valfitsAll ,par.h)

###This generates Figure4.pdf
Figure4(preds, version="main")

###This generates "ExtendedData_Figure3.pdf"
Figure4(preds, version="supp")

### This will generate the trend analysis pvalues table, and the plot of clinical variables dependence on tumor diameter used in Figure 1.
source("BCD_TrendAnalysis.R")




