source("A_Handy.R")
source("BCD_Models.R")
source("BCD_DataAccess.R")
library(survival)
lower= 2.5
upper = 3.4



RetrieveQuantSurv <- function( quant, ds ,pars ){
	pars$quant = quant	
	medSurvExpV( ds, pars)
}

CollectQuantiles<- function(quants, ds, pars){
	appr = NULL
	for (q in quants){
		qs = RetrieveQuantSurv( q, ds, pars )
		appr = cbind(appr, qs)
	}
	colnames(appr )=quants
	### Now average over the diameters ds 
	appr.d = apply(appr, 2, mean)
	appr.d
}

###### values taken from the plots @ Haeno et al - these are the frequencies in months
haeno.au = c(7, 12, 17,24,30,36,42,48)
haeno.au = haeno.au/12 #recalculate in years
names(haeno.au)=cumsum( c(0.15, 0.27, 0.23, 0.09, 0.075, 0.06, 0.04, 0.04 ) )

haeno.ad = c(7, 12, 17,24,30,36,42,48)
haeno.ad = haeno.ad/12
names(haeno.ad)=cumsum( c(0.07, 0.23, 0.18, 0.08, 0.09, 0.06, 0.05, 0.055 ) )


##### Collect the data for the pancreas patients from Haeno et al.
h = read.table("HaenoData.txt", sep = "\t", header  = T)

h =h[h$Survival >= 0 ,]
h$Survival = h$Survival/12
h = h[h$SizePrimary>=lower,] #where readable on the Figures in Haeno et al.
h = h[h$SizePrimary<=upper,]

h$Status = 1
h.autopsy.b = h[1:31, ]
h.adjuvant.b = h[32:nrow(h), ]
km.autopsy.b <- survfit(Surv(Survival, Status) ~ 1, data = h.autopsy.b, conf.type = "log-log")
km.adjuvant.b <- survfit(Surv(Survival, Status) ~ 1, data = h.adjuvant.b, conf.type = "log-log")
h =h[h$Resected != "No" ,]
h = h[h$ALIVE_DEAD%in%c("DEAD","dead","Dead"),]

h[is.na(h$NeoadjTx),"h$NeoadjTx"]<-"NA"
h =h[h$NeoadjTx != "Yes" ,]

h$DateDx = as.Date(as.character( h$DateDx_Surg ), "%m/%d/%Y")


xm = max( h$Survival)
xm = 5

h.autopsy = h[1:7, ]
h.autopsy$Cohort = "Autopsy"

h.adjuvant = h[8:nrow(h), ]
h.adjuvant$Cohort = "Adjuvant"

km.autopsy <- survfit(Surv(Survival, Status) ~ 1, data = h.autopsy, conf.type = "log-log")
km.adjuvant <- survfit(Surv(Survival, Status) ~ 1, data = h.adjuvant, conf.type = "log-log")

h.km = rbind(h.autopsy, h.adjuvant)
km<- survfit(Surv(Survival, Status) ~ Cohort, data = h.km)

### km$surv
### km$time

##### Our approximation of the survival functions from the quantiles
load(paste(fitFolder,"Res_PancreasAdenoCa.RData",sep=""))

### For a given d, what are the quantiles?
quants = seq(0.05, 0.95, by = 0.05)
pars.basic = list( a = NA, r = rs[cancer], q = q, tumor.model = "exp", distr = "lnorm")
pars = getPars(res.exp, pars.basic)



### For every d, retrieve the quantiles
### appr is a matrix with quantile survival for every d, with ds as rows and quants as colums


#h = h.adjuvant
nam.ad = "adjuvant"
nam.au = "autopsy"
ds.adjuvant = h.adjuvant$SizePrimary
ds.autopsy = h.autopsy$SizePrimary
appr.d.au = CollectQuantiles(quants, ds.autopsy, pars)
appr.d.ad = CollectQuantiles(quants, ds.adjuvant, pars)


### Now make a plot of the survival function and of the 1-quantiles



doPlot<- function(tit, appr.d, km, haeno){

	plot( x = km$time, y = km$surv,  ylab="Surviving fraction", xlab = "Years", col = "black",ylim=c(0,1),xlim=c(0,xm),  main=tit, type = "o", cex = 0.4, pch = 20)	
	lines( x = appr.d, y = (1-as.numeric(names(appr.d))),col = "red")
	lines( x = haeno, y = (1-as.numeric(names(haeno))),col = "orange", type = "o", cex = 0.5)
	lines( x = km.adjuvant.b$time, y  = km.adjuvant.b$surv, col = "gray", type = "o", cex = 0.4)
}

printPlot<- function(nam, appr.d, km, haeno){
	pdf(paste( "../HaenoData/HaenoValidationSurvival_",nam,"_v1.pdf",sep=""), width = 2.5, height = 2.6)
	par(ps=9, mar=c(2.6, 2.7,2.3,0.2), mgp=c(1.6,0.4,0))
	doPlot(nam, appr.d, km, haeno)
	legend("topright",legend=c("Data","MFB","Haeno et al."),fill = c("black","red","orange"), cex = 0.7)
	dev.off()
}

#printPlot(nam.ad, appr.d.ad, km.adjuvant, haeno.ad)
#printPlot(nam.au, appr.d.au, km.autopsy, haeno.au)

km.ad = data.frame(Years = km.adjuvant$time, Survival = km.adjuvant$surv, lci = km.adjuvant$lower, uci = km.adjuvant$upper, Cohort = rep("Adjuvant", length(km.adjuvant$time)))
km.au = data.frame(Years = km.autopsy$time, Survival = km.autopsy$surv, lci = km.autopsy$lower, uci = km.autopsy$upper, Cohort = rep("Autopsy", length(km.autopsy$time)))

km.data = rbind(km.ad,km.au)

appr.ad = data.frame(years = appr.d.ad, survival = 1-as.numeric(names(appr.d.ad)) , Cohort = rep("Adjuvant", length(appr.d.ad)), Model = rep("Prediction", length(appr.d.ad)))
haeno.ad.d = data.frame(years = haeno.ad, survival = 1-as.numeric(names(haeno.ad)), Cohort = rep("Adjuvant", length(haeno.ad)), Model = rep("Haeno et al.", length(haeno.ad)) )

appr.au = data.frame(years = appr.d.au, survival = 1-as.numeric(names(appr.d.au)), Cohort = rep("Autopsy", length(appr.d.au)), Model = rep("Prediction", length(appr.d.au)) )
haeno.au.d = data.frame(years = haeno.au, survival = 1-as.numeric(names(haeno.au)), Cohort = rep("Autopsy", length(haeno.au)), Model = rep("Haeno et al.", length(haeno.au)) )

appr = rbind(appr.ad, haeno.ad.d, appr.au, haeno.au.d)


