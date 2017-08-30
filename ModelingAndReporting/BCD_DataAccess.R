source("A_Handy.R")
source("BCD_Models.R")
fitFolder = "Fit_BCDf_largeInit/"


core.number = 61 #IMPORTANT! Here you fix the number of cores you wish to use for model fitting on your server! The larger no of cores you can use, the faster this computation will finish. Note that with this no of cores the fitting to 14 cancers took several days. So be patient.
dir.create(fitFolder)

getPars<-function(res.exp,pars){
	pars$T0 = as.numeric( res.exp["T0"] )
	pars$h=as.numeric( res.exp["h"] )
	pars$T1=as.numeric( res.exp["T1"] )
	pars$c1=as.numeric( res.exp["c1"] )
	pars$c2=as.numeric( res.exp["c2"] )
	pars$f=as.numeric( res.exp["f"] )
	if (pars$distr == "uniform"){
		pars$down  = pars$c1
		pars$up = pars$c2
	}else{
		pars$up = upper_l_normal 
		pars$down=0
		}
	pars
}







