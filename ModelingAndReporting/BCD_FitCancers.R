## The main model fitting code.
## To finish calculations in reasonable time, run this code on multiple cores. 
### The number of cores need to be specified as the core.number variable in H_Handy.R


source("A_Handy.R")
source("BCD_FittingLib.R")
source("BCD_DataAccess.R")
require(foreach)
require(doMC)


## Decide on the distribution of the bottleneck severity b
#distribution = "uniform"
#distribution = "gamma"
distribution = "lnorm"

actual.names = switch(distribution,
	lnorm= c("T0","h","T1","mu","sigma","f", "rmse","r2"),
	gamma = c("T0","h","T1","shape","scale","f", "rmse","r2"),
	uniform = c("T0","h","T1", "c1","c2","f","rmse","r2")
	)


### Fit parameters for a given cancer
FitCancer <- function(cancer, draw=FALSE, nam = ""){
	
	gd = getData(cancer, include=toinclude)
	surv.data = gd$surv.data
	data.gen = gd$data.gen
	
	#### generate initial parameters
	#### parameters that are evenly spread on a grid of values, to make sure some regions of parameter space are not avoided by later random sampling:
	ts = seq(from = 5, to = 15, by = 5)
	hs=c(2)
	T1s= c( 2, 5.5, 9, 12.5)	
	if (distribution=="uniform"){
		c1s = c(10, 15, 25, 40) 
		c2s = c1s  	
	}else{
		 ##gamma
		 c1s = c(10, 20, 30) 
		 c2s = c(0.5, 2, 5) 
		 if (distribution =="lnorm"){
			 c1s = log(c1s)
			 c2s = c(0.01, 0.15, 0.3) 
		 }
	  }
	fs = c(1, 10, 15, 20)
	pars1 <- expand.grid(ts, hs, T1s,  c1s, c2s, fs)
	pars1 <- pars1[ pars1[,3]<pars1[,1], ]
	
	#### parameters that are randomly uniformly sampled:  
	set.seed(1)
	parno = 5
	ts = runif(n=(parno-1), min = 5, max = 15)
	hs=runif(n=2, min = 0.5, max = 2.5)
	T1s= runif( n=(parno-1), 2,12.5)	
	if (distribution=="uniform"){
		c1s = c(10, 15, 25, 40) 
		c2s = c1s  	
	}else{
		#gamma
		c1s = runif(n=parno, 10,  30) 
		c2s = c(0.5, 2, 5) 
		if (distribution =="lnorm"){
			c1s = log(c1s)
			c2s = runif(n=parno,0.01, 0.3)
			#c2s = 1/sqrt(6)
			}
	}
	fs = runif(n=parno, 1, 20)
	pars <- expand.grid(ts, hs, T1s,  c1s, c2s, fs)
	pars <- pars[ pars[,3]<pars[,1], ]
	
	if (distribution=="uniform")
		pars <- pars[ pars[,4]<pars[,5], ]
		
	### Final initial parameterization set: Combine initial parameters on a grid with random initial parameters
	pars = rbind(pars, pars1)	
	
	### Run the fit
	ptm <- proc.time()
	res.exp<- FitSurvivalPMetCondObsV(surv.data=surv.data, pars=pars, data.gen=data.gen, paralelPars)
	elapsed=	proc.time() - ptm
	print(paste( "Fit for",cancer,"took",elapsed["elapsed"],"CPU sec"))
	
	### Report the fit
	nam = cancer
	nam = gsub(x=nam,replacement="", " ")
	nam=gsub(x=nam,replacement="_", "/")
	if(draw){
		 DrawCancerFit(surv.data, data.gen, res.exp, nam)		 
	}
	save( surv.data, data.gen,nam, res.exp, fitFolder,cancer, file=paste(fitFolder,"Res_",nam,".RData",sep=""))

	### Return fitted parameter values
	 if (is.null(res.exp)){
	 	res = rep(NA, length(names.out)) 
	 }else{
		res =  res.exp[names.out]
	}
	res
}

#######This is the main call for model fitting. Variable cancers is defined in A_Handy.R
res = MakeCancerFits( cancers )
