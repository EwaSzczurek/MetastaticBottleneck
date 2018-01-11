library("ggplot2")
library("minpack.lm")

### Clean the fit results from such results on initial parameters which did not converge
cleanFits <- function(fits, cn){
	
	if (!is.null(fits)){
		nulls = unlist( lapply(fits, is.null) )
		
		if (!all(nulls)){
			fits = fits[!nulls]
			fits = matrix( unlist( fits), byrow = TRUE, ncol = length(cn) )
			colnames(fits )= cn
			fits = unique(fits)
			fits <- fits[!is.na(fits[,"rmse"]),,drop=F]
			fits <- fits[order(as.numeric(fits[,"rmse"])), ,drop=F]
		}else
			fits<- NULL
	}

	fits
}

### Run the fit on multiple cores
FitSurvivalPMetCondObsV <-function(surv.data, pars, data.gen, paral = FALSE){
	if (paral){
		sim.no=core.number
		registerDoMC(cores = sim.no)
		np = nrow(pars)
		times = ceiling(np/sim.no)

		fits <- foreach (i = 0:sim.no, .combine = rbind) %dopar% {
  			sel = ( ( min( (i*times)+1, np ) ) : ( min( (i*times+times), np )) )
  			if (length(sel)==1){
	  			if(sel==np)
  					return(NULL)
  			}
  			pars.sel = pars[sel,]
  			fits.sel <- apply( pars.sel, 1, runSurvPMetCondObsV, data.gen=data.gen, dat=surv.data)
			cleanFits(fits.sel, names.run)
  		}
  	}else{
		fits <- apply( pars, 1, runSurvPMetCondObsV, data.gen=data.gen, dat=surv.data)
		fits <- cleanFits(fits, names.run)
	}
	
	res = NULL

	if (!is.null(fits)){
		colnames(fits)=names.run
		fits = fits[ order(as.numeric(fits[,"rmse"])), ,drop=FALSE]
		if (nrow(fits)>0){
	 	 	res <- fits[1, ]
	 		res = as.numeric( res[names.out]) 			
			names(res)= names.out
	 	}
	}
	
	res

}
