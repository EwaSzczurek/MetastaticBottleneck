
### Perform the model fit to fata for given set of initial parameters 
runSurvPMetCondObsV <- function( param.init, data.gen, dat ){
	T0.init = param.init[1]
	h.init = param.init[2]
	T1.init = param.init[3]
	c1.init= param.init[4]
	c2.init=param.init[5]
	f.init = param.init[6]
	dat <- dat[,c("Size","SurvMonths")]	

	out <- tryCatch(
        {	
        	### Call the nonlinear least squares optimization with the Levenberg-Marquardt algorithm
      suppressWarnings( nlsLM(SurvMonths ~ medSurvPmetCondObsV( Size, data.gen, T0, h, T1, c1, c2, f), data= dat, 
         start  = list(  T0 = T0.init, h = h.init, T1 = T1.init, c1= c1.init, c2=c2.init, f = f.init), 
        	 lower = c(0, 0, 0, 0.01, 0.001,  0.92), 
        	 upper = c(Inf, Inf, Inf,  Inf, Inf, Inf) ))  
        	 #the lower bound on f = 1/a (a being the expected age of met becoming irremovable)  assures that after 5 years the probability of removing a met is not larger than 1%.
        },
        error=function(cond) {
            return(NULL)
        },
        warning=function(cond) {
             return(NULL)
        },
        finally={
        }
    )
	### If the optimization converged, we get the set of the parameters found starting from the set of inital parameters. Here we return them together with the rmse and r^2 of the fit.
    if (!is.null(out)){
    		est = coef(out)
    	 	if (sum(is.na(est))>0){
    	 		est = param.init
    	 	}
    		y1 = predict(out)
    		n = length(y1)
    		y = dat[,"SurvMonths"]
    		rmse = sqrt( sum( (y1 - y )^2 )/n )
    		r2 = sqrt( cor( y1, y ))
 		out = c( param.init, est, rmse, r2 )
    		names(out) = names.run
    }
    return(out)
}

