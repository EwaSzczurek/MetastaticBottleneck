##############################
### Code implementing the model. Some very basic model functions are in A_Handy.R
##############################

##############################
### Utility functions
##############################
GetDens<-function(b, c1, c2, distr){
	dens <- switch(distr,
		uniform = 1/(c2-c1),
     	gamma = dgamma(b, shape=c1, scale=c2),
     	lnorm = dlnorm(b, meanlog = c1, sdlog = c2))
    dens
}

doIntegration<-function(integrando, pars){
	out <- tryCatch(
    {	
	    suppressWarnings(integrate( integrando, pars, lower = pars$down, upper  =pars$up) )
    },
    error=function(cond) {             
        return(NULL)
    },
	warning=function(cond) {
        return(NULL)
     },
     finally={}
     )
 	
	if ( !(is.null(out) )){      	
		val =(out$value)
	}else{       	
		val = NA
	}
	val
}
			

########################
###### Cancer death rate
########################

#################### Not Distributed bottleneck 
pmetcond.nd <- function( d, b, pars, delay = 0){ 
	Td = tdExp( d, pars$r )
	Td0= (Td - pars$T0)
	Td0 = max(Td0, 0)
	s1=s1SS(b)
	as.numeric(FT1T2ExpC(Td0, Td+delay, Td+delay, pars$f, pars$r, s1) ) 
}	

#################### Distributed bottleneck 
pmetcond.d <- function( d, pars, delay = 0){ 
	pars$Td = tdExp( d, pars$r )
	pars$Td0= (pars$Td - pars$T0)
	pars$Td0 = max(pars$Td0, 0)
	pars$delay = delay
	
	integrando <- function(b, pars){ 
		s1=s1SS(b)
		dens <- GetDens(b, pars$c1, pars$c2, pars$distr)
		as.numeric( dens*FT1T2ExpC(pars$Td0, pars$Td+pars$delay, pars$Td+pars$delay, pars$f, pars$r, s1) ) 
	} 
	val <-doIntegration( integrando, pars)
	val		
}	

PmetExpVCond <- function( D, pars){
	res=as.numeric( unlist( sapply(D, pmetcond.d, pars)) )
	res	
}

PmetExpVCondND <- function( D, pars, delay, fixed.b){
	res=as.numeric( unlist( sapply(D, pmetcond.nd, fixed.b, pars, delay = delay)) )
	res	
}

PmetExpVCondDiff <- function( D, pars, delay){
	resD=as.numeric( unlist( sapply(D, pmetcond.d, pars, delay=delay)) )
	res	=as.numeric( unlist( sapply(D, pmetcond.d, pars, delay=0)) )
	(resD - res)
}

PmetExpVCondDiffFixedB <- function( D, pars, pars1, delay, fixed.b){
	resD = PmetExpVCondND( D, pars, delay=delay, fixed.b= fixed.b)
	res = PmetExpVCondND( D, pars1, delay=0, fixed.b= fixed.b)
	(resD - res)
}

########################
# Quantile time to death
########################

medSurvExpV<- function( D, pars, delay = 0 ){
	
	solve.d <- function( d ){ 
		pars$Td = tdExp( d, pars$r )
		pars$Tdd= pars$Td + delay
		pars$Td0= pars$Td - pars$T0
		pars$Td0= max(pars$Td0, 0)
		
		medxm<- function(x){
			pars$x = x
			integrando <- function(b, pars){				
				s1= s1SS(b)				
				dens <- GetDens(b, pars$c1, pars$c2, pars$distr)
     		  	res = as.numeric(dens* FT1T2ExpC(pars$Td0, pars$x, pars$Tdd, pars$f, pars$r, s1)/FT1T2ExpC(pars$Td0, pars$Tdd, pars$Tdd, pars$f, pars$r, s1))
				res
			} 
			val = doIntegration( integrando, pars )		
			val - pars$quant
		}
	
		out <- tryCatch(
        	{	       
        		suppressWarnings(	uniroot(medxm, lower=(pars$Td0), upper=pars$Tdd) )   
        },
        	error=function(cond) {
            return(NULL)
        	},
        	warning=function(cond) {       
             return(NULL)
        	},
       		 finally={}
        	)
        	xm=NA
        	if ( !(is.null(out) )){      	
			xm= out$root
		}        	
	
		res= xm +pars$T0 - pars$Td+ pars$h
		res
	}
		
	meds=as.numeric( unlist( sapply(D, solve.d )) )
	meds
}

medSurvExpVDiff<- function( D, pars, delay ){	
	meds = medSurvExpV( D, pars, delay=delay )
	meds.wo.delay = medSurvExpV( D, pars, delay = 0 )
	meds.wo.delay - meds
}

#### We know, that the mets must have occurred before t(d)-T1
medSurvExpVObsMet<- function( D, pars ){	
	solve.d <- function( d ){ 
		pars$Td = tdExp( d, pars$r )
		pars$Td0= (pars$Td - pars$T0)
		pars$Td0 = max(pars$Td0, 0)
		pars$Td1 = max( (pars$Td-pars$T1), pars$Td0) # cannot be smaller than Td0
		
		medxm<- function(x){
			pars$x=x
			integrando <- function(b, pars){
				s1= s1SS(b)
				dens <- GetDens(b, pars$c1, pars$c2, pars$distr)
				res=as.numeric(  dens* FT1T2ExpC(pars$Td0, pars$x, pars$Td, pars$f, pars$r, s1)/FT1T2ExpC(pars$Td0, pars$Td1, pars$Td, pars$f, pars$r, s1)   )
				res
			} 
			
			val = doIntegration( integrando, pars)
			
			val - pars$quant
		}
	
		out <- tryCatch(
        	{	       
        	 suppressWarnings(	uniroot(medxm, lower=pars$Td0, upper=pars$Td1) )   
        },
        	error=function(cond) {
            return(NULL)
        	},
        	warning=function(cond) { 
            return(NULL)
        	},
       		 finally={}
        	)
        	xm=NA
        	if ( !(is.null(out) )){      	
			xm= out$root
		}        	
	
		res= xm +pars$T0 -pars$Td+ pars$h
		res
	}
		
	meds=as.numeric( unlist( sapply(D, solve.d )) )
	
	meds
}


############################
# Metastasis incidence rate
obspmetcond.d <- function( d, pars){ 
		pars$Td = tdExp( d, pars$r )
		pars$td0= (pars$Td - pars$T0)
		pars$td0 = max(pars$td0, 0)
		pars$td1 = pars$Td - pars$T1
		pars$td1 = max(pars$td1, pars$td0) ## Cannot be smaller than td0
		integrando <- function(b, pars){ 
			s1 = s1SS(b)
			dens <- GetDens(b, pars$c1, pars$c2, pars$distr)			
			as.numeric( dens * FT1T2Exp(pars$td0, pars$td1, pars$r, s1) )	
		}
		doIntegration( integrando, pars)		
}	


obsMetProbExpV<- function(D, pars){
	res=as.numeric( unlist( sapply(D, obspmetcond.d, pars)) )
	res
}


############################
## Metastasis removal probability (due to treatment)
############################
EPcured <- function( d, pars){ 

	pars$Td = tdExp( d, pars$r )
	pars$Td0= (pars$Td - pars$T0)
	pars$Td0 = max(pars$Td0, 0)
	
	integrando <- function(b, pars){ 
		s1 = s1SS(b)
		dens <- GetDens(b, pars$c1, pars$c2, pars$distr)		
		as.numeric(  dens*(1 - FT1T2ExpC(pars$Td0, pars$Td, pars$Td, pars$f, pars$r, s1)/FT1T2Exp(pars$Td0, pars$Td, pars$r,  s1)) )	
	}
	doIntegration( integrando, pars)	
}	

EPcure<- function(D, pars){
	res=as.numeric( unlist( sapply(D, EPcured, pars )) )
	res
}

Prem<- function(x, pars){
	exp( - pars$f * x )
}

############################
## Median no of mets
############################
MedNMets<- function( D, pars, nocuring = F, delay = 0, maxk = 1000){
	pars$up = 700
	
	solve.d <- function(d){
		pars$Td = tdExp( d, pars$r )
		pars$Tdd=pars$Td + delay
		pars$Td0= (pars$Td - pars$T0)
		pars$Td0 = max(pars$Td0, 0)
		
		solve.k <- function( k ){ 
			pars$k = k
			
			integrando <- function(b, pars){
				s1= s1SS(b)
				dens <- GetDens(b, pars$c1, pars$c2, pars$distr)
				if (nocuring){
					res = as.numeric(  dens * ( ppois( pars$k, s1* NT1T2Exp(pars$Td0, pars$Tdd, pars$r))))  
				}else{
					res = as.numeric(  dens * (ppois(pars$k, s1* NT1T2ExpC(pars$Td0, pars$Tdd, pars$Tdd, pars$f, pars$r))))
				}
				res
			}					
			doIntegration(integrando, pars)

		}
		
		ds=as.numeric( unlist( sapply(seq(0, maxk), solve.k )) )
		ds[ds==-Inf] = NA
		ds[ds == Inf] = NA
		if (all( is.na( ds)   ) || all(ds == Inf) || all(ds == -Inf)){
			med.k = NA
		}else{
			med.k = 0
			if ( length( which( ds <= 0.5) ) > 0 ){
				med.k = max( which(ds <= 0.5)  )
			}
		}
		med.k
	}
	kms=as.numeric( unlist( sapply(D, solve.d )) )
	
	kms
}

######
# This is the function fitted to data
#####
# data.gen - list of lists
# each element: D, function, basic params (list with names a, r, q, quant)
medSurvPmetCondObsV<- function(D,  data.gen,  T0, h, T1, c1, c2, f){
	res = NULL
	for (i in (1:length(data.gen))){	
		dg = data.gen[[i]]
		func = dg$func
		pars = dg$basic.pars
		pars$T0 = T0
		pars$h = h
		pars$T1= T1
		pars$c1 = c1
		pars$c2 = c2
		pars$f= f
		pars$down = 0
		pars$up = upper_l_normal
		if (pars$distr == "uniform"){
        		pars$down = as.numeric( pars$c1)
        		pars$up = as.numeric(pars$c2)
        	}
		res =  c(res, func( dg$D, pars))
	}
	
	res
}
