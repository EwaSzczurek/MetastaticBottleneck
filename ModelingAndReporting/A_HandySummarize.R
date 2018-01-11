## Collect estimated model parameters
getSimRes <- function(){
	SIM.RES = NULL
	cancersn = NULL
	for (rfile in res.files){
		#rfile = res.files[1]
		load(rfile)
	
		#surv.data, data.gen,nam, res.exp, fitFolder,cancer,
		nam = cancer
		cancersn=c(cancersn, nam)
		SIM.RES = rbind(SIM.RES, res.exp)
	}
	
	rownames(SIM.RES)=cancersn
	SIM.RES
}

### Organize model predictions for a given analyzed cancer
GetCancerPred <- function(surv.data, data.gen, res.exp, cancer){
	
	preds = NULL
	# columns of preds : x, y, pred, additional, Cancer
	
	dg = data.gen[[1]]	
	pars = dg$basic.pars
	pars = getPars(res.exp, pars)

	x = seq(1/(2^3), 10,by=1/(2^3))
	bs = c( 17, 19, 22)
		
	print(cancer)	
	print("Calculating cancer death risk increase due to surgery delay")
	delays = c(8,16)
		
	xl ="% met prob increase due to surg delay"
	del = NULL
	for (delay in delays){
		delayy = delay * 1/52.1786
		for (b in bs){
			y = PmetExpVCondDiffFixedB(x, pars, pars, delayy, fixed.b = b)
			#y = y*100
			del = rbind(del, data.frame(x = x, y = y, pred =rep(xl, length(x)), additional = rep(b, length(x)), additional2 = rep(delay, length(x) )) )
		}
		y = PmetExpVCondDiff(x, pars, delayy)
		del = rbind(del, data.frame(x = x, y = y, pred =rep(xl, length(x)), additional = rep("Marginal", length(x) ), additional2 = rep(delay, length(x) ))  ) 
	}
	preds = rbind(preds, del)

	print("Calculating cancer death risk decrease due to vaccines")
#### IMPACT OF IMMUNOTHERAPY - INCREASE OF MEDIAN BOTTLENECK	 SEVERITY
	bot = NULL				
	xl ="% met prob decrease due to immunotherapy"	
	times.inc = rev( c(1.1, 1.2) )
	for (timesBot in times.inc){
		percBot = paste( (timesBot-1)*100,"%",sep="")
		c1.new = log(timesBot) + pars$c1
		pars.tmp = pars
		pars.tmp$c1 = c1.new
		for (b in bs){
			y1 = PmetExpVCondND(x, pars, delay = 0, fixed.b = b)
			y2 = PmetExpVCondND(x, pars.tmp, delay = 0, fixed.b = b*timesBot)
			y = (y1-y2) 	
			bot = rbind(bot, 
				data.frame( x = x, y = y, pred =rep(xl, length(x)), additional = rep(b, length(x)), additional2 = rep(percBot, length(x) ) ) ) 
		}
		y1 = PmetExpVCond(x, pars)
		y2 = PmetExpVCond(x, pars.tmp) ### Here, the bottleneck is stronger so we should have less mets
		y = (y1-y2) # we want a positive number,we compute the decrease				
		bot = rbind(bot, data.frame( x = x, y = y, pred = rep(xl, length(x)), additional = rep("Marginal", length(x)), additional2 = rep(percBot, length(x))) )
	}
	preds = rbind(preds, bot)
	
	print("Calculating cancer death risk decrease due to  chemotherapy efficacy increase")
	
#### IMPACT OF CHEMOTHERAPY - INCREASE OF AVG TREATABLE MET AGE		DrawCancerChemoth
	chem = NULL
	xl ="% met prob decrease due to chemotherapy"	
	times.inc = rev( c(1.1, 1.2) )
	for (timeschem in times.inc){
		percA = paste( (timeschem-1)*100,"%",sep="")
		f.new = pars$f/timeschem	
		pars.tmp= pars
		pars.tmp$f = f.new
		for (b in bs){
			y1 = PmetExpVCondND(x, pars, delay = 0, fixed.b = b)
			y2 = PmetExpVCondND(x, pars.tmp, delay = 0, fixed.b = b)
			y = (y1-y2) 	
			chem = rbind(chem, 
				data.frame( x = x, y = y, pred =rep(xl, length(x)), additional = rep(b, length(x)), additional2 = rep(percA, length(x) ) ) ) 
		}
		y1 = PmetExpVCond(x, pars)
		y2 =PmetExpVCond(x, pars.tmp) ### Here, the bottleneck is larger so we should have less mets
		y = (y1-y2) # we want a positive number,we compute the decrease
		chem = rbind(chem,  	data.frame( x = x, y = y, pred = rep(xl, length(x)),additional = rep("Marginal", length(x)),  additional2 = rep(percA, length(x))) )
	}
	preds = rbind(preds, chem)
	preds$Cancer = rep( cancer, nrow(preds))
	colnames(preds)=c("x", "y", "pred","additional","additional2", "Cancer")
	preds
}

### Organize model predictions for epidemiological data (SEER) for analyzed cancers
GetPreds<-function(){
	preds = NULL
	for (rfile in res.files){
		load(rfile)
		gd =getData(cancer, include=c("medSurvExpV","obsMetProbExpV","PmetExpVCond","PmetExpVConda","PmetExpVCondc","medSurvExpVObsMet"))
		surv.data = gd$surv.data
		data.gen=gd$data.gen
		cancerPreds<- GetCancerPred(surv.data, data.gen, res.exp, cancer)
		preds = rbind(preds, cancerPreds)
	}
	preds
}

### Return the theoretical curves from the model for given data and parameters
GetFitFunc <- function(surv.data, data.gen, res.exp, funcnams, qs=c(0.5)){
	
	le= length(data.gen)

	ds=0
	got = NULL
	for (i in (1:le)){
		dg = data.gen[[i]]
		func = dg$func
		funcnam=dg$funcnam
		
		pars = dg$basic.pars
		x=dg$D
		xs = seq(0.1,10,by=0.1)
		ds2=ds+length(x)
		y = surv.data[(ds+1):ds2,2]
		ys = rep(NA, length(xs))
		have_val = which(as.character(xs)%in%as.character(x))
		#browser()
		ys[have_val] = y		
		ds=ds2
		xl = "Time to death"
		tit = ""
		pars = getPars(res.exp, pars)
		
		get = FALSE
		
		if (funcnam%in%funcnams){
			get=TRUE
			switch(funcnam,
			medSurvExpV = {
				if (pars$quant%in%qs){
					
					axtit = paste(pars$quant,"-th time to death quantile", sep="")
				}else
					get=FALSE
			},
			obsMetProbExpV = {
					
					axtit = "Met incidence"
			},
			PmetExpVCond={
				
				axtit = "Met probability"				
				},
			PmetExpVCondc={
				
				axtit = "Met probability down"
			},
			PmetExpVConda={
				
				axtit = "Met probability up"
			},
			medSurvExpVObsMet={
				
				axtit = "Med time to death w. mets"
			}
			)
			
		}
		if (get){
			preds = func( xs, pars)
			
			g = data.frame( cbind( diameter = xs, y = ys, pred  = preds))
			g$axtit = rep(axtit, length(xs))  
			got = rbind(got, g)
		}
		
		
	}
	
	if (all(c( "Met probability", "Met probability up","Met probability down") %in% got$axtit)){
		if (! "Avg met probability"%in% got$axtit){
			xs = got$diameter[got$axtit == "Met probability"]
			pred = got$pred[got$axtit == "Met probability"]
			mp = got$y[got$axtit == "Met probability"]
			mpu = got$y[got$axtit == "Met probability up"]
			mpd = got$y[got$axtit == "Met probability down"]
			mpavg = apply( data.frame(mp=mp, mpu = mpu, mpd=mpd), 1, mean, na.rm = T)
			g= data.frame( diameter=xs, y = mpavg, pred = pred)
			g$axtit = rep("Avg met probability",length(xs))
			got = rbind(got, g)
		}
	}
	got
}

### Return the theoretical curves from the model for given cancers
GetFits<- function(cancers, funcnams, quants){
	gots = NULL
	for (rfile in res.files){
		#rfile = res.files[1]
		load(rfile)
		if (cancer %in% cancers){
				
			gd =getData(cancer, include= funcnams )
			surv.data = gd$surv.data
			data.gen=gd$data.gen	
			got = GetFitFunc(surv.data, data.gen, res.exp, funcnams=funcnams, qs= quants)
			got = cbind(got, cancer=rep( cancer, nrow(got)))
			
			gots = rbind(gots, got)
		}
	}
	gots
}


### summarize how well theoretical predictions fit to data
MeasureCancerFits <- function(surv.data, data.gen, res.exp, funcnam.cancer, measure = "NRMSE"){
	
	le= length(data.gen)
	ds=0

	rmses=c()
	for (i in (1:le)){
		dg = data.gen[[i]]
		func = dg$func
		funcnam=dg$funcnam
		
		pars = dg$basic.pars
		x=dg$D ## diam
		ds2=ds+length(x)
		if (funcnam%in%funcnam.cancer){			 
			y = surv.data[(ds+1):ds2,2] ###data points
			if (!is.null(res.exp)){		
				pars = getPars(res.exp,pars)
				y1 = func( dg$D, pars)
				
			}
			is.na.y1=is.na(y1)
			is.na.y = is.na(y)
			is.nas = is.na.y1 | is.na.y
			y1 = y1[!is.nas]
			y = y[!is.nas]	
			if (all(y == -1)){
				rmse=NA
				nrmse = NA
			}else{
				n = length(y1) 
				ra = range(y)
				rmse = sqrt( sum( (y1 - y )^2 )/n )
				r2 =( cor( y1, y ))^2

			}
			switch(funcnam, 
				medSurvExpV={
					nam = paste(pars$quant, "-th time to death quantile", sep="" )
				},
				PmetExpVCond={
					nam =  "Metastasis probability"
				},
				obsMetProbExpV={
					nam="Met incidence rate"},
				medSurvExpVObsMet={
					nam = "0.5-th time to death quantile w. mets"
				}
			)

			nams=names(rmses)
			if (measure == "NRMSE")
				rmses=c(rmses, rmse)
			else
				rmses=c(rmses, r2)
			names(rmses)=c(nams, nam)
		}
		ds=ds2
	}
	#print((rmses))
	rmses
}

### Compute how well the model theoretical curves fit to data
MeasureFitsTog<- function(funcnam, measure = "NRMSE"){
	rmses= NULL
	
	for (rfile in res.files){
		#rfile = res.files[1]
		load(rfile)
	
		#surv.data, data.gen,nam, res.exp, fitFolder,cancer,

		gd =getData(cancer, include= toinclude)
		surv.data = gd$surv.data
		data.gen=gd$data.gen	
		res = MeasureCancerFits(surv.data, data.gen, res.exp, funcnam, measure)
		rmse = data.frame( rmse= res, data = names(res), cancer=cancer)
		rmses = rbind(rmses, rmse)
		
	}
	colnames(rmses)=c("NRMSE", "Data", "Cancer")
	
	rmses
}



