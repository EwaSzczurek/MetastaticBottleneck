

drawPars<-function(par, pars, tit ="",sor = FALSE, col="blue"){
	pars = pars[,par]
	if (sor)
		pars = sort(pars)
	titn = gsub(" ","",tit)
	pdf(paste(fitFolder,"Par_",par,"_", titn, ".pdf",sep=""), height=3.3,width=4.2)
	par(ps=9,mar=c(6.5, 3,0.8,0.1), mgp=c(1.7,0.4,0))
	ylim <- c(0, 1.2*max(pars))
	## Plot, and store x-coordinates of bars in xx
	xx <- barplot(pars, xaxt = 'n', xlab = '', width = 0.4, ylim = ylim,
              ylab = par, col=col)
	## Add text at top of bars
	text(x = xx, y = pars, label = round(pars,2), pos = 3, cex = 0.7, col = "black")
	## Add x-axis labels 
	axis(1, at=xx, labels=names(pars), tick=FALSE, las=2, line=-0.1, cex.axis=1)
	title(tit)
	dev.off()
}


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

DrawSummaries<- function(model = "BCD"){
	for (rfile in res.files){
		#rfile = res.files[1]
		load(rfile)
	
		#surv.data, data.gen,nam, res.exp, fitFolder,cancer,
		nam = cancer
		nam = gsub(x=nam,replacement="", " ")
		nam=gsub(x=nam,replacement="_", "/")
		gd =getData(cancer, include= c("medSurvExpV","obsMetProbExpV"))
		surv.data = gd$surv.data
		data.gen=gd$data.gen	
		DrawCancerFit(surv.data, data.gen, res.exp, paste("Summary_",nam,sep=""))
		gd =getData(cancer, include= c("PmetExpVCond","PmetExpVConda","PmetExpVCondc","medSurvExpVObsMet"))
		surv.data = gd$surv.data
		data.gen=gd$data.gen	
		DrawCancerVal(surv.data, data.gen, res.exp, paste("Summary_",nam,sep=""))
		
		if (model == "BCD"){
		gd =getData(cancer, include=c("medSurvExpV","obsMetProbExpV","PmetExpVCond","PmetExpVConda","PmetExpVCondc","medSurvExpVObsMet"))
		surv.data = gd$surv.data
		data.gen=gd$data.gen
		DrawCancerPred(surv.data, data.gen, res.exp, paste("Summary_",nam,sep=""))
		}
	}
}

GetCancerPred <- function(surv.data, data.gen, res.exp, cancer){
	
	preds = NULL
	# columns of preds : x, y, pred, additional, Cancer
	
	dg = data.gen[[1]]	
	pars = dg$basic.pars
	pars = getPars(res.exp, pars)

	x = seq(1/(2^3), 10,by=1/(2^3))
	bs = c( 10, 15, 17.5, 20, 22.5, 25)
	
#### BOTTLENECK DISTRIBUTION
	# xl ="Bottleneck density"
	# x = seq(0.1, 40,by=0.1)
	# if (distribution == "lnorm"){
		# plot( x =  x, y = dlnorm( x, meanlog = pars$c1, sdlog =pars$c2), cex=cex.p, ylab=xl, xlab = "b", col = "darkgreen", type="l")
	# }
	# title(xl, cex.main=1.1)
	# mtext("A", side = 3, line = 0.7, outer = F, at = -1.6, adj = 0, padj = NA, cex = 0.9, col = NA, font = NA)

#### E[# mets]
	# xl ="E[# mets]"
	# x = seq(0.1, 10,by=0.1)
	# yc = ENMets(x, pars, nocuring=F)
	# y=ENMets(x, pars, nocuring=T)

	# if (distribution == "lnorm"){
			# plot( x =  x, y = y, cex=cex.p, ylab=xl, xlab = "Tumor diameter (cm)", col = "forestgreen",  type="l")
			# lines(x = x, y = yc, col ="magenta")
	# }
	# title(xl)
	# legend("topright",c("at diagnosis","after curing"),fill=c("darkgreen","magenta"))
	# mtext("B", side = 3, line = 0.7, outer = F, at = -1.6, adj = 0, padj = NA, cex = 0.9, col = NA, font = NA) 


#### DISTRIBUTIONOF THE NUMBER OF METS
	# xl ="P(# mets <= k)"	
	# #tit = "Distribution of the # mets"
	# maxk = 1000
	# ds = c(0.5, 1 , 2, 4, 8)
	# mtpl = NULL
	# for (d in ds){
			# # This gives the distr over 0--maxk values. We need these values 1--maxk
			# distnm = DistrNMets( d, maxk, pars, nocuring = F, delay = 0 )
			# distnm = distnm[2:length(distnm)]
			# mtpl = rbind(mtpl, data.frame( x = seq(1, maxk), y = distnm, pred = rep(xl, maxk), additional=rep(d, maxk)) )
	# }

	# preds = rbind(preds, mtpl)
	
	print("Calculating median preds")
	maxk = 250
	xl ="Median(# mets) at diagnosis"	
	distnm = MedNMets( x,  pars, nocuring = T, delay = 0, maxk=maxk)
	mtpl = data.frame( x = x, y = distnm, pred = rep(xl, length(x)), additional=rep(NA, length(x)), additional2=rep(NA, length(x)) ) 
	preds = rbind(preds, mtpl)
	
	xl ="Median(# mets) after chemo"	
	distnm = MedNMets( x,  pars, nocuring = F, delay = 0, maxk=maxk)
	mtpl = data.frame( x = x, y = distnm, pred = rep(xl, length(x)), additional=rep(NA, length(x)), additional2=rep(NA, length(x))) 
	preds = rbind(preds, mtpl)
	
	
#### P CURE FRACTION OF METS REMOVED	
	
	# xl ="Fraction mets removed"
	# x = seq(0.1, 10,by=0.1)
	# plot( x =  x, y = EPcure( x, pars), cex=cex.p, ylab=xl, xlab = "Tumor diameter (cm)", col = "darkgreen", ylim = c(0,1), type="l")

	# title(xl, cex.main=1.1)
	# mtext("D", side = 3, line = 0.7, outer = F, at = -1.6, adj = 0, padj = NA, cex = 0.9, col = NA, font = NA) 
	
#### P REM PROBABILITY TO REMOVE METASTASES	

	# xl =tit="Probability to remove mets"
	# x = seq(0.1, 10,by=0.1)
	# plot( x =  x, y = Prem( x, pars), cex=cex.p, ylab=xl, xlab = "Met age (years)", col = "darkgreen", ylim = c(0,1), type="l")
	# title(tit, cex.main=1.1)
	# mtext("E", side = 3, line = 0.7, outer = F, at = -1.6, adj = 0, padj = NA, cex = 0.9, col = NA, font = NA) 	
	

#### IMPACT OF TREATMENT DELAY
	# delays = rev( c(4, 8, 16) )
	# delays = c(8)
	# no = 500
	# set.seed(100)
	# bs.draw = c( 10, 15, 18, 20, 22, 25)
	# bs = c( bs.draw, rlnorm(no, meanlog = pars$c1, sdlog = pars$c2) )
	# xl ="% met prob increase due to surg delay"
	# del = NULL
	# for (delay in delays){
		# for (b in bs){
			# delayy = delay * 1/52.1786			
			# y = PmetExpVCondDiffFixedB(x, pars, delayy, fixed.b = b)
			# #y = y*100
			# del = rbind(del, data.frame(x = x, y = y, pred =rep(xl, length(x)), additional = rep(b, length(x))) )
		# }
	# }
	# preds = rbind(preds, del)
	
	#delays = rev( c(4, 8, 16) )
	print("Calculating surg delay preds")
	delays = c(8,16)
	
	#no = 500
	#set.seed(100)
	#bs = c( bs.draw, rlnorm(no, meanlog = pars$c1, sdlog = pars$c2) )
	
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

	print("Calculating immunoth preds")
#### IMPACT OF IMMUNOTHERAPY - INCREASE OF MEDIAN BOTTLENECK	
	bot = NULL				
	xl ="% met prob decrease due to immunotherapy"	
	times.inc = rev( c(1.1, 1.2) )
	for (timesBot in times.inc){
		percBot = paste( (timesBot-1)*100,"%",sep="")
		c1.new = log(timesBot) + pars$c1
		pars.tmp = pars
		pars.tmp$c1 = c1.new
		for (b in bs){
			#y = PmetExpVCondDiffFixedB(x, pars, pars.tmp, delay=0, fixed.b = b)
			y1 = PmetExpVCondND(x, pars, delay = 0, fixed.b = b)
			y2 = PmetExpVCondND(x, pars.tmp, delay = 0, fixed.b = b*timesBot)
			y = (y1-y2) 	
			bot = rbind(bot, 
				data.frame( x = x, y = y, pred =rep(xl, length(x)), additional = rep(b, length(x)), additional2 = rep(percBot, length(x) ) ) ) 
		}
		y1 = PmetExpVCond(x, pars)
		y2 = PmetExpVCond(x, pars.tmp) ### Here, the bottleneck is larger so we should have less mets
		y = (y1-y2) # we want a positive number,we compute the decrease				
		bot = rbind(bot, data.frame( x = x, y = y, pred = rep(xl, length(x)), additional = rep("Marginal", length(x)), additional2 = rep(percBot, length(x))) )
	}
	preds = rbind(preds, bot)
	
	print("Calculating chemo preds")
	
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

MeasureFits<- function(funcnam){
	rmses= c()
	r2s = c()
	cans = c()
	for (rfile in res.files){
		#rfile = res.files[1]
		load(rfile)
	
		#surv.data, data.gen,nam, res.exp, fitFolder,cancer,

		gd =getData(cancer, include= toinclude)
		surv.data = gd$surv.data
		data.gen=gd$data.gen	
		res = MeasureCancerFit(surv.data, data.gen, res.exp, funcnam)
		rmses = c(rmses, res$rmse)
		r2s = c(r2s, res$r2)
		cans = c(cans, cancer)
	}
	df = cbind(rmses, r2s)
	rownames(df)=cans
	colnames(df)=c("RMSE","R2")
	df
}

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
			   #browser()
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
				#le = ra[2]-ra[1]
				#me = mean(y)
				rmse = sqrt( sum( (y1 - y )^2 )/n )
				#nrmse = rmse/me
				#nrmse = rmse/le
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



