load("Cancers.RData")


MeasureCancerFit <- function(surv.data, data.gen, res.exp, funcnam.cancer){
	
	le= length(data.gen)
	ds=0
	ys = c()
	y1s=c()
	for (i in (1:le)){
		dg = data.gen[[i]]
		func = dg$func
		funcnam=dg$funcnam
		pars = dg$basic.pars
		x=dg$D ## diam
		ds2=ds+length(x)
		if (funcnam==funcnam.cancer){
			 
			y = surv.data[(ds+1):ds2,2] ###data points
			ys = c(ys,y)
			
			if (!is.null(res.exp)){		
				pars = getPars(res.exp,pars)
				y1 = func( dg$D, pars)
				y1s = c(y1s, y1)
			}
		}
		ds=ds2
	}
	n = length(y1s)    
	#browser()	
	is.na.y1=is.na(y1s)
	is.na.ys = is.na(ys)
	is.nas = is.na.y1 | is.na.ys
	y1s = y1s[!is.nas]
	ys = ys[!is.nas]	
	if (all(ys == -1)){
		rmse= NA
		r2 = NA
	}else{
		if (funcnam.cancer=="medSurvExpV"){
			ra = range(ys)
			le = ra[2]-ra[1]
			y1s = y1s/le
			ys = ys/le
	#		print(y1s)
	#		print(ys)
		}
		rmse = sqrt( sum( (y1s - ys )^2 )/n )
		r2 = sqrt( cor( y1s, ys ))
	}
	
	#browser()
	list(r2=r2,rmse=rmse)
}


DrawCancerDelay <- function(funcnam,surv.data, data.gen, res.exp, delay, col = "green", add = F){
	cex.p = 0.7
	cex.l = 0.8
	dg = data.gen[[1]]
	delay = delay * 1/52.1786
	pars = dg$basic.pars

	pars = getPars(res.exp, pars)
	x = seq(1/(2^3), 10,by=1/(2^3))
	if ("medSurvExpVDiff" == funcnam){
		xl ="Med surv delay impact"	
		y = medSurvExpVDiff(x, pars, delay)
#	print(y)
	}else{
		xl ="met risk increase"	
		y = PmetExpVCondDiff(x, pars, delay)
		y = y
		#print(y)
	}
	
	if (add)
		points( x =  x, y = y, cex=cex.p, ylab=xl, xlab = "Tumor diameter (cm)", col = col,  type="l")
	else
		plot( x =  x, y = y, cex=cex.p, ylab=xl, xlab = "Tumor diameter (cm)", col = col,  type="l", ylim= c(0,max(y,na.rm = T)))

	
}

DrawCancerImmunoth <- function(funcnam,surv.data, data.gen, res.exp, med.bott.increase, col = "green", add = F){
	cex.p = 0.7
	cex.l = 0.8
	dg = data.gen[[1]]
	
	pars = dg$basic.pars
	pars = getPars(res.exp, pars)
	## median.curr = exp(pars$c1)
	## median.new = med.bott.increase*median.curr
	## pars$c1 = log( median.new ) = log( med.bott.increase*median.curr ) = log( med.bott.increase ) + pars$c1
	c1.new = log(med.bott.increase) + pars$c1
	x = seq(0.1, 10,by=0.1)
	if ("medSurvExpVDiff" == funcnam){
		xl ="Med surv immunoth impact"	
		y1 = medSurvExpV(x, pars)
		pars$c1 = c1.new
		y2 = medSurvExpV(x, pars)
		y = y1 - y2
#	print(y)
	}else{
		xl ="met risk decrease"	
		y1 = PmetExpVCond(x, pars)
		pars$c1 = c1.new
		y2 =PmetExpVCond(x, pars) ### Here, the bottleneck is larger so we should have less mets
		y = y1-y2 # we want a positive number,we compute the decrease
		#y = y*100
		#print(y)
	}
	
	if (add)
		points( x =  x, y = y, cex=cex.p, ylab=xl, xlab = "Tumor diameter (cm)", col = col,  type="l")
	else
		plot( x =  x, y = y, cex=cex.p, ylab=xl, xlab = "Tumor diameter (cm)", col = col,  type="l", ylim= c(0,max(y,na.rm = T)))

	
}



DrawCancerChemoth <- function(funcnam,surv.data, data.gen, res.exp, avg.th.age.increase, col = "green", add = F){
	cex.p = 0.7
	cex.l = 0.8
	dg = data.gen[[1]]
	
	pars = dg$basic.pars
	pars = getPars(res.exp, pars)
	## median.curr = exp(pars$c1)
	## median.new = med.bott.increase*median.curr
	## pars$c1 = log( median.new ) = log( med.bott.increase*median.curr ) = log( med.bott.increase ) + pars$c1
	f.new = pars$f/avg.th.age.increase
	
	x = seq(0.1, 10,by=0.1)
	if ("medSurvExpVDiff" == funcnam){
		xl ="Med surv immunoth impact"	
		y1 = medSurvExpV(x, pars)
		pars$f = f.new
		y2 = medSurvExpV(x, pars)
		y = y1 - y2
#	print(y)
	}else{
		xl ="met risk decrease"	
		y1 = PmetExpVCond(x, pars)
		pars$f = f.new
		y2 =PmetExpVCond(x, pars) ### Here, the bottleneck is larger so we should have less mets
		y = y1-y2 # we want a positive number,we compute the decrease
		#y = y*100
		#print(y)
	}
	
	if (add)
		points( x =  x, y = y, cex=cex.p, ylab=xl, xlab = "Tumor diameter (cm)", col = col,  type="l")
	else
		plot( x =  x, y = y, cex=cex.p, ylab=xl, xlab = "Tumor diameter (cm)", col = col,  type="l", ylim= c(0,max(y,na.rm = T)))

	
}



DrawCancerFuncnam <- function(surv.data, data.gen, res.exp, funcnam.cancer, cancer,model){
	
	cex.p = 0.7
	cex.l = 0.8
	dg = data.gen[[1]]
	
	pars = dg$basic.pars

	pars = getPars(res.exp, pars)
	
	tit = paste(cancer,",",model,sep="")
	
	if( "Pcure"%in%funcnam.cancer){

		xl ="Fraction mets removed"
		x = seq(0.1, 10,by=0.1)
		plot( x =  x, y = EPcure( x, pars), cex=cex.p, ylab=xl, xlab = "Tumor diameter (cm)", col = "forestgreen", ylim = c(0,1), type="l")

		title(tit)
		return()
	}
	
	if( "Prem"%in%funcnam.cancer){

		xl ="Probability to remove met"
		x = seq(0.1, 10,by=0.1)
		plot( x =  x, y = Prem( x, pars), cex=cex.p, ylab=xl, xlab = "Met age (years)", col = "forestgreen", ylim = c(0,1), type="l")
		title(tit)
		return()
	}
	
	if( "DBottleneck"%in%funcnam.cancer){
		xl ="Bottleneck density"
		x = seq(0.1, 40,by=0.1)
		if (distribution == "lnorm"){
			plot( x =  x, y = dlnorm( x, meanlog = pars$c1, sdlog =pars$c2), cex=cex.p, ylab=xl, xlab = "b", col = "forestgreen", type="l")
			
		}
		title(tit)
		return()
	}
	
	if( "WTDistr"%in%funcnam.cancer){
		xl ="Waiting time distribution"
		ts = seq(0.1, 15,by=0.1)
		ds = c(0.1, seq(1, 10))
		lambdac = ENMets(x, pars, nocuring=F)
		lambda = ENMets(x, pars, nocuring=T)
		wt.data = NULL
		for (diam in ds){
			wtd = 
			wt.data = cbind(wt.data, )
		}

		y=ENMets(x, pars, nocuring=T)

		if (distribution == "lnorm"){
			plot( x =  x, y = y, cex=cex.p, ylab=xl, xlab = "Tumor diameter (cm)", col = "forestgreen",  type="l")
			lines(x = x, y = yc, col ="magenta")
		}
		title(tit)
		legend("topright",c("without curing","with curing"),fill=c("forestgreen","magenta"))
		return()
	}
	
	if( "ENMets"%in%funcnam.cancer){
		xl ="E[# mets]"
		x = seq(0.1, 10,by=0.1)
		
		yc = ENMets(x, pars, nocuring=F)

		y=ENMets(x, pars, nocuring=T)

		if (distribution == "lnorm"){
			plot( x =  x, y = y, cex=cex.p, ylab=xl, xlab = "Tumor diameter (cm)", col = "forestgreen",  type="l")
			lines(x = x, y = yc, col ="magenta")
		}
		title(tit)
		legend("topright",c("at diagnosis","after curing"),fill=c("forestgreen","magenta"))
		return()
	}
	
	
	if( "MedNMets"%in%funcnam.cancer){
		xl ="Median # mets"
		x = seq(1/2^3, 10,by=1/2^2)
		
		yc = MedNMets(x, pars, nocuring=F)
		y=MedNMets(x, pars, nocuring=T)
#print(y)
		if (distribution == "lnorm"){
			plot( x =  x, y = y, cex=cex.p, ylab=xl, xlab = "Tumor diameter (cm)", col = "forestgreen",  type="l")
			lines(x = x, y = yc, col ="magenta")
		}
		title(tit)
		legend("topright",c("at diagnosis","after chemo"),fill=(c("forestgreen","magenta")))
		return()
	}
	
	if( "DistrNMets.curing"%in%funcnam.cancer){
			
		xl ="P(# mets <= k) after chemo"
		
		maxk = 1000
		
		ds = c(0.5, 1 , 2, 4, 8)
		cols = c("orange","red","firebrick","violet","darkblue")
		mtpl = NULL
		
		for (d in ds){
			distnm = DistrNMets( d, maxk, pars, nocuring = F, delay = 0 )
			
			mtpl = cbind(mtpl, distnm)
		}
		colnames(mtpl)=ds
		matplot(x = seq(0, maxk), y = mtpl, lty = 1, col = cols, 
        xlab = "k", ylab = xl, xlim = c(0, maxk), type = "b")
        legend("topright", legend=ds, fill=cols, title = "d")
		title(tit)
		
		return()
	}
	
	if( "DistrNMets.nocuring"%in%funcnam.cancer){
			
		xl ="P(# mets <= k) at diagnosis"
		maxk = 1000
		
		ds = c(0.5, 1 , 2, 4, 8)
		cols = c("orange","red","firebrick","violet","darkblue")
		mtpl = NULL
		
		for (d in ds){
			distnm = DistrNMets( d, maxk, pars, nocuring = T, delay = 0 )
			mtpl = cbind(mtpl, distrnm)
		}
		colnames(mtpl)=ds
		matplot(x = seq(0, maxk), y = mtpl, lty = 1, col = cols, 
        xlab = "k", ylab = xl, xlim = c(0, maxk), type = "b")
        legend("topright", legend=ds, fill=cols, title = "d")
		title(tit)
		
		return()
	}
	
	le= length(data.gen)
	ds=0
	for (i in (1:le)){
		dg = data.gen[[i]]
		func = dg$func
		funcnam=dg$funcnam
		#pars = dg$basic.pars
		x=dg$D
		ds2=ds+length(x)
		
		if (funcnam%in% funcnam.cancer){
			
			y = surv.data[(ds+1):ds2,2]
			tit = tit
			
			lcol = "red"
		#	if (funcnam.cancer %in% toinclude)
		#		lcol = "blue"
			switch(funcnam, 
				PmetExpVCond={
					xl =  "Metastasis probability"
					plot( x = x, y = y, cex=cex.p, ylab=xl, xlab = "Tumor diameter (cm)", col = "black", ylim = c(0,1),xlim=c(0.1,10))},
				PmetExpVCondc={
					points( x = x, y = y, cex=1, ylab=xl, xlab = "Tumor diameter (cm)", col = "gray", ylim = c(0,1),xlim=c(0.1,10), pch=15)},
				PmetExpVConda={	
					points( x = x, y = y, cex=1, ylab=xl, pch=18, col="gray")},
				obsMetProbExpV={
					xl = "Observed metastasis rate"
					plot( x = x, y = y, cex=cex.p, ylab=xl, xlab = "Tumor diameter (cm)", col = "black", ylim = c(0,1),xlim=c(0.1,10))},
				medSurvExpVObsMet={
					xl = "Med survival (cnc.dth.) w mets"
					if (all(y==-1)){
						lcol='forestgreen'
					}
					plot( x = x, y = y, cex=cex.p, ylab=xl, xlab = "Tumor diameter (cm)", col = "black",ylim=c(0,max(surv.data[,2] )),xlim=c(0.1,10))	
				}
			)
			title(tit)
			xs = seq(0.1,10,by=0.1)
			preds = func( xs, pars)
			lines (x = xs, y = preds, col=lcol)
						
		}
		
		ds=ds2
	}
	

}

DrawCancersFuncnam<-function(funcnam.cancer, model	){
	
	pt.size=0.7
		#hp=7
	hp = 7.5
	wp=8
	print(paste(fitFolder,"/Fit_Func_",funcnam.cancer,"_",model,".pdf",sep=""))
	pdf(paste(fitFolder,"/Fit_Func_",funcnam.cancer,"_",model,".pdf",sep=""), height=hp, width=wp)
	if (funcnam.cancer=="PmetExpVCond")
		funcnam.cancer = c("PmetExpVCond","PmetExpVConda","PmetExpVCondc")
	par(mfrow=c(4,4))
	par(ps=11, mar=c(2.6, 2.7,2.4,0.2), mgp=c(1.6,0.4,0))
	for (rfile in res.files){
	#rfile = res.files[1]
		load(rfile)

		gd =getData(cancer, include= toinclude)
		surv.data = gd$surv.data
		data.gen=gd$data.gen
		DrawCancerFuncnam(surv.data, data.gen, res.exp, funcnam.cancer, cancer,model)
					
	}
	dev.off()
}

DrawCancersDelay<-function(funcnam.cancer, model, delays	){
	cols = sample( colors(), length(delays))
	cols = rep("green", length(delays))
	cols= c("forestgreen","blue","orange")
	names(cols) = delays
	pt.size=0.7
	#hp=7
	hp = 7.5
	wp=8
	
	print(paste(fitFolder,"/Fit_Delay_",funcnam.cancer,"_",model, ".pdf",sep=""))
	pdf(paste(fitFolder,"/Fit_Delay_",funcnam.cancer,"_",model,".pdf",sep=""), height=hp, width=wp)
	#browser()
	par(mfrow=c(4,4))
	par(ps=11, mar=c(2.6, 2.7,2.4,0.2), mgp=c(1.6,0.4,0))
	for (rfile in res.files){
	#rfile = res.files[1]
			load(rfile)
			
			gd =getData(cancer, include= toinclude)
			surv.data = gd$surv.data
			data.gen=gd$data.gen
			for (delay in delays){
				DrawCancerDelay(funcnam.cancer,surv.data, data.gen, res.exp, delay, col = cols[names(cols)==delay], add = (delay != delays[1]))		
			}
			title(cancer)
			legend("topright",legend = paste("delay",delays, "weeks"), fill=cols)
	}	
	dev.off()	
}
	
	


