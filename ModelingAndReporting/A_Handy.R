load("Cancers.RData")
source("A_HandyDraw.R")
source("A_HandyFit.R")
source("A_HandySummarize.R")

#### Data access functions
cancers = c("Breast Ca (NST)", "Breast Ca lobular" ,"Ovarial Ca","Endometrial Ca",
  "Esophageal Ca", "Gastric AdenoCa", "Colon AdenoCa" , "Colon AdenoCa muc",  "Rectum AdenoCa",
  "Pancreas AdenoCa",
  "NSCLC", 
  "HNSCC" , 
    "Renal Cell Ca", 
     "Urothelial Ca")
     
newcancers = c("Breast", "Breast lob","Ovarian", "Endometrial","Esophageal", "Gastric", "Colon","Colon muc","Rectal","Pancreatic", "Lung","Head & neck","Renal", "Bladder")

cancerMat = cbind(cancers, newcancers)
# Fix r values
rs.days = c(rep(212,4), #Female tumors
rep(255,5), #GI tumors
144, # Pancreas
166.3, #"NSCLC", 
99, #HNSCC
603, #Renal cell ca
255 #Urothelial
)
rs = rs.days*log(2)/365.24
names(rs)=cancers
q = 11 # per cell per year release rate (for cells on tumor surface)
e = 0.8 #per cell extravasation rate
X.0 = 10^{-4} #initial condition for the tumor growth
const = (4*pi)^{1/3}*3^{2/3}*10^6*(X.0)^{2/3}

toinclude=c("medSurvExpV","obsMetProbExpV") ### these are the functions (model equations), which will be used for model fitting to data: quantile survival and metastasis detection rate



# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#################################	
#### Fit cancers functions
#################################	
MakeCancerFits<-function(cancers){
	sim.no = length(cancers)
	if (!paralelPars){
		registerDoMC(cores = sim.no)
	
		sim.res <- foreach (cancer = cancers , .combine = rbind) %dopar% {
			print(cancer)
  			FitCancer(cancer, draw=TRUE, nam = cancer)
		}
	}else{
		sim.res = NULL
		for (cancer in cancers){
			print(cancer)
			sim.res = rbind(sim.res, FitCancer(cancer, draw=FALSE, nam = cancer))	
		}
	}
	colnames(sim.res)=actual.names
	rownames(sim.res)=cancers	
	save(sim.res, fitFolder, file=paste(fitFolder,"Fits.RData",sep=""))
	write.table(sim.res, paste(fitFolder,"Fits.txt",sep=""), sep="\t", quote = F)
	sim.res
}


###### 
## Get cancers data functions
replaceCancers<- function( cancers ){
	for (i in 1:nrow(cancerMat)){
		c= cancerMat[i,1]
		nc = cancerMat[i,2]
		cancers = gsub( c, nc, cancers, fixed = T)
	}
	cancers
}


getQuants<- function( cancer, quantdata ){
	for (n in names(quantdata)){
		qs = quantdata[[n]]
		if ( cancer %in% qs[[1]]$Site )
			return(qs)
	}
	NULL
}

getData<- function(cancer, include=c("medSurvExpV","medSurvExpVCA","obsMetProbExpV","PmetExpVCond","medSurvExpVObsMet")){
	load("Cancers.RData")
	
	basic.pars = list( r = rs[cancer], tumor.model = "exp", distr = distribution)
	SurvMonthsPMetCond = NULL
	quants  = seq(from = 0.35, to = 0.65, by = 0.03)
	
	makemedpars <- function(quant, quant.data){
		qs = quant.data[[as.character(quant)]]
		qs = qs[qs$Site == cancer, ]
		basic.pars$quant = quant
		list(func = medSurvExpV, basic.pars = basic.pars, D = qs$sizeBin, funcnam="medSurvExpV")
	}
	
	quant.data = quant.surv.data
	quant.data = getQuants(cancer, quant.data)
	basic.pars$quant = NULL
	data.gen = list()
	
	### Constructing data.gen	
	if ("medSurvExpV" %in% include){
		data.gen = lapply(quants, makemedpars, quant.data=quant.data)
	}
	
	if ("PmetExpVCond" %in% include){
		cda = binsite.cda[binsite.cda$Site == cancer, ]
		data.gen[[ length(data.gen)+1]] <- list(func = PmetExpVCond, basic.pars = basic.pars, D = cda$SizeBin,funcnam="PmetExpVCond")
	}
	
	if ("PmetExpVCondc" %in% include){
		pmet.data.oc = pmet.data.oc[pmet.data.oc$Site == cancer, ]	
		data.gen[[ length(data.gen)+1]] <- list(func = PmetExpVCond, basic.pars = basic.pars, D = pmet.data.oc$sizeBin,funcnam="PmetExpVCondc")
	}	
	
	if ("PmetExpVConda" %in% include){
		pmet.data.oa = pmet.data.oa[pmet.data.oa$Site == cancer, ]
		data.gen[[ length(data.gen)+1]] <- list(func = PmetExpVCond, basic.pars = basic.pars, D = pmet.data.oa$sizeBin,funcnam="PmetExpVConda")	
	}
	
	if ("obsMetProbExpV" %in% include){
		obs.pmet.data = obs.pmet.data[obs.pmet.data$Site == cancer, ]	
		data.gen[[length(data.gen)+1]] <- list(func = obsMetProbExpV, basic.pars = basic.pars, D = obs.pmet.data$sizeBin, funcnam="obsMetProbExpV")
	}
	
	if ("medSurvExpVObsMet"%in%include){
		med.met.obs.data = med.met.obs.data[med.met.obs.data$Site == cancer, ]
		if (nrow(med.met.obs.data)==0){
			Ds= seq(1,10)
		}else{
			Ds=med.met.obs.data$sizeBin
		}
		basic.pars$quant = 0.5
		data.gen[[length(data.gen)+1]] <- list(func = medSurvExpVObsMet, basic.pars = basic.pars, D = Ds, funcnam="medSurvExpVObsMet")
	}
	
	makemeds <- function(quant, quant.data){
		qs = quant.data[[as.character(quant)]]
		qs = qs[qs$Site == cancer, ]
		(qs[,7])/12
	}
	makequantsizes <- function(quant, quant.data){
		qs = quant.data[[as.character(quant)]]
		qs = qs[qs$Site == cancer, ]
		qs$sizeBin
	}
 
 ## Preparing SurvMonths
 
	SurvMonths = NULL 
	quantsizes= NULL
	if ("medSurvExpV" %in% include){
		SurvMonths= unlist( lapply(quants, makemeds, quant.data) )
		quantsizes=unlist(lapply(quants, makequantsizes, quant.data))
	}
	
	pmetCond = NULL
	pmetCond.oc = NULL
	pmetCond.oa = NULL
	pmet.data.oc.size = NULL
	pmet.data.oa.size = NULL
	pmetCond.size=NULL

	
	if ("PmetExpVCond" %in% include){
		pmetCond = cda$CancerVsAll
		pmetCond.size = cda$SizeBin	
	}
	
	if ("PmetExpVConda" %in% include){
		pmetCond.oa = pmet.data.oa$Cancer
		pmet.data.oa.size = pmet.data.oa$sizeBin
	}
	
	if ("PmetExpVCondc" %in% include){
		pmetCond.oc = pmet.data.oc$Cancer	
		pmet.data.oc.size = pmet.data.oc$sizeBin
	}
	
	obsMet = NULL
	obs.pmet.data.size = NULL
	if ("obsMetProbExpV" %in% include){
		obsMet = obs.pmet.data$met
		obs.pmet.data.size = obs.pmet.data$sizeBin
	}
	
	med.met.obs =NULL
	med.met.obs.size = NULL
	if ("medSurvExpVObsMet"%in%include){
		if (nrow(med.met.obs.data) == 0){
			med.met.obs = rep(-1,10)
			med.met.obs.size = seq(1,10)
		}else{
			med.met.obs = med.met.obs.data$SurvivalMnths/12
			med.met.obs.size = med.met.obs.data$sizeBin
		}
	}
	
	SurvMonthsPMetCond = c( SurvMonths, pmetCond, pmetCond.oc, pmetCond.oa, obsMet, med.met.obs )

	surv.data = data.frame(Size=c(quantsizes, pmetCond.size, pmet.data.oc.size, pmet.data.oa.size, obs.pmet.data.size, med.met.obs.size), SurvMonths=SurvMonthsPMetCond)
	list(surv.data = surv.data, data.gen = data.gen)
}


##############
### Model functions
########
### Basic functions
########
NT1T2 <- function(T1, T2, rate){
	e*q*const*1/rate*( exp(rate*T2) - exp(rate*T1))
}

NT1T2Exp <- function(T1, T2, r){
	rate = 2/3*r
	NT1T2(T1, T2, rate)
}

FT1T2Exp<- function(T1, T2, r, s1){
	-expm1(-1*s1* NT1T2Exp(T1, T2, r))

}

NT1T2ExpC <- function(T1, T2, Td, f, r ){
	rate = (2/3*r + f)
	NT1T2Exp(T1, T2, r)-exp(-f*Td)*NT1T2(T1, T2, rate)
}

FT1T2ExpC<- function(T1, T2, Td, f, r, s1){
	-expm1(-1*s1* NT1T2ExpC(T1, T2, Td, f, r))
}

s1SS<-function(b){
	exp( -b )
}

tdExp <- function( D, r ){
	1/r *( log(pi/6) + 3*log(D) - log(X.0) )
}
