library("ggplot2")
library(calibrate)
library(gridExtra)
library(survival)
tumor.max = 10
binSize = 20
binVal = 0.5
min.year = 1988
mid.year = 2004
early.year = 1990
late.year = 1990

#censoring.cutoff = 2000
#censoring.cutoff = 1995
#censoring.cutoff = 2005
#censoring.cutoff =2003
#censoring.cutoff =1993
censoring.cutoffs = c(1993, rep(1998,10), rep(2003, 3))
names(censoring.cutoffs)=c("Breast Ca (NST)", "Breast Ca lobular", "Urothelial Ca", "Renal Cell Ca", "NSCLC","HNSCC", "Colon AdenoCa", "Rectum AdenoCa", "Colon AdenoCa muc", "Endometrial Ca", "Ovarial Ca", "Pancreas AdenoCa",  "Esophageal Ca" , "Gastric AdenoCa") 

 
#life.limit=12*10

last.year = 2013
surv.year = 5
pt.size=0.5
hp = 3
wp = 5

gr = "Grouped.Coding.1"
#cancer.lifetime = 10
cols12 = c('DXage','DXyear',"SEQ_NUM", 'grade','RegNodesPos','RegNodesExamined',  'met', 'metExamined', 'NO_SURG', 'DTH_CLASS', 'OTHER_DTH_CLASS', 'SurvivalMnths', 'SurvivalMnthsFlag', 'tumorsize', 'Site','cancertype',"SURGPRIM","RADIATN", "RAD_SURG","SS_SURG")


binByTumorSizes <- function(df, ts){
	
	binSite<- function(Site, df){
		df = df[ df$Site == Site, ]
		cases.l = seq(0,10, by = 0.5)
		cases.l[cases.l == 0]<- 0.1
		approx.ts = sapply(df[,ts], function(x){ min.d = min(abs(cases.l - x)); cases.l[ which(abs(cases.l - x)==min.d)]})
		approx.ts
	}
	sites <- unique(df$Site)
	ms = sapply(sites, binSite, df= df)
		
	unlist(ms)
}
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

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE, measure = "mean", qnt=0.5) {
    require(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm),
          med = median(xx[[col]], na.rm=na.rm),
          quantile = quantile(xx[[col]], c(qnt), na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column 
    if (measure == "mean"){
	    datac <- rename(datac, c("mean" = measurevar))
	}else{
		if(measure=="med"){
			datac <- rename(datac, c("med" = measurevar))
		}
	}
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

sizeVsSummary <- function( df, toSummarize, groupingCols, measure, qnt=0.5){
	df$sizeBin = binByTumorSizes(df, "tumorsize")
#	sizeGroupDf = summarySE(df, "tumorsize", groupingCols , measure = "med")
#	sizeGroups = sizeGroupDf$tumorsize
	
	summ = summarySE(df, toSummarize, groupingCols , measure = measure, qnt=qnt)	
	summ = summ[ (summ$N>=binSize),]
#	summ$sizeGroups = sizeGroups
	summ
	
}

initialize<- function(){
	pmet.data.oc = NULL #cancer/all assuming patients dead of other reasons all would die of cancer
	pmet.data.oa = NULL #cancer/all assuming patients dead of other reasons all would be alive
	obs.pmet.data=NULL
	quant.surv.data = list()
binsite.cda =NULL
	med.met.obs.data=NULL
	save(pmet.data.oc, pmet.data.oa, binsite.cda, obs.pmet.data, quant.surv.data, med.met.obs.data, file= "Cancers.RData" )
}
