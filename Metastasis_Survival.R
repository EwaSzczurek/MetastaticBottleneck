

analyzeSurvival <- function(inc12, folder, type){
	cat("\n2) Generating survival data with further filtering from the initial....\n")
	########################################################
	### LIMITING
	
	
	
	### SURVIVAL MONTHS FLAG 
	### 0   Complete dates are available and there are 0 days of survival 
	#1  Complete dates are available and there are more than 0 days of survival 
	#2 Incomplete dates are available and there could be zero days of follow-up 
	#3 Incomplete dates are available and there cannot be zero days of follow-up 
	#9 Unknown 
	### Limit to patients for the survival analysis: patients with known survival time 
	inc.surv.all = inc12[ inc12$SurvivalMnths != 9999, ]
	inc.surv.all = inc.surv.all[inc.surv.all$SurvivalMnthsFlag == 1,]	

	cat( nrow(inc.surv.all), "of patients for quantile survival fitting with known survival time\n")	
	
	#### Get median survival of all patients for barplotting
	attach(inc.surv.all,  warn.conflicts = FALSE)
	agg.survival.all = aggregate(inc.surv.all[,c("SurvivalMnths")], by = list(DXyear), FUN=median)
	agg.survival.all = t(agg.survival.all)
	colnames(agg.survival.all)= agg.survival.all[1,]
	agg.survival.all = agg.survival.all[-1,]
	detach(inc.surv.all)
	
	#### Get median survival of alive patients for barplotting
	inc.surv.alive  = inc.surv.all[(inc.surv.all$DTH_CLASS ==0)&(inc.surv.all$OTHER_DTH_CLASS ==0),]
	attach(inc.surv.alive,  warn.conflicts = FALSE)
	agg.survival.alive = aggregate(inc.surv.alive[,"SurvivalMnths"], by = list(DXyear), FUN=median)
	agg.survival.alive = t(agg.survival.alive)
	colnames(agg.survival.alive)= agg.survival.alive[1,]
	agg.survival.alive = agg.survival.alive[-1,]
	detach(inc.surv.alive)
	
	#### Get median survival of other-dead patients for barplotting
	inc.surv.other  = inc.surv.all[(inc.surv.all$DTH_CLASS ==0)&(inc.surv.all$OTHER_DTH_CLASS ==1),]
	attach(inc.surv.other,  warn.conflicts = FALSE)
	agg.survival.other = aggregate(inc.surv.other[,"SurvivalMnths"], by = list(DXyear), FUN=median)
	agg.survival.other = t(agg.survival.other)
	colnames(agg.survival.other)= agg.survival.other[1,]
	agg.survival.other = agg.survival.other[-1,]
	detach(inc.surv.other)	

	### Limit to patients with surgery
	#inc.surv = inc.surv.all[inc.surv.all$NO_SURG == 0,]
	inc.surv = inc.surv.all	
	
	### Limit to patients who died of cancer 
	inc.surv = inc.surv[inc.surv$DTH_CLASS ==1, ]
	cat( nrow(inc.surv), "of patients for quantile survival fitting with known survival time, with surgery, and who died of cancer\n")	
	
	#### CHECK IF RIGHT TIMING!!
	#### Get median survival of cancer-dead patients for barplotting
	attach(inc.surv)
	agg.survival.cancer = aggregate(inc.surv[,"SurvivalMnths"], by = list(DXyear), FUN=median)
	agg.survival.cancer = t(agg.survival.cancer)
	colnames(agg.survival.cancer)= agg.survival.cancer[1,]
	agg.survival.cancer = agg.survival.cancer[-1,]
	detach(inc.surv)
	
	inc.surv.cancer = inc.surv
	
	### Summarize the counts after filtering
	counts.3 =summarizeGrouping( "Site",inc.surv)
	drawCounts("CancerSubtypesAfterSurvivalFiltering",counts.3)
	
	inc.surv.met = inc.surv[ inc.surv$metExamined, ]
	cat( nrow(inc.surv.met), "of patients with known survival time, with surgery, and who died of cancer and who had mets examined\n")
	
	inc.surv.met = inc.surv.met[ inc.surv.met$met==1, ]
	cat( nrow(inc.surv.met), "of patients with known survival time, with surgery, and who died of cancer and who had mets examined and observed\n")
	counts.5 =summarizeGrouping( "Site",inc.surv.met)
	drawCounts("CancerSubtypesAfterSurvivalFilteringWithMetsObs",counts.5)
	
	########################################################
	### Divide into validation and training by sampling and into early and late by mid.year
	 # indeces = 1:nrow(inc.surv)
	 # tr = sample( indeces, round(nrow(inc.surv)/2))
	 # #te = indeces[ !indeces%in%tr]
	 # inc.surv$testtrain = "validation"
	 # inc.surv$testtrain[indeces%in%tr] = "training"
	# # ########ACHTUNG THIS NEEDS RECOMPUTINNG; TOO SMALL BIN SIZES
	 # inc.surv.testtrain = inc.surv[order(inc.surv$testtrain),]
	 # summ.12.surv.testtrain = sizeVsSummary(inc.surv.testtrain, "SurvivalMnths", c("sizeBin", "Site", "testtrain"), measure = "med")
	# print("Summarized survival test train")	

	
	load("Cancers.RData")
	### Making quantile summaries
	quants  = seq(from = 0.35, to = 0.8, by = 0.03)
	quants.ca  = seq(from = 0.35, to = 0.65, by =0.03)
	runSummQuants<- function(qnt, ds){
		sizeVsSummary(ds, "SurvivalMnths", c("sizeBin", "Site"), measure= "quantile",qnt=qnt)
	}
	quant.surv = lapply(quants, runSummQuants, inc.surv)
	names(quant.surv) = quants
	
	med.met.obs.data = sizeVsSummary(inc.surv.met, "SurvivalMnths", c("sizeBin", "Site"), measure= "med")
	
	print("Summarized survival")
	
	prev.qsd.names = names(quant.surv.data)
	quant.surv.data[[length(quant.surv.data)+1]] = quant.surv

	names(quant.surv.data)=c(prev.qsd.names,type)

	#med.surv.test.train.data = rbind(med.surv.test.train.data, summ.12.surv.testtrain)
	save(pmet.data.oc, pmet.data.oa, binsite.cda, obs.pmet.data, quant.surv.data, med.met.obs.data, file= "Cancers.RData" )


	
	inc.surv$sizeBin = binByTumorSizes(inc.surv, "tumorsize")
	inc.surv$SurvivalYrs = inc.surv$SurvivalMnths/12
	hp = 3.2
	wp = 5
	p2.66666 = ggplot(inc.surv, aes(as.factor(Site), SurvivalYrs)) + 
		 	geom_boxplot( size=pt.size, outlier.size = pt.size)+
			# geom_jitter( aes(shape = impurity), size = pt.size, alpha = 0.4)+
			theme_bw()+
			xlab("Cancer type")+
			ylab("Survival distribution (years)")+
			theme(	text = element_text(size=9), legend.position="bottom",plot.title=element_text(vjust=1, face="bold"),axis.text.x = element_text(angle=90, vjust=1))+ guides(col = guide_legend(nrow = 2))+
    ggtitle("Survival")	
	 pdf(paste(folder,type, "AvgSurvivalMnths.pdf",sep=""), height=hp, width=6)
	 par(ps=9)
	 grid.arrange(p2.66666)
	 dev.off()	


 # p2.666 = ggplot(summ.12.surv.testtrain, aes(x= sizeBin, y= SurvivalMnths, group = Site)) + 
	    # # geom_errorbar(aes(ymin= SurvivalMnths-se, ymax= SurvivalMnths +se), width=.5) +
	     # geom_line(aes(col = Site)) +
	     # theme_bw()+
	     # facet_wrap(~ testtrain)+
	     # xlab("Tumor  size (cm)")+
	     # ylab("Median survival (months)")+
	     # geom_point(aes( col = Site))+# guides(col = guide_legend(nrow = 2))+
	     # #theme(	text = element_text(size=9), legend.text=element_text(size=3), legend.position="bottom", plot.title=element_text(vjust=1,  face = "bold"), legend.text=element_text(size=4))
		 # theme(	text = element_text(size=9), legend.position="bottom",plot.title=element_text(vjust=1,  face = "bold"),legend.text=element_text(size=4))+ guides(col = guide_legend(nrow = 2))
	
	 # pdf(paste(folder,type, "AvgSurvivalMnths_1988_",last.year,"_Site_ValidationTraining_new.eps",sep=""), horizontal = FALSE, onefile = FALSE, paper = "special", height=hp, width=wp)
	 # par(ps=9)
	 # grid.arrange(p2.666)
	 # dev.off()	
	
 pdf(paste(folder,type, "ChangesCancerSurvivalOverDiagYears.pdf",sep=""),  height=5,width=12)

	
	p.aggsurvcancer = ggplot(inc.surv.cancer, aes(as.factor(DXyear), SurvivalMnths)) + 
		 	geom_boxplot( size=pt.size, outlier.size = pt.size, aes(fill= cancertype))+
			# geom_jitter( aes(shape = impurity), size = pt.size, alpha = 0.4)+
			theme_bw()+
			xlab("Diagnosis Year")+
			ylab("Survival (months)")+
			theme(	text = element_text(size=9), legend.position="bottom",plot.title=element_text(vjust=1,  face="bold"),axis.text.x = element_text(angle=90, vjust=1),legend.text=element_text(size=5))+ guides(col = guide_legend(nrow = 2))+
    ggtitle("Survival patients dead of cancer")
    

    
		#multiplot(p.aggsurvall, p.aggsurvcancer ,p.aggsurvother, p.aggsurvalive,cols=2)
	#grid.arrange(p.aggsurvdead+title("Dead patients"))
	grid.arrange(p.aggsurvcancer)
	 dev.off()

	
}