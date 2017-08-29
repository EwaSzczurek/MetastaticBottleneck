analyzeDeathCause<-function(inc12, folder, type){
	cat("\n1) Generating death cause data with further filtering from the initial....\n")	
	
	### Limit to patients with known survival time for the Kaplan Meier estimation of cancer death
	inc.da = inc12[ inc12$SurvivalMnths != 9999, ]
	inc.da = inc.da[inc.da$SurvivalMnthsFlag == 1,]
	cat( nrow(inc.da), "of patients for death cause analysis with known survival time\n")

	####################
	### Limiting to known death status for death cause analysis
	###
	
	###SEER CAUSE-SPECIFIC DEATH CLASSIFICATION ("DTH_CLASS")
	#SEER*Stat Name: SEER cause-specific death classification Code 
	#0 Alive or dead of other cause 
	#1 Dead 
	#9 N/A not first tumor 
	
	#SEER OTHER CAUSE OF DEATH CLASSIFICATION 
	#SEER*Stat Name: SEER other cause of death classification 
	#0 Alive or dead due to cancer 
	#1 Dead 
	#9 N/A not first tumor 
	inc.da = inc.da[, c("DXyear","DTH_CLASS","OTHER_DTH_CLASS","SurvivalMnths","SurvivalMnthsFlag","Site","tumorsize","NO_SURG","SURGPRIM","RAD_SURG")]
	colnames(inc.da) = c("DXyear","Cancer","Other","Survival","SurvivalFlag","Site", "tumorsize","NO_SURG","SURGPRIM","RAD_SURG")
	inc.da$Alive = (inc.da$Cancer!=1)&(inc.da$Other!=1)
	attach(inc.da,  warn.conflicts = FALSE)
	agg.deathcause = aggregate( inc.da[ ,c("Cancer","Other","Alive")], by = list(DXyear), FUN=sum )
	agg.deathcause = t(agg.deathcause)
	colnames(agg.deathcause)=agg.deathcause[1,]
	agg.deathcause = agg.deathcause[-1,]
	class(agg.deathcause)<-"numeric"	
	detach(inc.da)
	
	pdf(paste(folder,type, "DeathCausesOverYears.pdf",sep=""), height=3,width=4)
	par(ps=9, mar=c(2.6, 2.5,1.4,0.1), mgp=c(1.7,0.4,0))
	cols=c("darkblue","red","green")
	barplot(agg.deathcause,col=cols,las=2,main="Death causes",ylab="# cases")
	legend("topleft", legend=rownames(agg.deathcause),fill=cols,cex = 0.7)
	dev.off()
	
	attach(inc.da,  warn.conflicts = FALSE)
	agg.deathcause = aggregate( inc.da[ ,c("Cancer","Other","Alive")], by = list(Site), FUN=sum )
	agg.deathcause = t(agg.deathcause)
	colnames(agg.deathcause)=agg.deathcause[1,]
	agg.deathcause = agg.deathcause[-1,]
	class(agg.deathcause)<-"numeric"	
	detach(inc.da)
	
	pdf(paste(folder,type, "DeathCausesOverSite.pdf",sep=""), height=4,width=4)
	par(ps=9, mar=c(7, 2.5,0.1,0.1), mgp=c(1.7,0.4,0))
	cols=c("darkblue","red","green")
	barplot(agg.deathcause,col=cols,las=2,ylab="# cases")
	legend("topleft", legend=rownames(agg.deathcause),fill=cols,cex = 0.7)
	dev.off()

	su = apply(data.frame(agg.deathcause), 2, sum)
	agg.deathcause.r = t(apply(agg.deathcause, 1, function(x){x/su}))

	pdf(paste(folder,type, "DeathCausesOverSiteRescaled.pdf",sep=""), height=4,width=4)
	par(ps=9, mar=c(7, 2.5,0.1,0.1), mgp=c(1.7,0.4,0))
	cols=c("darkblue","red","green")
	barplot(agg.deathcause.r,col=cols,las=2,ylab="# cases")
	legend("topleft", legend=rownames(agg.deathcause.r),fill=cols,cex = 0.7)
	dev.off()
	# returns string w/o leading whitespace
trim.leading <- function (x)  sub("^\\s+", "", x)
# returns string w/o trailing whitespace
trim.trailing <- function (x) sub("\\s+$", "", x)

# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
load("Cancers.RData")
	## we treat the cancer death as the death event. Alive and dead of other reasons are treated as censoring
	inc.da$SurvObj <- with(inc.da, Surv(Survival, Cancer == 1))
	inc.da$sizeBin = binByTumorSizes(inc.da, "tumorsize")
	km.by.d <- survfit(SurvObj ~ sizeBin+Site, data = inc.da, conf.type = "log-log")
	binsite = names(km.by.d$strata)
	binsite.cda=NULL
	sites = NULL
	for (i in 1:length(binsite)){
		txt=binsite[i]
		ss.txt = strsplit(txt,",")[[1]]
		sb = as.numeric( strsplit(ss.txt[1],"=")[[1]][2])
		si = ( strsplit(ss.txt[2],"=")[[1]][2])
		si = trim(si)
		suri = km.by.d[i]$surv
		cdai = 1-suri[length(suri)] # we take 1-survival at the last timepoint as the cancer/all estimate
		len = length(suri)
		
		if (len>=binSize){
			sites = c(sites, si)
			binsite.cda<- rbind(binsite.cda, c(sb, len, cdai))
		}
	}
	
	binsite.cda<-data.frame(binsite.cda)
	colnames(binsite.cda)=c("SizeBin","Patient.no", "CancerVsAll")
	binsite.cda$Site = sites
	
	pdf(paste(folder, "KM_before.pdf",sep=""), height=3,width=4)
	par(ps=9, mar=c(2.6, 2.5,1.4,0.1), mgp=c(1.7,0.4,0))
	km.data <- list()
	for (sit in unique(inc.da$Site)){
		#print(sit)
		inc.surv.surg.s = inc.da[inc.da$Site==sit,]
		km.by.d <- survfit(SurvObj ~ sizeBin, data = inc.surv.surg.s, conf.type = "log-log")
		col.n = length(unique(inc.surv.surg.s$sizeBin))
		plot(km.by.d, col = topo.colors(col.n) , main=sit, xlab="Years", ylab="Survival", cex=0.9)
		legend("bottomright", legend=unique(inc.surv.surg.s$sizeBin), fill=topo.colors(col.n), cex=0.4)
		km.data[[sit]] <- km.by.d
		}
	dev.off()

	inc.da.surg.oc = inc.da
	inc.da.surg.oa = inc.da
	marked.other = (inc.da$Other==1)
	
	## Making a dataset with patients who died of other reasons dead due to cancer  
	inc.da.surg.oc$Cancer[marked.other] = 1
	
	## Making a dataset with patients who died of other reasons alive 
	inc.da.surg.oa$Other = 0
		
	
	#### Summarize the counts after filtering
	counts.2 =summarizeGrouping( "Site",inc.da.surg.oc)
	drawCounts("CancerSubtypesAfterDeathCauseFiltering",counts.2)
	
	
	summ.12.death.site.surg.oc= sizeVsSummary(inc.da.surg.oc, "Cancer", c("sizeBin","Site") ,measure = "mean")
	summ.12.death.site.surg.oa= sizeVsSummary(inc.da.surg.oa, "Cancer", c("sizeBin","Site") ,measure = "mean")
	print("Summarized death")
	
	
	pmet.data.oc = rbind(pmet.data.oc, summ.12.death.site.surg.oc)
	pmet.data.oa = rbind(pmet.data.oa, summ.12.death.site.surg.oa)
	#cancer.alive = rbind(cancer.alive, summ.12.life.site.surg) ### the prevalence of alive cases
save(pmet.data.oc, pmet.data.oa, binsite.cda , obs.pmet.data, quant.surv.data, med.met.obs.data, file= "Cancers.RData" )
	
	pdf(paste(folder,type, "AvgDeathCause.pdf",sep=""), height=10, width=7)
	par(ps=9, mar=c(2, 1,0.4,0.1), mgp=c(1.7,0.4,0))
			    
	p55.1 = ggplot(summ.12.death.site.surg.oc, aes(x= sizeBin, y= (Cancer), group=Site)) + 
	  #  geom_errorbar(aes(ymin= met-se, ymax= met +se), width=.5) +
	    geom_line(aes(col= Site)) +
	    theme_bw()+
	    ylim(0, 1)+
	    xlab("Tumor  size (cm)")+
	    ylab("% of patients dying of cancer")+
	    geom_point(aes(col= Site))+
	    theme(	text = element_text(size=9), legend.position="bottom",plot.title=element_text(vjust=1,  face = "bold"),legend.text=element_text(size=4))+ guides(col = guide_legend(nrow = 4))+
	    ggtitle(paste("Cancer deaths (<=",(censoring.cutoff),", others-->cancer)"))
#browser()
	p55.2 = ggplot(summ.12.death.site.surg.oa, aes(x= sizeBin, y= (Cancer), group=Site)) + 
	  #  geom_errorbar(aes(ymin= met-se, ymax= met +se), width=.5) +
	    geom_line(aes(col= Site)) +
	    theme_bw()+
	    ylim(0, 1)+
	    xlab("Tumor  size (cm)")+
	    ylab("% of patients dying of cancer")+
	    geom_point(aes(col= Site))+
	    theme(	text = element_text(size=9), legend.position="bottom",plot.title=element_text(vjust=1,  face = "bold"),legend.text=element_text(size=4))+ guides(col = guide_legend(nrow = 4))+
	    ggtitle(paste("Cancer deaths (<=",(censoring.cutoff),", others-->alive)"))
	multiplot(p55.1, p55.2,cols=1)
	dev.off()
	
}
