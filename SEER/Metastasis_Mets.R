analyzeMets<-function(inc12, folder, type){
	cat("\n3) Generating observed met rate data with further filtering from the initial....\n")
	########################################################
	### Limit to patients for  average dist metastasis fitting 
	

	inc.met = inc12[ inc12$metExamined, ]
	cat( nrow(inc.met), "of patients for average dist metastasis fitting with mets examined\n")
counts.4 =summarizeGrouping( "Site",inc.met)
drawCounts("CancerSubtypesAfterFilteringForMetsExam",counts.4)


	#summ.12.met= summarySE(inc.met, "met", c("sizeBin","Site") )
		
	summ.12.met.site= sizeVsSummary(inc.met, "met", c("sizeBin","Site") ,measure = "mean")
	print("Summarized met site")
	#write.table(summ.12.met.site, paste(folder, "size_met_site.txt", sep=""), sep="\t")
	load("Cancers.RData")
		
	obs.pmet.data=rbind(obs.pmet.data, summ.12.met.site)
	#med.surv.early.late.data = rbind(med.surv.early.late.data,summ.12.surv.earlylate)
	
#	save(pmet.data, obs.pmet.data, med.surv.data, med.surv.test.train.data, med.surv.early.late.data, file= "Cancers.Rdata" )
save(pmet.data.oc, pmet.data.oa, binsite.cda, obs.pmet.data, quant.surv.data, med.met.obs.data, file= "Cancers.RData" )


 
 hp = 4
	pdf(paste(folder,type, "AvgDistMetastasis.pdf",sep=""),  height=hp, width=wp)
	par(ps=9)
	
	p44 = ggplot(summ.12.met.site, aes(x= sizeBin, y= (met), group=Site)) + 
	  #  geom_errorbar(aes(ymin= met-se, ymax= met +se), width=.5) +
	    geom_line(aes(col= Site)) +
	    theme_bw()+
	    ylim(0, 1)+
	    xlab("Tumor  size (cm)")+
	    ylab("% of patients with dist. metastasis")+
	    geom_point(aes(col= Site))+
	    theme(	text = element_text(size=9), legend.position="bottom",plot.title=element_text(vjust=1,  face = "bold"),legend.text=element_text(size=3))+ guides(col = guide_legend(nrow = 4))
	
	grid.arrange(p44)
	dev.off()

	 


}
