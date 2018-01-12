########################
# Prforms trend analysis - checking whether the clinical variables monotonically decrease/increase with tumor diameter
########################

library(trend)


doTrendSigniff<-function(r.valfitsAll){


# we now compute pvalues for different cancers and datatypes
datat = c( "Met incidence" , "Met probability","0.5-th time to death quantile", "Med time to death w. mets"   )
alltests=expand.grid(cancers, datat)

testTrends<-function( cd){
	data = r.valfitsAll[ r.valfitsAll$axtit == cd[2] & as.character( r.valfitsAll$cancer ) == cd[1], ]
	data = data[, c(1,2)]
	nas =( is.na(data[,2]))
	data = data[!nas,, drop=FALSE]
	if(!is.null(data) & nrow(data)>2){
		ret=mk.test(data[,2], alternative="two.sided")$p.value
	}else
		ret=NA
	ret
}

pvals=unlist( apply(alltests, 1, testTrends) )
alltests$p.values=pvals
alltests$p.values.bh = p.adjust( pvals, method="BH")
alltests$signifficance.adj.bh = alltests$p.values.bh<0.05
 
 ## Reformat
 tests.pvals = matrix(alltests$p.values, nrow=14, byrow=F)
 tests.pvals.bh = matrix(alltests$p.values.bh, nrow=14, byrow=F)
 tests.ppvals = cbind(tests.pvals.bh,tests.pvals)
 rownames(tests.ppvals) = cancers
 tests.ppvals=formatC(tests.ppvals, format = "e", digits = 2)
 write.table( tests.ppvals, file="trendPvalues.txt", sep=" & ", quote=F) ## outpput to latex table in paper..
 }
 
 drawTrends<-function(r.valfitsAll){
# Now we draw it

nas = is.na(r.valfitsAll$y)
r.valfitsAll= r.valfitsAll[!nas,]

minus1 = (r.valfitsAll$y==-1)
r.valfitsAll= r.valfitsAll[!minus1,]

n = 14
cols = gg_color_hue(n)
cols[10]="gray30"


nas = is.na(r.valfitsAll$y)
r.valfitsAll= r.valfitsAll[!nas,]
r.valfitsAll$Cancer = factor( replaceCancers(( r.valfitsAll$cancer)),levels=replaceCancers(levels( r.valfitsAll$cancer)))  
dt.medTD = r.valfitsAll[r.valfitsAll$axtit%in%c("0.5-th time to death quantile","Med time to death w. mets"),]
dt.medTD$axtit[dt.medTD$axtit=="0.5-th time to death quantile"] = "All patients"
dt.medTD$axtit[dt.medTD$axtit=="Med time to death w. mets"] = "Patients with mets"
dt.medTD$y[dt.medTD$y == -1]<-NA
dt.metP = r.valfitsAll[r.valfitsAll$axtit%in%c("Met probability","Met incidence"),]
dt.metP$axtit[dt.metP$axtit=="Met incidence"] = "Met detection"
dt.metP$axtit[dt.metP$axtit=="Met probability"] = "Cancer death"
dt.metP$axtit = factor(dt.metP$axtit, levels=c("Met detection", "Cancer death"))
medTD <- ggplot(dt.medTD, aes(x=diameter, y = y, color= Cancer, group=Cancer)) +
	theme_classic()+ 
	geom_point(shape=1, size = 0.7)+#,color= Cancer)+
	geom_line()+
	scale_color_manual( values=cols) + 
	#scale_y_continuous(name="Median time to death")+#, limits = c(0,7)
	scale_y_continuous(name="")+
	scale_x_continuous(name="Diameter (cm)", breaks=c(0,2,4,6,8, 10), labels=c("0","2","4","6","8","10"))+
	facet_grid(. ~ axtit, scales = "fixed", space = "fixed")+
	theme(axis.title=element_text(size=9), 
	axis.text.y = element_text(size=8), 
	plot.margin = 	margin(1, 1, 1, 1),
	axis.text.x = element_text(size=8), 
	panel.background = element_rect(fill="white", color = "gray60"),  
	strip.text.y = element_text(size=8), 
	strip.text.x = element_text(size=8 , lineheight=0.8), 
	strip.background = element_rect(colour="gray60", fill="white",   size=0.5, linetype="solid"), 
	plot.title = element_blank(), 
	legend.position="none")

cancerdth<- ggplot(dt.metP, aes(x=diameter, y = y, color= Cancer, group=Cancer)) +
	theme_classic()+ 
	geom_point(shape=1, size = 0.7)+#,color= Cancer)+
	geom_line()+
	scale_color_manual( values=cols) +  
	scale_y_continuous(name="")+#, limits = c(0,7)
	scale_x_continuous(name="Diameter (cm)", breaks=c(0,2,4,6,8, 10), labels=c("0","2","4","6","8","10"))+
	facet_grid(.~ axtit, scales = "fixed", space = "fixed")+
	 theme(axis.title=element_text(size=9), 
	 axis.text.y = element_text(size=8), 
	 plot.margin = 	margin(1, 1, 1, 1),
	axis.text.x = element_text(size=8), 
	panel.background = element_rect(fill="white", color = "gray60"),  
	strip.text.y = element_text(size=8), 
	strip.text.x = element_text(size=8 , lineheight=0.8), 
	strip.background = element_rect(colour="gray60", fill="white",   size=0.5, linetype="solid"),  
	plot.title = element_blank(), 
	legend.position="none")
	
	
 pdf("Figure1_trend.pdf", height=1.5, width=5.3)
  multiplot( cancerdth,  medTD, layout=matrix( 
 c(  rep( c( rep(1,30), rep(2,30)), 1)
   ), byrow=T, nrow=1))
 dev.off()


}
