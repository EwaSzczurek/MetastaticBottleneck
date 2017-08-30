library(trend)


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}




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
alltests$p.values.bonf = alltests$p.values*50 # no of non na pvalues= 50
alltests$signifficance.adj.bonf = (alltests$p.values*50 < 0.05)
alltests$p.values.bh = p.adjust( pvals, method="BH")
alltests$signifficance.adj.bh = alltests$p.values.bh<0.05
 
 ## Reformat
 tests.pvals = matrix(alltests$p.values, nrow=14, byrow=F)
 tests.pvals.bh = matrix(alltests$p.values.bh, nrow=14, byrow=F)
 tests.pvals.bonf = matrix(alltests$p.values.bonf, nrow=14, byrow=F)
 tests.ppvals = cbind(tests.pvals.bh,tests.pvals.bonf,tests.pvals)
 rownames(tests.ppvals) = cancers
 write.table( tests.ppvals, file="trendPvalues.xls", sep="\t", quote=F)
 
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
dt.medTD$axtit[dt.medTD$axtit=="0.5-th time to death quantile"] = "Overall"
dt.medTD$axtit[dt.medTD$axtit=="Med time to death w. mets"] = "With metastases"
dt.medTD$y[dt.medTD$y == -1]<-NA
dt.metP = r.valfitsAll[r.valfitsAll$axtit%in%c("Met probability","Met incidence"),]
dt.metP$axtit[dt.metP$axtit=="Met incidence"] = "Metastasis detection"
dt.metP$axtit[dt.metP$axtit=="Met probability"] = "Cancer mortality"
dt.metP$axtit = factor(dt.metP$axtit, levels=c("Metastasis detection", "Cancer mortality"))
medTD <- ggplot(dt.medTD, aes(x=diameter, y = y, color= Cancer, group=Cancer)) +
	theme_classic()+ 
	geom_point(shape=1, size = 0.8)+#,color= Cancer)+
	geom_line()+#data = r.valfitsAll[r.valfitsAll$axtit=="0.5-th time to death quantile",], aes(x=diameter, y = y, color= Cancer))+
#	geom_line(color= Cancer)+
	scale_color_manual( values=cols) + 
	scale_y_continuous(name="Median time to death")+#, limits = c(0,7)
	scale_x_continuous(name="Diameter", breaks=c(0,2,4,6,8, 10), labels=c("0","2","4","6","8","10"))+
	facet_grid(. ~ axtit, scales = "fixed", space = "fixed")+
	theme(axis.title=element_text(size=9), 
	axis.text.y = element_text(size=8), 
	plot.margin = 	margin(1, 1, 1, 1),
	axis.text.x = element_text(size=8), 
	panel.background = element_rect(fill="white", color = "gray60"),  
	strip.text.y = element_text(size=8), 
	strip.text.x = element_text(size=8 , lineheight=0.8), 
	strip.background = element_rect(colour="gray60", fill="white",   size=0.5, linetype="solid"), 
#	plot.title = element_text(lineheight = 0.05, hjust = -0.25, vjust = -73, face="bold",size=11), 
	plot.title = element_blank(), 
	legend.position="none")

cancerdth<- ggplot(dt.metP, aes(x=diameter, y = y, color= Cancer, group=Cancer)) +
	theme_classic()+ 
	geom_point(shape=1, size = 0.8)+#,color= Cancer)+
	geom_line()+#data = r.valfitsAll[r.valfitsAll$axtit=="0.5-th time to death quantile",], aes(x=diameter, y = y, color= Cancer))+
#	geom_line(color= Cancer)+
	scale_color_manual( values=cols) +  
	scale_y_continuous(name="Rate")+#, limits = c(0,7)
	scale_x_continuous(name="Diameter", breaks=c(0,2,4,6,8, 10), labels=c("0","2","4","6","8","10"))+
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
	
	
 pdf("Figure1_trend.pdf", height=1.5, width=6.3)
  multiplot( cancerdth,  medTD, layout=matrix( 
 c(  rep( c( rep(1,14), NA, rep(2,14)), 1)
   ), byrow=T, ncol=29))
 dev.off()


 metinc<- ggplot(r.valfitsAll[r.valfitsAll$axtit=="Met incidence",], aes(x=diameter, y = y, color= Cancer, group=Cancer)) +
	theme_classic()+ 
	geom_point(shape=1, size = 0.8)+
	geom_line()+
		scale_color_manual( values=cols) + 
	scale_y_continuous(name="Metastasis incidence rate")+#, limits = c(0,7)
	scale_x_continuous(name="Diameter", breaks=c(0,2,4,6,8, 10), labels=c("0","2","4","6","8","10"))+
		 guides(colour = guide_legend(nrow = 2))+
	 theme(axis.title=element_text(size=9), 
	 axis.text.y = element_text(size=8), 
	 plot.margin = 	margin(3, 5, 1, 0),
	axis.text.x = element_text(size=8), 
	#panel.background = element_rect(fill="white", color = "gray60"),  
	#strip.text.y = element_text(size=8), 
	#strip.text.x = element_text(size=8 , lineheight=0.8), 
	strip.background = element_rect(colour="gray60", fill="white",   size=0.5, linetype="solid"), 
#	plot.title = element_text(lineheight = 0.05, hjust = -0.25, vjust = -73, face="bold",size=11), 
	plot.title = element_blank(), legend.title= element_blank(),
	legend.position = "bottom", legend.justification="right", legend.text=element_text(size=8), legend.key.size = unit(0.25, "cm"),
panel.spacing.y=unit(0.3, "lines"), legend.margin=margin(-1, -1,0,-1))
 pdf("Figure1_trend_legend.pdf", height=1.5, width=6.5)
 metinc
 dev.off()

