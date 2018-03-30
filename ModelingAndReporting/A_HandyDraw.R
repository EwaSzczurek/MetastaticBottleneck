library(devtools)
library(ggplot2)
library(reshape2)
library(scales)
library(RColorBrewer)
library(grid) # for unit()
library(gridExtra)
library(ggfortify)
library(survival)



makeReplacement<- function(cancer){
	cancer = gsub("Endometrial", "Endo-\nmetrial", cancer )
    cancer = gsub("Esophageal", "Eso-\nphageal", cancer )
   cancer = gsub("Head & neck", "Head\n& neck",cancer )
	cancer = gsub("Colon muc", "Colon\nmuc", cancer )
	cancer = gsub("Breast lob", "Breast\nlob", cancer)
	cancer = gsub("Pancreatic", "Pan-\ncreatic", cancer)
}

reverseReplacement<-function(cancer){
		cancer = gsub("-\n", "", cancer )
		cancer	= gsub("\n"," ", cancer)
		cancer
}


Figure2<- function(r.fits, r.val, r.valfitsAll, km.data, appr, par.h){
	red="red3"
	blue="royalblue3"
	orange="darkorange2"


	r.fits$Cancer= factor( replaceCancers(( r.fits$Cancer)),levels=replaceCancers(levels( r.fits$Cancer)))
	r.fits$Cancer = makeReplacement(r.fits$Cancer)
	r.fits = r.fits[!r.fits$Cancer =="Ovarian",]
		
	tot = aggregate(r.fits$NRMSE, by=list( Cancer= r.fits$Cancer), FUN=function(x){sum(x, na.rm = T)/sum(!is.na(x))})
	Cancers = as.character( tot$Cancer[order(tot$x)])
	r.fits$Cancer = factor( r.fits$Cancer, levels = Cancers)

	r.fits$Data = gsub(" time to death quantile", "", r.fits$Data )
	#r.fits$Data = gsub("w. mets", "with metastases", r.fits$Data )
	r.fits$Data = gsub(" rate", "", r.fits$Data )
	r.fits$Data = gsub("0.5-th", "0.50-th", r.fits$Data )
	r.fits$Data[r.fits$Data=="Met incidence"] = "Met detection"
	r.fits$Cancer2=r.fits$Cancer
	r.fits.b = r.fits[!r.fits$Data%in%c("Met detection"),]
	r.fits.a = r.fits[r.fits$Data%in%c("Met detection"),]
 
	r.fits$Data = factor(r.fits$Data, levels = unique(r.fits$Data))
	
	
	
	Atotal1.a<- ggplot(r.fits.a, aes(x = Data, y = NRMSE, fill=Type)) +
  	geom_bar(stat = "identity")+scale_fill_manual( values=c(blue)) + 
 	 facet_grid(Cancer2 ~ ., scales = "fixed", space = "fixed") + 
 	 theme_classic()+
 	 scale_y_continuous(name="Root mean squared error (probability)", breaks=c(0,0.02, 0.04), labels=c("0","0.02","0.04"), limits=c(0, 0.045))+
 	 scale_x_discrete(name="")+
 	 # ggtitle("a")+
 	 # guides(fill=guide_legend(title=NULL))+
 	 theme(axis.title.y=element_text(size=9),
 	 	axis.text.y = element_text(size=8),
		axis.text.x = element_text(angle=60, hjust = 1, size=8), 
		panel.background = element_rect(fill="white", color = "gray60"),
           panel.grid.major = element_blank(), 
           plot.margin = 	margin(5, 3, 0, 0),
			strip.text.y = element_blank(),
			strip.background = element_blank(),
			axis.title.x = element_text(hjust=0, vjust = 14, color="gray30", size = 8),
			#plot.title = element_text(lineheight = 0.05, hjust = -0.25, vjust = -67, face="bold",size=12),
			plot.title= element_blank(),
			 legend.position = "none")

Atotal1.b<- ggplot(r.fits.b, aes(x = Data, y = NRMSE, fill=Type)) +
  	geom_bar(stat = "identity")+scale_fill_manual( values=c(blue,blue)) + 
 	 facet_grid(Cancer2 ~ ., scales = "fixed", space = "fixed") + 
 	 theme_classic()+
 	 scale_y_continuous(name="Root mean squared error (years)")+#, breaks=c(0, 0.5, 1, 1.5), labels=c("0","0.5", "1.0", "1.5"), limits=c(0, 1.52))
 	 scale_x_discrete(name="Time to death quantiles")+
 	 # ggtitle("a")+
 	  #guides(fill=guide_legend(title=NULL))+
 	 theme(axis.title.y=element_text(size=9),
 	 	axis.text.y = element_text(size=8),
		axis.text.x = element_text(angle=60, hjust = 1, size=8), 
		panel.background = element_rect(fill="white", color = "gray60"),
           panel.grid.major = element_blank(), 
           plot.margin = 	margin(5, 0, 16, 7),
			strip.text.y = element_text(angle=0,hjust=0, size=8), 
			strip.background = element_rect(colour="gray60", fill="white",   size=0.5, linetype="solid"), 
			axis.title.x = element_text(hjust=0.05, vjust = 2, color="gray30", size = 8),
			#plot.title = element_text(lineheight = 0.05, hjust = -0.22, vjust = -68, face="bold",size=11), 
			plot.title = element_blank(), 
			legend.position = "none")

 pdf("ExtendedData_Figure2.pdf", height=6.7, width=3.4)
  multiplot(Atotal1.a, Atotal1.b, layout=matrix( 
  c( rep(1,4), rep(2,17) ), byrow=T, ncol=21))
 dev.off()


    r.valfitsAll$cancer = factor( replaceCancers(( r.valfitsAll$cancer)),levels=replaceCancers(levels( r.valfitsAll$cancer)))
	r.valfitsAll$cancer = makeReplacement(r.valfitsAll$cancer)
	r.valfitsAll = r.valfitsAll[! r.valfitsAll$cancer =="Ovarian",]
    r.valfitsAll$cancer = factor( r.valfitsAll$cancer)                                	
	r.valfitsAll$axtit[r.valfitsAll$axtit=="Met incidence"] = "Met detection"
	r.valfitsAll$y[r.valfitsAll$y==-1]<-NA
   
   	drawAtotal2<- function(r.fits1a){
		                                              
		axtita = c("0.5-th time to death quantile","Med time to death w. mets")
		
		r.fits1a = r.fits1a[r.fits1a$axtit %in% axtita,]
		r.fits1a$axtit[r.fits1a$axtit=="0.5-th time to death quantile"]="All patients"
		r.fits1a$axtit[r.fits1a$axtit=="Med time to death w. mets"]="Patients with mets"
		
		
		Atotal2 <- ggplot(r.fits1a, aes(x=diameter, y = y)) +
		theme_classic()+ 
		geom_point(shape=20, size = 0.8)+
		geom_line(data = r.fits1a, aes(x=diameter, y = pred, color= axtit))+
		scale_color_manual( values=c(blue,red)) + 
		scale_y_continuous(name="Median time to death (years)")+#, limits = c(0,7)
		scale_x_continuous(name="Diameter (cm)", breaks=c(0,2,4,6,8, 10), labels=c("0","2","4","6","8","10"))+
		facet_grid(cancer ~ axtit, scales = "fixed", space = "fixed")+
		 theme(axis.title=element_text(size=9), 
		 axis.text.y = element_text(size=8), 
		 plot.margin = 	margin(3, 1, 0, 0),
		axis.text.x = element_text(size=8), 
		panel.background = element_rect(fill="white", color = "gray60"),  
		strip.text.y = element_text(size=8), 
		strip.text.x = element_text(size=8 , lineheight=0.8), 
		strip.background = element_rect(colour="gray60", fill="white",   size=0.5, linetype="solid"), 
		plot.title = element_blank(), 
		legend.position="none")
		Atotal2
	}
	
	canc.sel=c("Gastric","Colon" ,"Renal")   
	r.fits1a.main = r.valfitsAll[r.valfitsAll$cancer %in% canc.sel,]

	Atotal2.main<-drawAtotal2(r.fits1a.main)
	
	r.fits1a.ext = r.valfitsAll[!r.valfitsAll$cancer %in% canc.sel,]
	Atotal2.ext<-drawAtotal2(r.fits1a.ext)
	
	
	
	drawAtotal3<- function(r.fits1){

		axtitb = "Met detection"
		r.fits1b = r.fits1[r.fits1$axtit==axtitb,]
	
		axtitcs =c( "Met probability", "Met probability up","Met probability down")
		axtitc = "Met probability"
		r.fits1cpoint = r.fits1[r.fits1$axtit == axtitc,]
		r.fits1c = r.fits1[r.fits1$axtit %in% axtitcs,]
		r.fits1c.c = rbind(r.fits1c, r.fits1b)
		
		
		r.fits1c.c$datatype = rep("Met detection", nrow(r.fits1c.c))	
		r.fits1c.c$datatype[r.fits1c.c$axtit =="Met probability"]= "Cancer death"
		r.fits1c.c$datatype[r.fits1c.c$axtit =="Met probability up"]= "Cancer death"
		r.fits1c.c$datatype[r.fits1c.c$axtit =="Met probability down"]= "Cancer death"
		r.fits1c.c$datatype = factor( r.fits1c.c$datatype, levels = c("Met detection","Cancer death") )
		r.fits1cpoint$datatype = factor( rep("Cancer death", nrow(r.fits1cpoint)), levels = c("Met detection","Cancer death"))
		r.fits1b$datatype = factor( r.fits1b$axtit, levels = c("Met detection","Cancer death"))
		Atotal3 <- ggplot(r.fits1c.c, aes( x=diameter, y =  y, color=axtit, shape = axtit  )) +
		theme_classic()+
		scale_color_manual( values= c("black","black","gray20","gray20" ))+
		scale_shape_manual( values = c(20,20,6,2))+
		geom_point( size = 0.8 )+
		facet_grid(cancer ~ datatype, scales = "fixed", space = "fixed" )+
		geom_line(data = r.fits1cpoint, aes(x=diameter, y = pred), color = red )+
		geom_line(data = r.fits1b, aes(x=diameter, y = pred), color = blue )+
		scale_y_continuous(name= "Probability", breaks= c( 0, 0.2, 0.4, 0.6, 0.8, 1.0 ), 
		labels=c( "0.0", "0.2", "0.4", "0.6", "0.8", "1.0" ) ) +
		scale_x_continuous(name="Diameter (cm)", breaks=c( 0, 2, 4, 6, 8, 10 ), labels=c( "0", "2", "4", "6", "8", "10" )) +
		 theme( axis.title=element_text(size=9),
		 axis.text = element_text(size=8),
		 plot.margin = margin(3, 1, 0 ,0),
		panel.background = element_rect(fill="white", color = "gray60"), 
		strip.text = element_text(size=8), 
		strip.background = element_rect(colour="gray60", fill="white",   size=0.5, linetype="solid"), 
		#plot.title = element_text(lineheight = 0.05, hjust = -0.25, vjust = -65, face="bold",size=11), 
		plot.title=element_blank(),
		legend.position="none")
		Atotal3
	}
	
	Atotal3.main = drawAtotal3(r.fits1a.main)
	Atotal3.ext = drawAtotal3(r.fits1a.ext)
			
	km.data$PSurvival = km.data$Survival*100
	km.data$Cohort = factor(km.data$Cohort, levels = c("Autopsy", "Adjuvant"))
	appr$Model1 = factor(appr$Model, levels = c("Prediction","Haeno et al.","Data"))
	Atotal4 = ggplot(km.data, aes(x = Years, y = PSurvival))+geom_line()+
	geom_ribbon(data= km.data,aes(ymin=lci*100,ymax=uci*100),alpha=0.17)+
	geom_line(data = appr, aes(x=years, y = survival*100, color = Model1))+
	scale_color_manual(values=c(red,orange,"black"))+
	#scale_fill_manual(name="Parameter",values=c("darkgreen","yellowgreen"),labels=expression(t[0],t[1]))+
	facet_grid( .~Cohort, scales = "fixed", space = "fixed")+
	theme_classic()+ggtitle("e")+
	guides(color=guide_legend(title=NULL,ncol=1))+
	scale_y_continuous(name="% Survival")+
	scale_x_continuous(name="Years", limits = c(0,4.5))+
	theme(plot.margin = margin(10, 45, 1, 0), 
	axis.title=element_text(size=9),
	axis.text.y = element_text(size=8) ,axis.text.x = element_text(size=8),
	plot.title = element_blank(),
	panel.background = element_rect(fill="white", color = "gray60"), 
	strip.text = element_text(size=8), 
	legend.position ="none",
	strip.background = element_rect(colour="gray60", fill="white", size=0.5, linetype="solid"))
	
	 
 
	 pdf("Figure2.pdf", height=8, width=2.9)
	  multiplot( Atotal3.main, Atotal2.main ,  Atotal4, layout=matrix( 
	 c(  rep( 1, 22),
	 	 rep( NA, 4), 
	     rep( 2, 21),
	     rep( NA, 4),
	     rep(3, 10)     
	   ), byrow=T, ncol=1))
	 dev.off()
	 
	 	 pdf("ExtendedData_Figure1.pdf", height=7.5, width=6.1)
	  multiplot( Atotal3.ext, Atotal2.ext ,layout=matrix( 
	 c(  rep( 1, 20),NA, rep( 2, 20)), byrow=T, nrow=1))
	 dev.off()

}

## Draw Figure 3 with fitted parameters
Figure3 <- function(parameters, preds, r.valfitsAll, par.h){
	
	parameters$Cancer = factor( replaceCancers(( parameters$Cancer)),levels=replaceCancers(levels( parameters$Cancer)))
	parameters = parameters[!parameters$Cancer=="Ovarian",]
	# order by the values
	tot = parameters[ parameters$Parameter=="t0", c("Time","Cancer")]
	# aggregate(parameters[parameters$Parameter%in%c("t0","h","a"),]$Time, by=list( Cancer= parameters[parameters$Parameter%in%c("t0","h","a"),]$Cancer), FUN=sum)
	
	Cancers = as.character( tot$Cancer[order(tot$Time)])
	
	parameters$Cancer = factor( parameters$Cancer, levels = Cancers)
	
	parameters1a = parameters[parameters$Parameter%in%c("t0","t1"),]
	parameters1b = parameters[parameters$Parameter%in%c("a"),] 
	
	
	
	##### Time interval parameters -- panel a
  p1a <- ggplot(parameters1a,aes(x=Cancer,y=Time,fill=Parameter))+
  geom_bar(stat="identity",position="dodge")+
  theme_classic()+
  coord_flip()+
   	labs(y="")+
  scale_y_continuous(name="")+ 
  scale_fill_manual(name="Parameter" ,values=c("yellowgreen","darkgreen"),labels= expression(delta[0],delta[1]) )+
  guides(fill=guide_legend(title=NULL, ncol =1))+
  theme(
  axis.title.x=element_text(size=8),
  axis.title.y=element_blank(),
  axis.text.y = element_text(size=8),
  axis.text.x = element_text(size=8),
  panel.background = element_rect(fill="white", color = "gray60"), 
  panel.grid.major.x = element_line( size=.1, color="gray" ) ,
  panel.grid.major.y = element_line( size=.1, color="gray" ) , 
  legend.position = "bottom", 
  legend.title =element_text(size = 7), 
  legend.text=element_text(size=7), 
  legend.background = element_rect(fill=alpha('white', 0.1)),
  legend.key.size = unit(0.25, "cm") ,
  legend.margin=margin(-1, 0,1,-1),
  plot.margin = 	margin(1, 1, 1, 1) )


	#### Analyze what the values in med time to death w mets come from -parameter h
	r.valfitsAll$Cancer = factor( replaceCancers(( r.valfitsAll$cancer)),levels=replaceCancers(levels( r.valfitsAll$cancer)))
	par.h$Cancer = factor( replaceCancers(( par.h$Cancer)),levels=replaceCancers(levels( par.h$Cancer)))
		
	r.w.mets = r.valfitsAll[r.valfitsAll$axtit=="Med time to death w. mets"  ,]
	avg.w.mets.pred = aggregate(r.w.mets$pred, by=list( Cancer= r.w.mets$Cancer), FUN=mean, na.rm = T)
	avg.w.mets.pred$type= rep("Prediction",nrow(avg.w.mets.pred))
	
	summ.h  = data.frame( rbind(avg.w.mets.pred, data.frame(Cancer = par.h[, "Cancer"], x = par.h[,"Time"], type = rep("Prolongation h", nrow(par.h)))) )
	summ.h$type = factor(summ.h$type)
	summ.h = summ.h[summ.h$Cancer!="Ovarian",]
	summ.h$Cancer = factor(summ.h$Cancer, levels = Cancers)
	Atotalh<- ggplot(summ.h, aes(x = Cancer, y = x, fill=type)) +
  	geom_bar(stat = "identity",position="dodge")+
  	scale_fill_manual( values=c( "red", "steelblue")) + 
 	theme_classic()+
 	coord_flip()+
 	labs(y="")+
 	guides(fill=guide_legend(title=NULL, ncol = 1))+
 	 theme(axis.title.y=element_blank(),
 	 	axis.text.y = element_blank(),
		axis.text.x = element_text(size=8), 
		panel.background = element_rect(fill="white", color = "gray60"),
		  panel.grid.major.x = element_line( size=.1, color="gray" ) ,
  		panel.grid.major.y = element_line( size=.1, color="gray" ) , 
        panel.grid.major = element_blank(), 
        plot.margin = 	margin(1, 1, 1, 2),
		strip.text.y = element_text(angle=0,hjust=0, size=8), 
		strip.background = element_rect(colour="gray60", fill="white",   size=0.5, linetype="solid"), 
		axis.title.x = element_text(size=8),
		legend.position="bottom",
		legend.background = element_rect(alpha('white', 0.0)),
		legend.text=element_text(size=7), 
		legend.key.size = unit(0.25, "cm"),
		legend.margin=margin(-1, 0,1,-1))


### Generate data points from the b distribution for each cancer 

setaside = c("Ovarian")
parameters1 = parameters[!parameters$Cancer%in%setaside, ]
stds = parameters1[parameters1$Parameter == "std", "Time"]
medians = parameters1[parameters1$Parameter == "median", "Time"]

b_dist = NULL
locs = parameters[parameters$Parameter == "location", ]
scals=parameters[parameters$Parameter == "scale", ]

for (cancer in Cancers){
	c1 = locs[locs$Cancer ==cancer, "Time"]
	c2 = scals[scals$Cancer == cancer, "Time"]
	std = parameters[(parameters$Parameter == "std"&parameters$Cancer == cancer),"Time"]
	me = parameters[(parameters$Parameter == "median"&parameters$Cancer == cancer),"Time"]
	
	if (std>5)
		var = "High"
	else
		var = "Low"
	no = 100000
	b_dist = rbind(b_dist, data.frame( LogNorm = rlnorm(no, meanlog = c1, sdlog = c2), Cancer = rep(cancer, no), Var = rep(var,no ), Std = rep(std, no), Med = rep(me, no)))
}


b_dist1 = b_dist[!b_dist$Cancer %in% setaside,]
Cancers2 = as.character( unique(b_dist1$Cancer) )
Stds2 = unique(b_dist1$Std)
Meds2 = unique(b_dist1$Med)

Cancers2 = Cancers2[order(Meds2)]

b_dist1$Cancer = factor(b_dist1$Cancer, levels = Cancers2)


p2a = ggplot(b_dist1, aes(x=Cancer, y=LogNorm)) +
 #geom_violin(trim=FALSE, aes(fill=Var), color="black", alpha=0.5)+scale_colour_gradient(limits=c(3, 4), low="red")
 #scale_fill_manual(values=c("purple","orange"))+
 geom_violin(trim=FALSE, color="black", alpha=0.6)+
 scale_fill_gradient(limits=c(min(b_dist1$Std), max(b_dist1$Std)), low="yellow", high="purple4")+
 labs(y="Bottleneck b")+
 geom_boxplot(outlier.size = 0.1, fill=NA,alpha=0.1,width=0.8)+
 theme_classic()+
 theme(
 axis.text.x = element_text( size=8),
 axis.title.y=element_blank(),
 axis.title.x = element_text( size=8),
 axis.text.y = element_text(size=8),
 plot.title = element_blank(),
  panel.background = element_rect(fill="white", color = "gray60") ,
 plot.margin = 	margin(1, 1, 27, 1) ,
  legend.position = "bottom",
 #legend.justification="right",
 legend.text=element_text(size=7), 
 legend.title=element_text(size=7),
 legend.key.size = unit(0.25, "cm"),
 legend.margin=margin(-1, -1,-1,-1)
 )+	coord_flip(ylim = c(0,35))#+ coord_cartesian(ylim = c(0,35))

par2a = parameters1b
par2a$Cancer = factor(par2a$Cancer, levels = Cancers)

p1b.a <- ggplot(par2a,aes(x=Cancer,y=Time))+#,fill="yellow"
  geom_bar(stat="identity",position="dodge", width = 0.7)+
  theme_classic()+
  coord_flip()+
 scale_y_continuous(name="",breaks = c(0, 0.5, 1), labels=c("0.0", "0.5", "1.0"),limits=c(0, 1.1))	+
 # scale_fill_gradient(limits=c(min(b_dist1$Std), max(b_dist1$Std)), low="yellow", high="purple4")+
  guides(fill=guide_legend(direction="horizontal"))+
  theme(axis.text.x = element_text(size=8), 
 axis.title.x = element_text( size=8),
  axis.title.y=element_blank(),
  axis.text.y = element_blank(),
 legend.position="bottom",
 		legend.key.size = unit(0.25, "cm"),
  panel.background = element_rect(fill="white", color = "gray60"), 
  panel.grid.major.x = element_line( size=.1, color="gray" ) ,
  panel.grid.major.y = element_line( size=.1, color="gray" ) ,
plot.margin = 	margin(1, 1, 27, 2) )



pdf("Figure3.pdf", height=2.5, width=5.7)
#multiplot(p1a, p1b.a, p1b.b, p1b.c,  p2a, p2b , p2c, pred.f.a, pred.f.b, pred.c,  pred.e, pred.d,
multiplot(p1a, Atotalh, p1b.a, p2a, layout=matrix( 
 c( rep(1, 52),  rep(2, 24), rep(3,18), NA,NA,NA, NA,rep(4, 47) ), byrow=T, nrow = 1) )
dev.off()

}


################################################################################################
################################################################################################
#### Predictions. For Figure 4 equivalent in the paper, use version = "main", 
#### for Extended data figure use version = "supp"
################################################################################################
################################################################################################

Figure4<-function(preds, ver="main"){

hth=3.4
#selcanc = c("Head & neck" ,"Colon","Breast lob")
selcanc = c("Renal","Pancreatic", "Breast lob")
#selcanc = c("Colon","Gastric", "Renal")

tit="Figure4.pdf"

if (ver=="supp"){
	hth=9.2
	#selcanc = c( "Endometrial","Bladder", "Rectal",'Lung',"Breast","Colon muc","Renal","Esophageal","Gastric", "Pancreatic")
	#selcanc = c( "Head & neck" ,"Colon","Endometrial","Bladder", "Rectal",'Lung',"Breast","Colon muc","Esophageal","Gastric")
	selcanc = c( "Head & neck" ,"Pancreatic","Endometrial","Bladder", "Rectal",'Lung',"Breast","Colon muc","Esophageal","Breast lob")
	tit="ExtendedData_Figure3.pdf"
}
	
selcanc=sort(selcanc)
preds1 = data.frame(preds)
preds1$Cancer = factor(preds1$Cancer)
preds1$Cancer = factor( replaceCancers(( preds1$Cancer)),levels=replaceCancers(levels(  preds1$Cancer)))

#### Printouts before cancers get selected for Figure 4.
prints = preds1[preds1$pred%in%c("% met prob increase due to surg delay","% met prob decrease due to immunotherapy","% met prob decrease due to chemotherapy"),]

doPrint<- function( sel, tit, add="Marginal" ){
	prints = 	prints[prints$pred==sel,]
	prints.marg = prints[prints$additional ==add,]
	prints.indv = prints[prints$additional !=add,]
	print(tit)
	cat("Maximum", add, "\n")
	print( prints.marg[which.max(prints.marg$y),])

	print("Maximum individual")
	print( prints.indv[which.max(prints.indv$y),])


}

doPrint( "% met prob increase due to surg delay","Surgery delay")
doPrint( "% met prob decrease due to chemotherapy","Chemotherapy boost")
doPrint("% met prob decrease due to immunotherapy", "Vaccine")


preds1 = preds1[preds1$Cancer%in%selcanc,]
preds1$Cancer = factor(preds1$Cancer, levels = selcanc)


preds11 = preds1[preds1$pred%in%c("% met prob increase due to surg delay","% met prob decrease due to immunotherapy","% met prob decrease due to chemotherapy"),]
colnames(preds11) = c("Diameter","y","pred","Mb","Change","Cancer")
#### For the sake of clarity of the figure in the main text, we remove the weaker change 

if (ver=="main"){
	preds11 <- preds11[!preds11$Change%in%c("8", "10%"),]
}else{
	preds11$Change= factor(preds11$Change, levels=c("16","8","20%","10%"))
}
preds11$Mb = factor(preds11$Mb)
preds11$y[preds11$y == -Inf] <- NA
colores = rev( c("black","#B79F00"  , "#00BA38" ,"#619CFF"  ) )

axtit = "% met prob increase due to surg delay"
preds.c = preds11[preds11$pred==axtit,]

pred.c <- ggplot(preds.c, aes(x=Diameter, y =(y), color=Mb, linetype = Change)) +
	 theme_classic()+
	 geom_line()+
	 guides(color=guide_legend(ncol=3),linetype=guide_legend(ncol=2)) +
	scale_color_manual(values= colores)+#, , "#F564E3"#"#F8766D",  "#00BFC4"
	 scale_y_continuous(name="Increase of cancer death probability")+   #scale_x_discrete(name="k")+
	 scale_x_continuous(name="Diameter (cm)", breaks=c(0,2,4,6,8, 10), labels=c("0","2","4","6","8","10")) +
	 facet_grid(Cancer ~ ., scales = "fixed", space = "fixed")+
	  theme( axis.title=element_text(size=8),
	  axis.text.y = element_text(size=8),
	  plot.margin = 	margin(10, 1, 50,1),
		 axis.text.x = element_text(size=8),
		 panel.background = element_rect(fill="white", color = "gray60"), 
		 strip.text = element_text(size=8), 
		 strip.background = element_rect(colour="gray60", fill="white",   size=0.5, linetype="solid"), 
		plot.title = element_text(size=8),
		 legend.position = c(0.25,-0.3),
  legend.text=element_text(size=7), 
   legend.title=element_blank(),
   legend.background = element_rect(alpha('white', 0.0)),
  legend.key.size = unit(0.3, "cm"),
  legend.margin=margin(-1, -1,-1,-1))+ggtitle("By surgery delay")



axtit = "% met prob decrease due to immunotherapy"
preds.d = preds11[preds11$pred==axtit,]

pred.d <- ggplot(preds.d, aes(x=Diameter, y =(y), color= Mb , linetype = Change)) +
	theme_classic()+
	geom_line()+
	guides(color=guide_legend(ncol=6), linetype=guide_legend(ncol=2)) +
	scale_color_manual(values= colores)+#"#F8766D",
	scale_y_continuous(name= "Decrease of cancer death probability")+ 
	scale_x_continuous(name="Diameter (cm)", breaks=c(0,2,4,6,8, 10), labels=c("0","2","4","6","8","10")) + #scale_x_discrete(name="k")+
	facet_grid(Cancer ~ ., scales = "fixed", space = "fixed")+
	guides(color=guide_legend(ncol=4)) +
	 theme( axis.title=element_text(size=8),
	 axis.text.y = element_text(size=8),
	 plot.margin = 	margin(10, 1, 50,1),
		axis.text.x = element_text(size=8),
		panel.background = element_rect(fill="white", color = "gray60"), 
		strip.text = element_text(size=8), 
		strip.background = element_rect(colour="gray60", fill="white",   size=0.5, linetype="solid"), 
			plot.title = element_text(size=8),
		legend.position = c(0.25,-0.3),
		 legend.text=element_text(size=7), 
		 legend.background = element_rect(alpha('white', 0.0)),
  legend.title=element_blank(),
 legend.key.size = unit(0.3, "cm"),
 legend.margin=margin(-1, -1,-1,-1))+ggtitle("By bottleneck boost")

 
 axtit = "% met prob decrease due to chemotherapy"
preds.e = preds11[preds11$pred==axtit,]
if (sum(preds.e$y<0, na.rm=T)>0){
preds.e[which(preds.e$y<0) ,]$y<-NA
}
pred.e <- ggplot(preds.e, aes(x=Diameter, y =(y), color= Mb, linetype=Change )) +
	theme_classic()+
	geom_line()+
	#scale_color_manual(values=rev( c("chartreuse","cyan2","brown1","darkorchid4") ))+
	guides(color=guide_legend(ncol=6), linetype=guide_legend(ncol=2)) +
	scale_color_manual(values= colores )+#"#F8766D",
	scale_y_continuous(name="Decrease of cancer death probability",breaks=c(0,0.01,0.02,0.03), labels=c("0.00","0.01","0.02","0.03"),limits=c(0,0.035) )+ 
	scale_x_continuous(name="Diameter (cm)", breaks=c(0,2,4,6,8, 10), labels=c("0","2","4","6","8","10")) + #scale_x_discrete(name="k")+
	facet_grid(Cancer ~ ., scales = "fixed", space = "fixed")+
	guides(color=guide_legend(ncol=4)) +
	 theme( axis.title=element_text(size=8),
	 axis.text.y = element_text(size=8),
	 plot.margin = margin(10, 1, 50,1),
		axis.text.x = element_text(size=8),
		panel.background = element_rect(fill="white", color = "gray60"), 
		strip.text = element_text(size=8), 
		strip.background = element_rect(colour="gray60", fill="white",   size=0.5, linetype="solid"), 
		plot.title = element_text(size=8),
	legend.position = c(0.25,-0.3),
	legend.background = element_rect(alpha('white', 0.0)),
 legend.text=element_text(size=7), 
  legend.title=element_blank(),
 legend.key.size = unit(0.3, "cm"),
 legend.margin=margin(-1, -1,-1,-1))+ggtitle("By chemotherapy boost")
 


	pdf(tit, height=hth, width=5.5)
	multiplot( pred.c,  pred.e, pred.d,
	 layout=matrix( c( rep(1, 30), rep(NA,4), rep(2,29), rep(NA,4), rep(3, 29) ), byrow=T,nrow=1) )
	dev.off()
	
	if (ver=="supp"){
		small=preds.c[preds.c$Cancer%in%c("Colon","Colon muc", "Head & neck")&preds.c$Mb=="Marginal"&preds.c$Change==16,]
		print("Largest from small marginal impacts of surgery delay")
		print(small[which.max(small$y),])
	}
}

### This combines several summaries in one figure - not used in the paper in the end.
Figure2_Large<- function(r.fits, r.val, r.valfitsAll, km.data, appr, par.h){
	palered = "tomato"
	red="red3"
	darkred = "red4"
	blue="royalblue3"
	orange="chartreuse4"
	orange="darkorange2"

	r.valfits = rbind(r.fits, r.val)

	r.valfits$Cancer = factor( replaceCancers(( r.valfits$Cancer)),levels=replaceCancers(levels( r.valfits$Cancer)))
	r.valfits$Cancer = makeReplacement(r.valfits$Cancer)
	r.fits$Cancer= factor( replaceCancers(( r.fits$Cancer)),levels=replaceCancers(levels( r.fits$Cancer)))
	r.fits$Cancer = makeReplacement(r.fits$Cancer)
	
	####Ovarian cancer data can be fitted but we do not trust nor want to report the results- it seems like Ovarian data does not follow model assumptions!
	r.valfits = r.valfits[!r.valfits$Cancer =="Ovarian",] 
	r.fits = r.fits[!r.fits$Cancer =="Ovarian",]
		
	##tot = aggregate(r.valfits$NRMSE, by=list(Cancer= r.valfits$Cancer), FUN=function(x){sum(x, na.rm = T)/sum(!is.na(x)) })
	# order by the fits
	tot = aggregate(r.fits$NRMSE, by=list( Cancer= r.fits$Cancer), FUN=function(x){sum(x, na.rm = T)/sum(!is.na(x))})
	Cancers = as.character( tot$Cancer[order(tot$x)])
	r.valfits$Cancer = factor( r.valfits$Cancer, levels = Cancers)
	#r.valfits$Data = gsub("0.5-th time to death quantile", "Median time to death", r.valfits$Data )
	#r.valfits$Data = gsub("Met incidence", "Metastasis incidence", r.valfits$Data )
	r.valfits$Data = gsub("Min metastasis probability", "Cancer death", r.valfits$Data )

	r.valfits$Data = gsub(" time to death quantile", "", r.valfits$Data )
	r.valfits$Data = gsub(" rate", "", r.valfits$Data )
	r.valfits$Data = gsub("0.5-th", "0.50-th", r.valfits$Data )
	r.valfits$Data = gsub("0.50-th w. mets", "Median time to\n death w. mets", r.valfits$Data )


	r.valfits$Data = factor(r.valfits$Data, levels = unique(r.valfits$Data))
	
	r.valfits$Cancer2=r.valfits$Cancer
	
	r.valfits.b = r.valfits[!r.valfits$Data%in%c("Met incidence","Cancer death", "Metastasis probability"),]
	r.valfits.a = r.valfits[r.valfits$Data%in%c("Met incidence","Cancer death"),]
 
	#With faceting
#	r.valfits.sel = r.valfits[r.valfits$Cancer%in%c("Gastric","Colon","Renal"),]
	
	
	Atotal1.a<- ggplot(r.valfits.a, aes(x = Data, y = NRMSE, fill=Type)) +
  	geom_bar(stat = "identity")+scale_fill_manual( values=c(blue,palered, red, darkred)) + 
 	 facet_grid(Cancer2 ~ ., scales = "fixed", space = "fixed") + 
 	 theme_classic()+
 	 scale_y_continuous(name="Root mean squared error", breaks=c(0,0.1, 0.2), labels=c("0","0.1","0.2"), limits=c(0, 0.24))+
 	 scale_x_discrete(name="Rates")+
 	 # ggtitle("a")+
 	  guides(fill=guide_legend(title=NULL))+
 	 theme(axis.title.y=element_text(size=9),
 	 	axis.text.y = element_text(size=8),
		axis.text.x = element_text(angle=60, hjust = 1, size=8), 
		panel.background = element_rect(fill="white", color = "gray60"),
           panel.grid.major = element_blank(), 
           plot.margin = 	margin(3, 8, 35, 0),
			strip.text.y = element_blank(),
			strip.background = element_blank(),
			axis.title.x = element_text(hjust=0, vjust = 14, color="gray30", size = 8),
			#plot.title = element_text(lineheight = 0.05, hjust = -0.25, vjust = -67, face="bold",size=12),
			plot.title= element_blank(),
			 legend.position = c(1.2, -0.195)
, legend.justification="right",legend.text=element_text(size=7), legend.key.size = unit(0.25, "cm"),
panel.spacing.y=unit(0.3, "lines"), legend.margin=margin(-1, -1,0,-1))

	Atotal1.b<- ggplot(r.valfits.b, aes(x = Data, y = NRMSE, fill=Type)) +
  	geom_bar(stat = "identity")+scale_fill_manual( values=c(blue,red)) + 
 	 facet_grid(Cancer2 ~ ., scales = "fixed", space = "fixed") + 
 	 theme_classic()+
 	 scale_y_continuous(name="Root mean squared error")+#, breaks=c(0, 0.5, 1, 1.5), labels=c("0","0.5", "1.0", "1.5"), limits=c(0, 1.52))
 	 scale_x_discrete(name="Time to death quantiles")+
 	 # ggtitle("a")+
 	  guides(fill=guide_legend(title=NULL))+
 	 theme(axis.title.y=element_blank(),
 	 	axis.text.y = element_text(size=8),
		axis.text.x = element_text(angle=60, hjust = 1, size=8), 
		panel.background = element_rect(fill="white", color = "gray60"),
           panel.grid.major = element_blank(), 
           plot.margin = 	margin(3, 0, 26, 0),
			strip.text.y = element_text(angle=0,hjust=0, size=8), 
			strip.background = element_rect(colour="gray60", fill="white",   size=0.5, linetype="solid"), 
			axis.title.x = element_text(hjust=0.0, vjust = 13, color="gray30", size = 8),
			#plot.title = element_text(lineheight = 0.05, hjust = -0.22, vjust = -68, face="bold",size=11), 
			plot.title = element_blank(), 
			legend.position = c(0.4, -0.18)
, legend.justification="right",legend.text=element_text(size=7), legend.key.size = unit(0.25, "cm"),
panel.spacing.y=unit(0.3, "lines"), legend.margin=margin(-1, -1,0,-1))



     r.valfits1 = r.valfitsAll[r.valfitsAll$cancer %in% c("Gastric AdenoCa","Colon AdenoCa" ,"Renal Cell Ca"),]                       
    r.valfits1$cancer = factor( replaceCancers(( r.valfits1$cancer)),levels=replaceCancers(levels( r.valfits1$cancer)))
	 r.valfits1$cancer = makeReplacement(r.valfits1$cancer)	
      lev = Cancers[Cancers%in%r.valfits1$cancer] 
      r.valfits1$cancer = factor( r.valfits1$cancer, levels = lev)                                	
	
	axtita = c("0.5-th time to death quantile","Med time to death w. mets")
	r.fits1a = r.valfits1[r.valfits1$axtit %in% axtita,]
	r.fits1a$axtit[r.fits1a$axtit=="0.5-th time to death quantile"]="Overall"
	r.fits1a$axtit[r.fits1a$axtit=="Med time to death w. mets"]="With mets"
	#r.fits1apoint = r.fits1a[r.fits1a$axtit=="With mets",]
	#r.fits1apoint$h[r.fits1apoint$Cancer == "Gastric"] = par.h$Time[par.h$Cancer =="Gastric"]
	#r.fits1apoint$h[r.fits1apoint$Cancer == "Colon"] = par.h$Time[par.h$Cancer =="Colon"]
	#r.fits1apoint$h[r.fits1apoint$Cancer == "Renal"] = par.h$Time[par.h$Cancer =="Renal"]
	
	
	Atotal2 <- ggplot(r.fits1a, aes(x=diameter, y = y)) +
	theme_classic()+ 
	#ggtitle("b")+
	geom_point(shape=1, size = 0.8)+
	geom_line(data = r.fits1a, aes(x=diameter, y = pred, color= axtit))+
	#geom_line(data = r.fits1apoint, aes(x=diameter, y = h), color = "darkseagreen2",linetype="dotted" )+
	scale_color_manual( values=c(blue,red)) + 
	scale_y_continuous(name="Median time to death (years)")+#, limits = c(0,7)
	scale_x_continuous(name="Diameter", breaks=c(0,2,4,6,8, 10), labels=c("0","2","4","6","8","10"))+
	facet_grid(cancer ~ axtit, scales = "fixed", space = "fixed")+
	 theme(axis.title=element_text(size=9), 
	 axis.text.y = element_text(size=8), 
	 plot.margin = 	margin(3, 0, 0, 0),
	axis.text.x = element_text(size=8), 
	panel.background = element_rect(fill="white", color = "gray60"),  
	strip.text.y = element_text(size=8), 
	strip.text.x = element_text(size=8 , lineheight=0.8), 
	strip.background = element_rect(colour="gray60", fill="white",   size=0.5, linetype="solid"), 
#	plot.title = element_text(lineheight = 0.05, hjust = -0.25, vjust = -73, face="bold",size=11), 
	plot.title = element_blank(), 
	legend.position="none")
		#,panel.grid.minor.y = element_line(color = "gray80"))
	
	
#  plot.margin = 	margin(-13, 0, -13, 0),
	
	
	#plot.margin = 	margin(7, 1, 90, 8),
	axtitb = "Met incidence"
	r.fits1b = r.valfits1[r.valfits1$axtit==axtitb,]

	axtitcs =c( "Met probability", "Met probability up","Met probability down")
	axtitc = "Met probability"
	r.fits1cpoint = r.valfits1[r.valfits1$axtit == axtitc,]
	r.fits1c = r.valfits1[r.valfits1$axtit %in% axtitcs,]
	r.fits1c.c = rbind(r.fits1c, r.fits1b)
	r.fits1c.c$datatype = rep("Met detection", nrow(r.fits1c.c))
	r.fits1c.c$datatype[r.fits1c.c$axtit =="Met probability"]= "Cancer death"
	r.fits1c.c$datatype[r.fits1c.c$axtit =="Met probability up"]= "Cancer death"
	r.fits1c.c$datatype[r.fits1c.c$axtit =="Met probability down"]= "Cancer death"
	r.fits1c.c$datatype = factor( r.fits1c.c$datatype, levels = c("Met incidence","Cancer death") )
	r.fits1cpoint$datatype = factor( rep("Cancer death", nrow(r.fits1cpoint)), levels = c("Met incidence","Cancer death"))
	r.fits1b$datatype = factor( r.fits1b$axtit, levels = c("Met incidence","Cancer death"))
	Atotal3 <- ggplot(r.fits1c.c, aes( x=diameter, y =  y, color=axtit, shape = axtit  )) +
	theme_classic()+
	#ggtitle("b")+
	scale_color_manual( values= c("black","black","gray20","gray20" ))+
	scale_shape_manual( values = c(1,20,1,1))+
	geom_point( size = 0.8 )+
	facet_grid(cancer ~ datatype, scales = "fixed", space = "fixed" )+
	geom_line(data = r.fits1cpoint, aes(x=diameter, y = pred), color = red )+
	geom_line(data = r.fits1b, aes(x=diameter, y = pred), color = blue )+
	scale_y_continuous(name= "Rate", breaks= c( 0, 0.2, 0.4, 0.6, 0.8, 1.0 ), 
	labels=c( "0.0", "0.2", "0.4", "0.6", "0.8", "1.0" ) ) +
	## scale_y_continuous(name="Metastasis probability") +
	scale_x_continuous(name="Diameter", breaks=c( 0, 2, 4, 6, 8, 10 ), labels=c( "0", "2", "4", "6", "8", "10" )) +
	#facet_grid(cancer ~ ., scales = "fixed", space = "fixed")+
	 theme( axis.title=element_text(size=9),
	 axis.text = element_text(size=8),
	 plot.margin = margin(3, 0, 0 ,0),
	panel.background = element_rect(fill="white", color = "gray60"), 
	strip.text = element_text(size=8), 
	strip.background = element_rect(colour="gray60", fill="white",   size=0.5, linetype="solid"), 
	#plot.title = element_text(lineheight = 0.05, hjust = -0.25, vjust = -65, face="bold",size=11), 
	plot.title=element_blank(),
	legend.position="none")
#	plot.margin = 	margin(1, 3, 0,-10)

				
	#	 
	#Atotal5 = ggsurvplot(km.adjuvant,palette = c("#E7B800", "#2E9FDF"),ggtheme = theme_bw())
	km.data$PSurvival = km.data$Survival*100
	appr$Model1 = factor(appr$Model, levels = c("Prediction","Haeno et al.","Data"))
	Atotal4 = ggplot(km.data, aes(x = Years, y = PSurvival))+geom_line()+
	geom_ribbon(data= km.data,aes(ymin=lci*100,ymax=uci*100),alpha=0.17)+
	geom_line(data = appr, aes(x=years, y = survival*100, color = Model1))+
	scale_color_manual(values=c(red,orange,"black"))+
	#scale_fill_manual(name="Parameter",values=c("darkgreen","yellowgreen"),labels=expression(t[0],t[1]))+
	facet_grid( .~Cohort, scales = "fixed", space = "fixed")+
	theme_classic()+ggtitle("e")+
	guides(color=guide_legend(title=NULL,ncol=1))+
	scale_y_continuous(name="% Survival")+
	scale_x_continuous(name="Years", limits = c(0,4.5))+
	theme(plot.margin = margin(3, 10,  1,0), 
	axis.title=element_text(size=9),
	axis.text.y = element_text(size=8) ,axis.text.x = element_text(size=8),
	plot.title = element_blank(),
	panel.background = element_rect(fill="white", color = "gray60"), 
	legend.position ="none",
	# c(1.2, -0), legend.direction="horizontal",legend.box = "horizontal",
	#legend.text=element_text(size=7), 
	#legend.key.size = unit(0.3, "cm"), strip.text = element_text( size=8), ,
	# legend.margin=margin(-1, 0,1,-1),
	#legend.background = element_rect(fill=alpha('white', 0.0)),
	strip.background = element_rect(colour="gray60", fill="white", size=0.5, linetype="solid"))
	#margin(-3, 20, -13, 11)
	

 
 
 pdf("Figure2_large.pdf", height=6.9, width=6)
  multiplot(Atotal1.a, Atotal1.b, Atotal3, Atotal2 ,  Atotal4, layout=matrix( 
 c(  rep( c( rep(1,4), rep(2,12), NA, rep(3,12) ), 19),
 	 rep( c( rep(1,4), rep(2,12), rep(NA, 13)  ), 1), 
     rep( c( rep(1,4), rep(2,12), NA, rep(4,12) ), 19),
     rep( c( rep(1,4), rep(2,12), rep(NA, 13)  ), 1),
     rep( c( rep(1,4), rep(2,12), NA, rep(5,12) ), 9)     
   ), byrow=T, ncol=29))
 dev.off()
}

