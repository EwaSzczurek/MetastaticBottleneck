source("Metastasis_handy.R")
source("Metastasis_Survival.R")
source("Metastasis_Mets.R")
source("Metastasis_DeathCause.R")
source("Metastasis_Coding.R")

drawCounts<-function(nam,counts){
	
	pdf(paste(folder,type,"_",nam,".pdf",sep=""), height=3.5)
	par(ps=9,mar=c(6.5, 3,0.4,0.1), mgp=c(1.7,0.4,0))
	ylim <- c(0, 1.2*max(counts))
	## Plot, and store x-coordinates of bars in xx
	xx <- barplot(counts, xaxt = 'n', xlab = '', width = 0.4, ylim = ylim,
              ylab = "# patients")
	## Add text at top of bars
	text(x = xx, y = counts, label = counts, pos = 3, cex = 0.7, col = "red")
	## Add x-axis labels 
	axis(1, at=xx, labels=names(counts), tick=FALSE, las=2, line=-0.1, cex.axis=1)
	dev.off()
}

### The file CancersSel.RData is available upon request, only to users with access to SEER granted
type = "Cancers"

### Initializing summary data frames with NULL values
initialize()

cat("Analyzing",type, "data with bn size", binSize,"..\n")
load(paste(type,"Sel.RData",sep=""))
dir.create(folder, showWarnings=F)

counts = summarizeGrouping( gr,inc)
#counts= counts[counts>10000]
counts = counts[names(counts)!=""]

cat( nrow(inc), "of tumors in the database\n")
inc = inc[inc[,gr] %in%names(counts), ] 
inc$cancertype = inc$Site
inc$Site = inc[,gr]

inc <- inc[!inc$Site%in% c("thyroid Ca","cervix SCC","HCC", "SCLC"),] #have less thn 500 patients in the end.-too few datapoints to fit tumor dimeter dependence.
inc$Site[ inc$Site == "esophageal AdenoCa"] <-"Esophageal Ca"
inc$Site[ inc$Site ==  "esophageal SSC"] <-"Esophageal Ca"


cat( nrow(inc), "of patients in seleced groups\n")
drawCounts(nam = "CancerSubtypesBeforeFiltering",counts = counts)
	
cat("\nInitial filtering and preprocessing...\n")
	
########## Renaming sites to nice names
inc[inc$Site=="breast NST","Site"]<-"Breast Ca (NST)"
inc[inc$Site=="breast lobular Ca","Site"]<-"Breast Ca lobular"
inc[inc$Site=="endometrial ca","Site"]<-"Endometrial Ca"
#inc[inc$Site=="esophageal AdenoCa","Site"]<-"Esophageal AdenoCa"
#inc[inc$Site=="esophageal SSC","Site"]<-"Esophageal SSC"
inc[inc$Site=="gastric AdenoCa","Site"]<-"Gastric AdenoCa"
inc[inc$Site=="Colon-AdenoCa","Site"]<-"Colon AdenoCa"
inc[inc$Site=="Colon mucinous AdenoCa","Site"]<-"Colon AdenoCa muc"
inc[inc$Site=="Rectum Adenoca","Site"]<-"Rectum AdenoCa"
inc[inc$Site=="Renal cell Ca","Site"]<-"Renal Cell Ca"
inc[inc$Site=="urothelial cancer","Site"]<-"Urothelial Ca"



########################################################
# Limiting only to patients with a single primary
# Sequence Number-Central describes the number and sequence of all reportable malignant, in situ, benign, and borderline primary tumors, which occur over the lifetime of a patient. The sequence number may change over the lifetime of the patient. If an individual previously diagnosed with a single reportable malignant neoplasm is subsequently diagnosed with a second reportable malignant neoplasm, the sequence code for the first neoplasm changes from 00 to 01. 
#This sequence number counts all tumors that were reportable in the year they were diagnosed even if the tumors occurred before the registry existed, or before the registry participated in the SEER Program. The purpose of sequencing based on the patient’s lifetime is to truly identify the patients for survival analysis who only had one malignant primary in their lifetimes. 
# In Situ/Malignant as Federally Required based on Diagnosis Year 
# 00 One primary only in the patient’s lifetime 
# 01 First of two or more primaries 
# 02 Second of two or more primaries 
	
inc = inc[ !is.na(inc$SEQ_NUM),]
inc = inc[ inc$SEQ_NUM==0,]
cat( nrow(inc), "of patients with a single diagnosed primary tumor\n")


########################################################
##### Restricting to patients diagnosed after 1988 (before tumor size was not coded)
inc = inc[(inc$DXyear >= min.year),]
cat( nrow(inc), "of patients with diagnosed after", min.year,"\n")
	
########################################################
##### Recording met status
##### Need to analyze metastasis status in two different codings used, one before 2004 and one after
	
inc1 = inc[(inc$DXyear < mid.year),] 
inc2 = inc[inc$DXyear >= mid.year,]
	
inc1$met = ( inc1$EODextension == 85 )
inc1$metExamined = ( ( inc1$EODextension != 99 )&( !is.na(inc1$EODextension) ) )
inc1$tumorsize = inc1$EODtumorsize
	
inc2$met = ( inc2$CSMetsDx >0  )
inc2$metExamined = (inc2$CSMetsDx !=99)&( !is.na(inc2$CSMetsDx) )
inc2$tumorsize = inc2$CSTumorSize	
	
##### Putting the patients from different time bins back together
inc12 = rbind( inc1[, cols12], inc2[,cols12])
inc12 = inc12[ !is.na(inc12$tumorsize),]
cat( nrow(inc12), "of patients with non-NA tumorsize\n")

	
########################################################
inc12 = inc12[ !( inc12$tumorsize%in%c(0, 980:999) ), ]
cat( nrow(inc12), "of patients with precisely known, nonzero tumorsize in mm\n")
	
### Making size in cm and limiting to max size
inc12$tumorsize = inc12$tumorsize/10
inc12 = inc12[ inc12$tumorsize <= tumor.max, ]
cat( nrow(inc12), "of patients with tumorsize less or equal", tumor.max,"cm\n")
	
########################################################
#inc12 = inc12[inc12$DXage<500,]
inc12 = inc12[inc12$DXage>18,]
cat( nrow(inc12), "of adults\n")
	

#######################################################
### Limiting to those who had surgery
### Limiting according to #RX SUMM—SURG PRIM SITE ("SURGPRIM")
# RX SUMM—SURG PRIM SITE NAACCR Item #: 1290
# Type of surgery field
# Nur 20-80 nehmen (site specific codes)
# 00 – no surgery, 10-19 no pathology or unknown, - exclude?
### RX SUMM—SURG TYPE ("SS_SURG")  Nur 20-80 nehmen
inc12 = inc12[(inc12$SURGPRIM %in% seq(from = 20, to=80))|(inc12$SS_SURG %in% seq(from = 20, to=80)), ]
cat( nrow(inc12 ), "of patients after filtering according to the type of surgery (RX SUMM—SURG PRIM SITE or RX SUMM—SURG TYPE )\n")
	
#######################################################
### Limiting to those who did not have radiation before surgery

# RX SUMM—SURG/RAD SEQ ("RAD_SURG") NAACCR Item #: 1380
# 0 no radiation and/or surgery as above -> nehmen!
# 2 radiation before surgery -> raus
# 3 radiation after -> nehmen!
# 4 radiation before and after -> raus
# 5 intraoperative radiation -> ok
# 6 Intraoperative radiation with other radiation given before or after surgery -> ?
# 9 sequence unknown, both given -> ?
inc12  = inc12[inc12$RAD_SURG %in% c(0,3,5), ]
cat( nrow(inc12 ), "of patients after filtering according to radiation/surg order (RX SUMM—SURG/RAD SEQ)\n")

## Taking only patients with defined death status
inc12 = inc12[inc12$DTH_CLASS!=9, ]
cat( nrow(inc12), "of patients with known death status\n")


######################### Avoid censoring
#inc.surv= inc.surv[inc.surv$DXyear<=censoring.cutoff,]	
for (cn in names(censoring.cutoffs)){
	censoring.cutoff = censoring.cutoffs[cn]
	
	to.remove = (( inc12$Site == cn ) & (inc12$DXyear > censoring.cutoff))
	
	inc12 = inc12[!to.remove,]		
}

cat( nrow(inc12), "of patients with early diagnosis before or at censoring cutoffs to avoid censoring\n")



### ORDERING by tumorsize and by site -  this is assumed by any call to binByTumorSizes
inc12 = inc12[ order( inc12$tumorsize ), ]
inc12 = inc12[ order( inc12$Site ), ]

cat("Done initial filtering and data prepocessing...\n")

#### Summarize the counts after filtering
counts.1 =summarizeGrouping( "Site",inc12)
drawCounts("CancerSubtypesAfterFirstFiltering",counts.1)


########################################################
analyzeDeathCause(inc12, folder, type)	
analyzeSurvival(inc12, folder, type)
analyzeMets(inc12, folder, type)
########################################################


		 


