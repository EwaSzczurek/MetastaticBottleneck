summarizeGrouping<-function( coding,inc ){
	histo= unique(inc[,coding])
	counts = sapply(histo, function(x){sum(inc[,coding]==x)})
	names(counts)=histo
	counts = sort(counts, decreasing = T)
	write.table(counts, file = paste(coding, ".xls",sep=""),sep="\t",quote=F)
	counts
}

checkCodings<- function(inc){
	type = "Cancers"
	cat("Analyzing",type, "data ..\n")
	load(paste(type,"Sel.Rdata",sep=""))
	cat( nrow(inc), "of all patients in the database\n")
	
	gr = c("Grouped.Coding.1", "Grouped.Coding.2",  "Grouped.Coding.3")
		
	counts1 = summarizeGrouping( gr[1],inc)
	counts2 = summarizeGrouping( gr[2],inc)
	counts3 = summarizeGrouping( gr[3],inc)
}
