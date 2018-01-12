source("A_Handy.R")

bs = c(17, 19, 21)

getjumppis<-function(b,is){
	jumppi(is,b)
}

getsis <- function( b, is){
	s1 = s1SS(b)
	sapply(is, si, b=b, s1=s1)
}

  cols=c("black","red","green")
  

 
ressi = sapply(bs, getsis, 1:40)
 
 pdf("Figure1e.pdf", height = 1.45, width=1.6)
par(ps=8,  mar=c(2,2,0.35,0.2), mgp=c(1.1,0.3,0))
matplot(ressi,type="l", lty=1, pch = 20, ylab="P(survival from i cells)", xlab="Number of cells i" , cex=0.4, las=1, tck=-0.04 )#,log="y")
#abline(h=0.5)
abline(v  = bs, lty=2, col=cols)
abline(h=0.5, col = "gray") 
dev.off()

 pdf("Figure1e_legend.pdf", height = 1.7, width=1.9)
par(ps=8, mar=c(2,2,0.2,0.2), mgp=c(1.1,0.3,0))
matplot(ressi,type="o", lty=1, pch = 20, ylab="",  cex=0.4,las=1, tck=-0.04)
#abline(h=0.5)
abline(v  = bs, lty=2, col=cols)
abline(h=0.5) 
legend("bottomright", legend=bs, pch=20, pt.cex=0.7, cex = 0.9, col =cols,bty="n",y.intersp=0.7, x.intersp=0.8)
 dev.off()

 
 

