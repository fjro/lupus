#' Function to draw screeplots of cumulative eigenvalues in principal component analysis
#' 
#' @param mydata The data.
#' @param cor T if correlation matrix should be used.
#' @param maxcomp The number of PCs to plot.
screeplot<-function(mydata, cor=F, maxcomp=10) {
  my.pc<-prcomp(mydata, scale = cor)
  k<-min(dim(mydata),maxcomp)
  x<-c(0:k)
  y<-my.pc$sdev[1:k]*my.pc$sdev[1:k]
  y<-c(0,y)
  z<-100*cumsum(y)/sum(my.pc$sdev*my.pc$sdev)
  
  plot(x,z,type="l",xlab="number of dimensions",
       cex.main=1.5, lwd=3, col="red",
       ylim=c(0,100),
       ylab="cumulative percentage of total variance",
       main="Scree plot of variances",
       xaxt="n", yaxt="n")
  
  axis(1, at=x, lwd=2)
  axis(2, at=c(0,20,40,60,80,100), lwd = 2)
  abline(a=100, b=0, lwd=2, lty="dashed", col="orange")
  text(x, z, labels=x, cex=0.8, adj=c(1.2,-.1), col="blue")
}

