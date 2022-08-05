##############################################################################
# This function estimates a linear model of coregionalziation to coefficients
# of Fourier basis functions (used for smoothing observed data).
# All single (direct) variograms and  cross-variograms of coefficients 
# of Fourier basis are modeled as a linear combination of nugget and   
# exponential models.     
##############################################################################

fitlmc<-function(coord,datafd,k)
{
   coefi<-matrix(datafd$coefs,nrow=k) 
   

   coef<-t(coefi) 
   coef<-cbind(coord,coef)
   coef<-as.data.frame(coef)
   n2<-paste("var",1:k, sep="")
   names(coef)<-c("x","y",n2)
   coordinates(coef)= ~x+y
   g<-NULL

   for (i in 1:k)
     { 
     g <- gstat(g,formula= as.formula(paste(n2[i],"~1")), data=coef)
     }
   v <- variogram(g)
   sigma2.0<-quantile(v$gamma,0.75) 
   phi.0<-quantile(v$dist,.75)
   #
   g <- gstat(g, model=vgm(psill = sigma2.0,model = "Mat",#cambiar modelos dependiendo de aplderiv.R
                           range = phi.0,nugget = 0,kappa = 0.2), fill.all=TRUE)
   g.fit <- fit.lmc(v,g,fit.lmc=TRUE,correct.diagonal=1.01)
   return(g.fit)
}
