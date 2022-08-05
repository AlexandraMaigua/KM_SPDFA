#Verificar teoría multivgm lmc
#en cov.spatial solo funciona con modelos comunes del mismo y de vgm

multiv<-function(Eu.d,coord,data,cov.model="Gau",max.dist=3,Kappa=NULL){
  k=dim(data)[1]
  varis=matrix(0,nc=120,nr=120)
  for (i in 1:k) {
    dat=cbind(coord,value=data[i,])
    coordinates(dat) = ~x+y
    emp.trace.vari=variogram(value~1, dat,dX=0,cutoff=max.dist)
    
    sigma2.0=quantile(emp.trace.vari$gamma,0.75)
    phi.0=quantile(emp.trace.vari$dist,0.2)
    nt=mean(emp.trace.vari$gamma)/4
    if(is.null(Kappa)==TRUE){
      trace.vari = fit.variogram(emp.trace.vari, vgm(sigma2.0, cov.model, phi.0,nt))
      
    }
    else{
      trace.vari = fit.variogram(emp.trace.vari, vgm(psill = sigma2.0,model =  cov.model, 
                                                     range =  phi.0, nugget =  nt,kappa = Kappa),
                                 fit.kappa = TRUE)
    }
    
    sigma2=trace.vari$psill[2]
    nugget=trace.vari$psill[1]
    tra.vari.mat <- sigma2+nugget - cov.spatial(Eu.d,cov.model=tolower(trace.vari$model[2]),
                                                cov.pars=c(sigma2,trace.vari$range[2]),
                                                kappa=trace.vari$kappa[2])
    
    
    
    varis=varis+tra.vari.mat
  }
  
  return(varis)
  
  
}
