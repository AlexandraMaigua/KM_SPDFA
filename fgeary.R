source("functions.R")
load("./1_Datos_MF/data.RData")
load("./1_Datos_MF/coord.RData")

GearyIndex<-function(i){

  grupos <- results[[i]]
  weig.mat=adya2
  
  a<-grupos$cluster
  njs<-data.frame(table(a))
  xmean<-func.mean(fdata(t(data)))
  
  lxmean<-grupos$centers$data
  
  ss<-c()
  s2<-0
  for(i in 1:length(a)){
    gi<-a[i]
    xmgi<-lxmean[gi,] 
    s=0
    for (j in 1:length(a)) {
      gj<-a[j]
      xmgj<-lxmean[gj,]  
      integra<-(xmgi-xmgj)^2
      integral<-int.simpson(fdata(integra))
      numerador<-weig.mat[i,j]*integral
      s=s+numerador
      ss[i]<-s    
    }
    integra1<-(xmgi-xmean$data)^2
    integral2<-int.simpson(fdata(integra1))
    s2<-s2+integral2
  }
  IGeary<-(1/2)*(length(a)*sum(ss))/(sum(sum(weig.mat))*s2)
  return(IGeary)
}




## PARALELIZACION

cov.model = "Gau"
Kappa = NULL
multivgm = FALSE

testm<-function(ncl=5){
  a<-kmeans.fdas(fdataobj = t(data),ncl = ncl,Vmdist = TRUE,coord = coord,
                 cov.model = cov.model,Kappa = Kappa,multivgm = multivgm)
  return(a)
}


cl <- parallel::makeCluster(detectCores())# Run parallel computation
clusterEvalQ(cl, source("functions.R"))
ncl=seq(2,33)
clusterExport(cl = cl, varlist = c("coord","cov.model",
                                   "Kappa","multivgm","ncl","data"))#, envir = .GlobalEnv)


system.time({
  results <- parallel::parLapply(cl,ncl,testm)
})

parallel::stopCluster(cl)

### paralelizacion  de Geary


GearyIndex(2)

cl1<- parallel::makeCluster(detectCores())
clusterEvalQ(cl1,source("functions.R"))
i=seq(1,length(results))
clusterExport(cl = cl1, varlist = c("results","adya2","i","data"))

system.time({
  resultsm <- parallel::parLapply(cl1,i,GearyIndex)
})

moranind=do.call(rbind,resultsm)
stopCluster(cl1)
moranind=as.data.frame(moranind)

moranind$gr=c(2:33)

ggplot(moranind,aes(x=gr,y=V1))+
  geom_point()+
  geom_line()+
  theme_bw()

