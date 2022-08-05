source("functions.R")
load("./1_Datos_MF/data.RData")
load("./1_Datos_MF/coord.RData")

MoranIndex<-function(i){

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
      integra<-(xmgi-xmean$data)*(xmgj-xmean$data)
      integral<-int.simpson(fdata(integra))
      numerador<-weig.mat[i,j]*integral
      s=s+numerador
      ss[i]<-s    
    }
    
    integra1<-(xmgi-xmean$data)^2
    integral2<-int.simpson(fdata(integra1))
    s2<-s2+integral2
  }
  IMoran<-(length(a)*sum(ss))/(sum(sum(weig.mat))*s2)
  return(IMoran)
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

### paralelizacion  de moran


MoranIndex(2)

cl1<- parallel::makeCluster(detectCores())
clusterEvalQ(cl1,source("functions.R"))
i=seq(1,length(results))
clusterExport(cl = cl1, varlist = c("results","adya2","i","data"))

system.time({
  resultsm <- parallel::parLapply(cl1,i,MoranIndex)
})

moranind=do.call(rbind,resultsm)
stopCluster(cl1)
moranind=as.data.frame(moranind)

moranind$gr=c(2:33)

ggplot(moranind,aes(x=gr,y=V1))+
  geom_point()+
  geom_line()+
  theme_bw()


### moran simulaciones------------------------------- corregir tendencia

set.seed(123)

c1=c(runif(30,-2,0) ,runif(30,-2,0), runif(30,0,2), runif(30,0,2))
c2=c(runif(30,-2,0) ,runif(30,0,2), runif(30,-2,0), runif(30,0,2))

coord=data.frame(cbind(c1,c2))
names(coord)=c("x","y")

Eu.d <-as.matrix(dist(coord,method="euclidian"))

# adya2<-ifelse(Eu.d<=2,1,0) 
# diag(adya2)<-0

adya2=1/Eu.d
diag(adya2)=0



ggplot(coord,aes(c1,c2))+
  geom_point(size=2)+
  geom_line(aes(x=0),linetype="dashed")+
  geom_line(aes(y=0),linetype="dashed")+
  theme_bw()


cov.model="gaussian"
Kappa=NULL
multivgm=FALSE
nbasis=65
cov.pars=c(2,4)
media=c(5,10)

set.seed(999)
data=Datsim(Eu.d,cov.model = cov.model,cov.pars = cov.pars,Kappa=Kappa,
         media = media) %>% data.frame()

cl <- parallel::makeCluster(detectCores())# Run parallel computation
clusterEvalQ(cl, source("functions.R"))
ncl=seq(2,33)
clusterExport(cl = cl, varlist = c("coord","cov.model",
                                   "Kappa","multivgm","ncl","data"))#, envir = .GlobalEnv)


system.time({
  results <- parallel::parLapply(cl,ncl,testm)
})


parallel::stopCluster(cl)


cl1<- parallel::makeCluster(detectCores())
clusterEvalQ(cl1,source("functions.R"))
i=seq(1,length(results))
clusterExport(cl = cl1, varlist = c("results","adya2","i","data"))

system.time({
  resultsm <- parallel::parLapply(cl1,i,MoranIndex)
})

moranind=do.call(rbind,resultsm)
stopCluster(cl1)
moranind=as.data.frame(moranind)

moranind$gr=c(2:33)

ggplot(moranind,aes(x=gr,y=V1))+
  geom_point()+
  geom_line()+
  theme_bw()



