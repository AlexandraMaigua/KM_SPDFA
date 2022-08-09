source("functions.R")


set.seed(123)

c1=c(runif(30,-2,-0) ,runif(30,-2,-0), runif(30,0,2), runif(30,0,2))
c2=c(runif(30,-2,-0) ,runif(30,0,2), runif(30,-2,-0), runif(30,0,2))

coord=data.frame(cbind(c1,c2))
names(coord)=c("x","y")

Eu.d <-as.matrix(dist(coord,method="euclidian"))

ggplot(coord,aes(c1,c2))+
  geom_point(size=5,color="#0099CC")+
  geom_line(aes(x=0),linetype="dashed")+
  geom_line(aes(y=0),linetype="dashed")+
  theme_bw()+
  theme(axis.text=element_text(size=30,face = "bold"))+
  labs(x="",y="")




simus=function(Eu.d=NULL,nsim=10,cov.model="exponential",cov.pars=c(0.5,1.5),
               Kappa=NULL,media=c(5,15),coord=NULL){
  
  res=list()
  datalist=list()
  for (i in 1:nsim) {
    datalist[[i]]=t(Datsim(Eu.d,cov.model = cov.model,cov.pars = cov.pars,
                           media = media,Kappa = Kappa))
    
  }
  cl <- parallel::makeCluster(detectCores())
  clusterEvalQ(cl, source("functions.R"))
  
  
  a<- parallel::parLapply(cl,datalist,kmeans.fdas,ncl = 4,Vmdist=TRUE,coord=coord,
                          cov.model=cov.model,Kappa=Kappa,multivgm=TRUE)
  
  # para sacar resultados mtv, tvm demutear b 
  # b<- parallel::parLapply(cl,datalist,kmeans.fdas,ncl = 4,Vmdist=FALSE,coord=NULL,
  #                         cov.model=cov.model,Kappa=Kappa,multivgm=TRUE)
  b<-a
  
  parallel::stopCluster(cl)
  for (i in 1:nsim) {
    c=cbind(tvm=a[[i]]$cluster,mtv=b[[i]]$cluster)
    c=data.frame(c)
    c$simu=i
    c$x=coord$c1
    c$y=coord$c2
    c$cuad=c(rep(3,30),rep(2,30),rep(4,30),rep(1,30))
    c=arrange(c,tvm)
    
    res[[i]]=c  
  }
  datat=do.call(rbind,res)
  
  z1=metod(datat,metodo = "tvm")
  z1=do.call(rbind,z1) %>% summarise(tvm=mean(perc))
  
  z2=metod(datat,metodo = "mtv")
  z2=do.call(rbind,z2) %>% summarise(mtv=mean(perc))
   
  return(data.frame(z1,z2))
}


# para multi: trabaja con vgm -> cov.model con Mayus 

system.time(
  k1<- simus(Eu.d,nsim = 100,cov.model = "Mat",cov.pars = c(6,2.5),
             Kappa = 2,media = c(5,15),coord = coord)
)
k1

