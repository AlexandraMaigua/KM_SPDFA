library(fda)
library(geoR)
library(fda.usc)
library(fpc)
library(factoextra)
source("functions.R")


MoranGfda<-function(b,weig.mat,data,coord,tipo="entre"){
  centros<-b$centers
  n<-length(b$centers$data[,1])
  if(tipo=="entre"){
      xmean<-func.mean(fdata(t(data)))
      ss<-c()
      s2<-0
      for (i in 1:n) {
        s=0
        for (j in 1:n) {
          num<-(centros[i,]-xmean$data)*(centros[j,]-xmean$data)
          
          integral<-int.simpson(fdata(num$data))
          numerador<-weig.mat[i,j]*integral
          s=s+numerador
      }
        
        ss[i]<-s   
        den<-(centros[i,]-xmean$data)^2
        integral2<-int.simpson(fdata(den$data))
        s2<-s2+integral2
      }  
      IMoran<-(n*sum(ss))/(sum(sum(weig.mat))*s2)
      }
  else{
    IMoran<-c()
    for (i in 1:length(centros$data[,1])) {
      centroide<-centros$data[i,]
      miembros<-which(b$cluster==i)
      miembros<-data[,miembros]  
      ss<-c()
      s2<-0
      #Adya
      coords<-coord[which(b$cluster==i),1:2]
      Eu.d <-as.matrix(dist(coords,method="euclidian"))
      weig.mat<-1/Eu.d
      diag(weig.mat)<-0
      
      for (j in 1:length(miembros[1,])) {
        s=0
        for (k in 1:length(miembros[1,])) {
          num<-(miembros[,j]-centroide)*(miembros[,k]-centroide)
          
          integral<-int.simpson(fdata(num))
          numerador<-weig.mat[j,k]*integral
          s=s+numerador
            
        }
          ss[j]<-s
          den<-(miembros[,j]-centroide)^2
          integral2<-int.simpson(fdata(den))
          s2<-s2+integral2
      }
      IMoran[i]<-(n*sum(ss))/s2
    }
    
  }
  
  
  return(IMoran)
}

adya<-c(0,1,0,0,0,1,0,0,1,0,0,0,0,1,1,0,1,1,0,0,0,0,1,0,0)
# adya<-matrix(1,ncol = 3,nr=3)
# diag(adya)<-0
adya<-matrix(adya,ncol=5,nrow=5)
Madya<-function(adya){
  for (i in 1:ncol(adya)) {
    adya[i,]<-adya[i,]/sum(adya[i,])
  }
  return(adya)
}
adya=Madya(adya)

#b=a
#
adya <-Ml2_dist(t(b$centers$data))
adya<-1/adya
diag(adya)<-0
#
MoranGfda(b,adya,data,coord,"entre")
MoranGfda(b,adya,data,coord,"o")

