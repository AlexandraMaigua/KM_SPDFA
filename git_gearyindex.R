library(fda)
library(geoR)
library(fda.usc)
library(fpc)
library(factoextra)
source("functions.R")

GearyGfda<-function(b,weig.mat,data,coord,tipo="entre"){
  centros<-b$centers
  n<-length(b$centers$data[,1])
  if(tipo=="entre"){
    xmean<-func.mean(fdata(t(data)))
    ss<-c()
    s2<-0
    for (i in 1:n) {
      s=0
      for (j in 1:n) {
        num<-(centros[i,]-centros[j,])^2
        
        integral<-int.simpson(fdata(num$data))
        numerador<-weig.mat[i,j]*integral
        s=s+numerador
      }
      
      ss[i]<-s   
      den<-(centros[i,]-xmean$data)^2
      integral2<-int.simpson(fdata(den$data))
      s2<-s2+integral2
    }  
    IGeary<-((n-1)*sum(ss))/(2*sum(weig.mat)*s2)
  }
  else{
    IGeary<-c()
    sss<-c()
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
          num<-(miembros[,j]- miembros[,k])^2
          
          integral<-int.simpson(fdata(num))
          numerador<-weig.mat[j,k]*integral
          s=s+numerador
          
        }
        ss[j]<-s
        den<-(miembros[,j]-centroide)^2
        integral2<-int.simpson(fdata(den))
        s2<-s2+integral2
      }
      sss[i]<-sum(ss)
      IGeary[i]<-((n-1)*sss[i])/(2*s2*sum(weig.mat))
    }
    
  }
  
  
  return(IGeary)
}

adya<-c(0,1,0,0,0,1,0,0,1,0,0,0,0,1,1,0,1,1,0,0,0,0,1,0,0)

boutG_entre<-GearyGfda(bout,adya,data,coord,"entre")

boutG_intra<-GearyGfda(bout,adya,data,coord,"e")

