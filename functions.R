source("fit.lmc.R")
source("parameters.lmc.R")
source("covariances.R")
source("variograms.R")
source("multiv.R")
library(sp)
library(operators)
library(fda.usc)
library(mvtnorm)
library(sp)
library(gstat)
library(geoR)
library(ggplot2)
library(plyr)
library(parallel)
library(data.table)
library(geofd)



Ml2_dist<-function(data,nbasis=65)
{
  nbasis<-nbasis
  n <-  dim(data)[1]
  s <-  dim(data)[2] 
  argvals<-seq(1,n, by=1)
  range  <- c(1,n)
  period <- n
  basis <- create.fourier.basis(range, nbasis, period)
  datafd <- Data2fd(data,argvals,basis)  
  L2norm<-matrix(0,nrow=s,ncol=s)
  coef<-datafd$coef 
  M<-fourierpen(basis,Lfdobj=0)
  ## se puede parallel/ no se obtuvo mejora
  for (i in 1:(s-1))
  {
    coef.i<-coef[,i]
    for (j in (i+1):s)
    {
      coef.j<-coef[,j]
      L2norm[i,j]<-t(coef.i-coef.j)%*%M%*%(coef.i-coef.j)
      L2norm[j,i]<-L2norm[i,j]
    }
  }
  names<-names(data.frame(data))
  dimnames(L2norm)<-list(names,names)
  return(L2norm)
}
library(gdata)

kmeans.center.iniL2<-function (fdataobj,Vmdist=FALSE,coord=NULL, ncl = 2,
                               cov.model="spherical", Kappa=NULL,multivgm=multivgm,
                               method = "sample", max.iter = 10000, max.comb = 1e+06, ...) 
{
  if (!is.fdata(fdataobj)) 
    fdataobj = fdata(fdataobj)
  
  z <- fdataobj[["data"]]
  if(!isTRUE(Vmdist) && is.null(coord))
  {
    mdist<-Ml2_dist(t(z),nbasis = 65) #cambiar nbasisi para data en funcion
  }
  else 
  {
    mdist<-MVari(t(z),coord,nugget.fix=NULL, max.dist.variogram=NULL,
                 cov.model = cov.model,Kappa = Kappa,multivgm = multivgm)
    print("entro al variograma con el cov.model")
    print(cov.model)
  }
  
  tt <- fdataobj[["argvals"]]
  rtt <- fdataobj[["rangeval"]]
  names <- fdataobj[["names"]]
  nr = nrow(fdataobj)
  nc = ncol(fdataobj)
  if (is.vector(ncl)) {
    len.ncl = length(ncl)
    ngroups = ncl
    max.combn <- choose(nr, ngroups)
    if (len.ncl == 1) {
      ind = 1
      if (method == "exact") {
        if (max.combn > max.comb) 
          warning(paste0(max.combn, " samples are required, it has been limited to a random sample of size ", 
                         max.comb))
        method = "sample"
      }
      if (method == "sample") {
        max.iter <- min(max.iter, max.combn)
        vec <- array(NA, dim = c(ngroups, max.iter))
        vec.d <- rep(NA, nc)
       
        for (i in 1:max.iter) {
          #i=1
          vec[, i] <- sample(1:nr, ngroups, replace = FALSE)#Tipo Macqueen / centros ini
          
          vec.d[i] <- sum(mdist[vec[, i], vec[, i]])
        }
        ind.max <- which.max(vec.d)# centros más separados
        lxm <- vec[, ind.max]# indices de los centros
      }
      else if (method == "exact") {
        co <- combn(1:nr, ngroups)
        nco <- ncol(co)
        vec <- rep(NA, nco)
        ## SE PUEDE Parallel/no se obtuvo mejora
        for (i in 1:nco) {
          vec[i] <- sum(mdist[co[, i], co[, i]])
        }
        max.vec <- which.max(vec)
        lxm <- co[, max.vec]
      }
      else stop("Center initialization method unknown")
      xm = z[lxm, ]# curvas de los centros
    }
    else stop("Argument 'ncl' is expected the number of groups to detect")
  }
  else stop("Argument 'ncl' is expected the number of groups to detect")
  d = rbind(mdist, mdist[lxm, ])
  centers = fdata(xm, tt, rtt, names)
  if (isTRUE(draw)) {
  #aumentar codigo para gráficos
  }
  out <- list(centers = centers, lcenters = lxm, z.dist = mdist, 
              fdataobj = fdataobj)
  class(out) <- "kmeans.fd"
  return(invisible(out))
}

library(geoR)
library(gstat)


# añadir más parámetros a la función
MVari<-function(data,coord,nugget.fix=NULL, max.dist.variogram=NULL,
                multivgm=multivgm ,nbasis=65,cov.model="Gau",Kappa=NULL)
{

  dia=1:length(data[,1])
  fourier.basis=create.fourier.basis(
    rangeval = range(dia), nbasis = nbasis
  )
  
  temp2fd=Data2fd(argvals = dia,
                  y=data,
                  basisobj = fourier.basis)
  
  

  M2d=dist(t(temp2fd$coefs))^2
  
  
  
  Eu.d <-as.matrix(dist(coord,method="euclidian"))
 
  if(isFALSE(multivgm)){
   
     emp.trace.vari<-variog(coords=coord, data=apply(temp2fd$coefs,2,sum),option="bin",message=FALSE,
                            trend = "cte",max.dist =2.8)#poner max.dist, trend=cte-cnd #para temp tnd=1st mxd=3
                                                        #maxdist=2.8 mejores resultados en simus out
    
    
    datap=cbind(coord,valor=apply(data, 2,sum))
    
    
    plot(as.geodata(datap),trend="2nd",cex.lab=2)
    
    axis(1, col = "white", col.ticks = "grey61", 
         lwd.ticks = 0.5, tck = -0.025, cex.axis = 2, col.axis = "black")
   
    
     axis(2, col = "white", col.ticks = "grey61", 
         lwd.ticks = 0.5, tck = -0.025, cex.axis = 2, col.axis = "black")
    grid(col = "grey75", lwd = 0.3)
   
    
   
    if (is.null(max.dist.variogram))
     
      if(is.null(Kappa)){
     
      data=cbind(coord,value=apply(temp2fd$coefs,2,sum))
      coordinates(data) = ~x+y
      emp.trace.vari=variogram(value~1, data,dX=0,cutoff=3)
      plot(emp.trace.vari)
      sigma2.0=quantile(emp.trace.vari$gamma,0.75)
      phi.0=quantile(emp.trace.vari$dist,0.2)
      nt=mean(emp.trace.vari$gamma)/4
      trace.vari = fit.variogram(emp.trace.vari, vgm(sigma2.0, cov.model, phi.0,nt))
      print("con vgm")
      g=plot(emp.trace.vari,trace.vari)
      show(g)
      
      
    }
    else{
     
      data=cbind(coord,value=apply(temp2fd$coefs,2,mean))
      coordinates(data) = ~x+y
      emp.trace.vari=variogram(value~1, data,dX=0,cutoff=3.5)
      plot(emp.trace.vari)
      sigma2.0=quantile(emp.trace.vari$gamma,0.75)
      phi.0=quantile(emp.trace.vari$dist,0.2)
      nt=mean(emp.trace.vari$gamma)/4
      trace.vari = fit.variogram(emp.trace.vari, vgm(psill =  sigma2.0,model =  cov.model,
                                                     range = phi.0,nugget = nt,kappa = Kappa),
                                 fit.kappa = TRUE)
      g=plot(emp.trace.vari,trace.vari)
      show(g)
      
      print("matern con kappa")
    }
    
   
    sigma2=trace.vari$psill[2]
    nugget=trace.vari$psill[1]
    tra.vari.mat <- sigma2+nugget - cov.spatial(Eu.d,cov.model=tolower(trace.vari$model[2]),
                                                           cov.pars=c(sigma2,trace.vari$range[2]),
                                                           kappa=trace.vari$kappa[2])

   
    weig.mat<-sqrt(as.matrix(M2d))*tra.vari.mat/1000 #revisar 
  } 
  else{
    
    multvariogram=multiv(Eu.d,coord,temp2fd$coefs,cov.model,3.5,Kappa = Kappa)
    L2norm<-M2d
    
    
    
    weig.mat<-sqrt(as.matrix(L2norm))*multvariogram/1000 #revisar
    print("matriz multivariograma")
  }
  
  return(weig.mat)
  
}

library(fda)
library(fda.usc)


kmeans.centers.update=function(out,group
                               ,dfunc=func.trim.FM,draw=FALSE
                               ,par.dfunc=list(trim=0.05)
                               ,...){
  if (class(out)!="kmeans.fd") 
    stop("Error: incorrect input data")
  z = out$fdataobj[["data"]]
  tt = out$fdataobj[["argvals"]]
  rtt <- out$fdataobj[["rangeval"]]
  names = out$fdataobj[["names"]]
  mdist = out$z.dist
  centers = out$centers
  xm = centers[["data"]]
  nr = nrow(z)
  nc = ncol(z)
 
  grupo = group #group=c(1:5) 

  ngroups = length(unique(group))
 
  d = out$d
  ncl = nrow(xm)
  for (j in 1:ngroups){
   jgrupo <- grupo==j
    dm=z[jgrupo,]
    ind=which(jgrupo)
    if (is.vector(dm) || nrow(dm)==1) {
      k=j
      stat=dm
    }
    else   {
      par.dfunc$fdataobj<-centers
      par.dfunc$fdataobj$data<-dm
      stat=do.call(dfunc,par.dfunc)
    }
    if (is.fdata(stat)) xm[j,]=stat[["data"]]
    else  xm[j,]=stat
    
  }#hasta aqui actualiza los centros por la media
  centers$data=xm
  row.names(centers$data) <- paste("center ",1:ngroups,sep="")
  if (isTRUE(draw)){
    }
  return(list("centers"=centers,"cluster"=grupo))
}


kmeans.assig.groups=function(out,draw=FALSE,...){
  if (!is.null(out$lcenters))  
    lxm=out$lcenters else   lxm=NULL
    nr = nrow(out$fdataobj)
    nc = ncol(out$fdataobj)
    xm = out$centers[["data"]]
    grupo = rep(0,nr)
    d = out$d
    ngroups=nrow(d)-nrow(out$fdataobj[["data"]])
    grupo <- apply(d[(nr+1):(nr+ngroups),],2,which.min)#asignacion de grupos
    if (isTRUE(draw)){
      if (nr==2){
        plot(out$fdataobj,main="Assigning groups")
        for (i in 1:ngroups){points(xm[i,1],xm[i,2],col=i+1,pch=8,cex=1.5)}
      }
      else{
        plot(out$fdataobj,col=grupo+1,main="Assigning groups",lwd=.3,lty=2)
        lines(out$centers,col=2:(length(grupo+1)),lwd=3,lty=1)    
      }
    }
    if (nc==2) { for (j in 1:nc){points(xm[j,1],xm[j,2],col=j+1,pch=7,cex=1.5)}}
    return(list("centers"=out$centers,"cluster"=grupo))
}


kmeans.fdM<-function (fdataobj, ncl = 2, metric = metric.lp,Vmdist=FALSE,coord=NULL,
                      cov.model="spherical", Kappa=NULL, multivgm=multivgm,
                      dfunc = func.trim.FM, 
                      max.iter = 10000, par.metric = NULL, par.dfunc = list(trim = 0.05), 
                      method = "sample", cluster.size = 1, draw = FALSE, ...) 
{
  if (!is.fdata(fdataobj)) 
    fdataobj = fdata(fdataobj)
  nas1 <- is.na(fdataobj)
  if (any(nas1)) 
    stop("fdataobj contain ", sum(nas1), " curves with some NA value \n")
  z <- fdataobj[["data"]]
  tt <- fdataobj[["argvals"]]
  rtt <- fdataobj[["rangeval"]]
  nr = nrow(z)
  nc = ncol(z)
  if (is.vector(ncl)) {
    len.ncl = length(ncl)
    if (len.ncl == 1) {
      par.ini <- list()
      par.ini$fdataobj = fdataobj
      par.ini$method = method
      par.ini$ncl = ncl
      par.ini$draw = draw
      par.ini$max.comb = 1e+06
      par.ini$max.iter = max.iter
      par.ini$Vmdist = Vmdist
      par.ini$coord = coord
      par.ini$cov.model=cov.model
      par.ini$Kappa=Kappa
      par.ini$multivgm=multivgm
      if (!is.null(par.metric)) 
        par.ini$par.metric <- par.metric
      par.ini$... <- par.metric
      out1 = do.call(kmeans.center.iniL2, par.ini)
      lxm <- out1$lcenters 
      out1$d = rbind(out1$z.dist, out1$z.dist[lxm, ])
     }
    else {
      ngroups <- length(ncl)
      lxm <- ncl
      xm <- z[lxm, ]
      out1 <- list()
      out1$fdataobj <- fdataobj
      out1$ncl <- len.ncl
      if (is.null(par.metric)) 
        par.metric <- list(p = 2, w = 1)
      par.metric$fdata1 <- fdataobj
       if(!isTRUE(Vmdist) && is.null(coord))
      {
        mdist<-Ml2_dist(t(z))
      }
      else 
      {
        mdist<-MVari(t(z),coord,nugget.fix=NULL, max.dist.variogram=NULL,cov.model = cov.model) 
        
      }
      out1$z.dist <- mdist
      out1$d <- rbind(mdist, mdist[lxm, ])
      out1$centers <- fdataobj[ncl, ]
      out1$lcenters <- ncl
      class(out1) <- "kmeans.fd"
    }
  }
  else if (is.fdata(ncl)) {
    lxm = NULL
    xm = ncl[["data"]]
    if (is.null(par.metric)) 
      par.metric = list(p = 2, w = 1)
    par.metric$fdata1 <- fdataobj
    if(!isTRUE(Vmdist) && is.null(coord))
    {
      mdist<-Ml2_dist(z)
    }
    else 
    {
      mdist<-MVari(t(z),coord,nugget.fix=NULL, max.dist.variogram=NULL)
    }
    par.metric2 <- par.metric
    par.metric2$fdata2 <- ncl
     if(!isTRUE(Vmdist) && is.null(coord))
    {
      mdist<-Ml2_dist(z)
    }
    else 
    {
      mdist<-MVari(t(z),coord,nugget.fix=NULL, max.dist.variogram=NULL)
    }
    out1 = list()
    out1$fdataobj <- fdataobj
    out1$centers = ncl
    out1$lcenters <- NULL
    ngroups = nrow(ncl)
    ncl = nrow(ncl)
    out1$d = rbind(mdist, t(mdist2))
    class(out1) = "kmeans.fd"
  }
  ngroups = nrow(out1$centers[["data"]])
  a = 0
  aa <- i <- 1
  same_centers = FALSE
  if (is.null(colnames(out1$d))){ 
    cnames <- colnames(out1$d) <- 1:NCOL(out1$d)}
  else{ cnames <- colnames(out1$d)}
  while ((i < max.iter) && (!same_centers)) {
    iterar <- FALSE
    out3 = kmeans.assig.groups(out1, draw = draw)
    names(out3$cluster) <- cnames
    tab <- table(out3$cluster)
    imin <- which.min(tab)[1]
    if (cluster.size > tab[imin]) {
      warning(paste0(" One of the clusters only has ", 
                     tab[imin], " curves and the minimum cluster size is ", 
                     cluster.size, ".\n The cluster is completed with the closest curves of the other clusters."))
      iclust <- out3$cluster == imin
      dist.aux <- out1$d[imin, ]
      icambios <- as.numeric(names(sort(out1$d[imin, ])[1:cluster.size]))
      out1$d[imin, icambios] <- 0
      out1$z.dist <- out1$d
      out3$cluster[icambios] <- imin
      out2 <- out3
      out1$cluster <- out3$cluster
      par.dfunc$fdataobj <- fdataobj[c(icambios)]
      out1$centers[imin] = do.call(dfunc, par.dfunc)
      iterar <- TRUE
      i = i + 1
    }
    out2 = kmeans.centers.update(out1, group = out3$cluster, 
                                 dfunc = dfunc, draw = draw, par.dfunc = par.dfunc)
    if (!iterar) {
      same_centers <- out2$centers$data == out3$centers$data
      out1$centers <- out2$centers
      i = i + 1
    }
  }
  out <- list(cluster = out2$cluster, centers = out2$centers)
  return(out)
}


kmeans.fdas<-function(fdataobj, ncl = 2,Vmdist=FALSE,coord=NULL,
                      cov.model=cov.model,Kappa=NULL,
                      multivgm=FALSE){
  if(!isTRUE(Vmdist) && is.null(coord))
  {
    return(kmeans.fd(fdataobj,ncl = ncl,draw = TRUE,cluster.size = 1))
  }
  else
  {
    return(kmeans.fdM(fdataobj,ncl = ncl,Vmdist=Vmdist,coord=coord,cov.model = cov.model,
                      Kappa=Kappa,multivgm=multivgm ))
  }
}




## necesita la data para correr, exportar data al nucleos
SSBindex1<-function(a)
{
  centros<-t(a$centers$data)
  xmean<-func.mean(fdata(t(data))) 
  centros<-cbind(centros,t(xmean$data))
  
  distancias<-Ml2_dist(centros,nbasis = 65)
 
  njs<-data.frame(table(a$cluster))
  n1<-length(distancias[,1])
  dnj<-sum(distancias[n1,1:(n1-1)]*njs$Freq)/length(data) 
  return(dnj)
}


SSWindex<-function(a,data)
{
  s1<-c()
  for(i in 1:max(a$cluster))
  {
    #i=1
    b<-which(a$cluster==i)
    c<-data[,b]
    cj<-a$centers$data[i,]
    c<-cbind(c,(cj))
    
    distancias<-Ml2_dist(c,5)
    n1<-length(distancias[,1])
    dnj<-sum(distancias[n1,1:(n1-1)])
    s1[i]<-dnj
  }
  ss<-sum(s1)/length(data) 
  return(ss)
}


Siluetindex<-function(a){
  si<-list()
  for(k in 1:length(a$centers))
  {
    
    b<-which(a$cluster==k)
    if(length(b)==1){
      si[[k]]<-1
    }
    else{
      c<-data[,b]
      distancias<-Ml2_dist(c)
      ai<-c()
      for (i in 1:length(b)) {
        
        ai[i]<-sum(distancias[i,])/(length(distancias[i,])-1)  
      }  
      distcent<-Ml2_dist(t(a$centers$data))
      
      
      d<-which(min(distcent[k,-k])==distcent[k,])
      b1<-which(a$cluster==d[1])
      c1<-data[,b1]
      
      cc1<-cbind(c,c1)
      distancias2<-Ml2_dist(cc1)
      bi<-c()
      for (i in 1:length(ai)) {
        
        bi[i]<-sum(distancias2[i,-(1:length(ai))])/(length(distancias2[i,-(1:length(ai))]))  
      }
      si[[k]]<-(bi-ai)/max(ai,bi)
    }
  }
  return(si)
}


Dunnindex<-function(a){
  
  centros<-a$centers$data
  distcent<-Ml2_dist(t(centros),15)
  dk<-c()
  for(i in 1:length(a$centers)){
    i=1
    b<-which(a$cluster==i)
    if(length(b)==1){
      dk[i]<-1
    }
    else{
      c<-data[,b]
      distancias<-Ml2_dist((c),15) #estaba sin t
      dk[i]<-max(distancias) 
    }
  }
  mdk<-max(dk)
  djt<-c()
  for (j in 1:length(a$centers)) {
    for (t in 1:length(a$centers)) {
      if(j!=t){
        djt[j]<-min(distcent[j,t]/mdk)
      }
    }
    
  }
  return(min(djt))
}


#------------------------------------------------

library(dplyr)

metod=function(datat,metodo="tvm"){

  metodo=ifelse(metodo=="tvm","mtv","tvm")

  nsin=unique(datat$simu)
  res=list()
  for (i in nsin) {
    z0=datat %>% dplyr::filter(simu==i) %>% dplyr::select(-metodo)
    names(z0)[1]="Grupo"

    z1=z0 %>% 
      dplyr::group_by(simu,Grupo,cuad) %>% dplyr::summarise(n=n()) %>%
      mutate(perc=round(n/30*100,2)) %>%
      dplyr::group_by(simu,cuad) %>% dplyr::summarise(perc=max(perc))


    z1.1=z0 %>% 
      dplyr::group_by(simu,Grupo,cuad) %>% dplyr::summarise(n=n()) %>%
      mutate(perc=round(n/30*100,2)) %>%
      dplyr::group_by(cuad)

    z1.2=merge(z1,z1.1,by=c("cuad","perc"))
    z1.2=z1.2[,-c(3,4)]
    z1.2=z1.2 %>% distinct(cuad,.keep_all = TRUE)

    z1.2.1=z1.2 %>% group_by(Grupo) %>% summarise(perc=max(perc))
    z1.2.1=merge(z1.2.1,z1.2,by=c("Grupo","perc")) %>% distinct(Grupo,.keep_all = TRUE)
    
    inde=z1.2.1$cuad
    indg=z1.2.1$Grupo
    
    ind1=which(1:4 %!in% inde)
    ind2=which(1:4 %!in% indg)
    
    if(sum(1:4 %in% z1.2$Grupo)<4){
       for (j in ind2) {
        w=z1.1 %>% dplyr::filter(Grupo==j)
        w=w[,-1]
        w=w %>% filter(cuad %in% ind1)
        ind1.1=which(w$perc==max(w$perc))  
        w=w[ind1.1[1],]
        ind1=ind1[-which(ind1==w$cuad)]
        w=w[,c(2,4,1,3)]
        ind2.1=which(w$cuad==z1.2$cuad)
        z1.2[ind2.1,]=w
      }
      res[[i]]=z1.2[,1:2]
      print("arriba")
    }
    else{
      z1.2= z1.1 %>% dplyr::group_by(cuad) %>% dplyr::summarise(perc=max(perc))
      res[[i]]=z1.2
      print("abajo")
    }
  }  
  return(res)
}


Datsim=function(Eu.d,cov.model="exponential",cov.pars=c(0.5,1.5),Kappa=NULL,
                media=c(5,6)){
  if(is.null(Kappa)){
    cov.model=tolower(cov.model)
     Sigma1=cov.spatial(Eu.d,cov.model = cov.model,cov.pars = cov.pars)
     cuad1<- rmvnorm(365,mean=c(rep(media[1],30)), sigma=Sigma1[91:120,91:120])
     cuad2<- rmvnorm(365,mean=c(rep(media[2],30)), sigma=Sigma1[31:60,31:60])
     cuad3<- rmvnorm(365,mean=c(rep(media[1],30)), sigma=Sigma1[1:30,1:30])
     cuad4<- rmvnorm(365,mean=c(rep(media[2],30)), sigma=Sigma1[61:90,61:90])

     sim1<-cbind(cuad3,cuad2,cuad4,cuad1)
   }
  else{
    cov.model=tolower(cov.model)
    Sigma1=cov.spatial(Eu.d,cov.model = cov.model,cov.pars = cov.pars,kappa = Kappa)
    sim1=rmvnorm(365,mean=c(rep(media[1],30),rep(media[2],30),
                            rep(media[2],30),rep(media[1],30)), sigma=Sigma1)
  }
  return(sim1)
}



Datsim1=function(Eu.d=NULL,cov.model="exponential",cov.pars=c(0.5,1.5),Kappa=NULL,
                media=c(5,15)){
  if(is.null(Kappa)){
    Sigma1=diag(120)
    sim1=rmvnorm(365,mean=c(rep(media[1],30),rep(media[2],30),
                            rep(media[2],30),rep(media[1],30)), sigma=Sigma1)
  }
  else{
    Sigma1=cov(matrix(rnorm(14400,0,1),nrow = 120,ncol = 120))
    sim1=rmvnorm(365,mean=c(rep(media[1],30),rep(media[2],30),
                            rep(media[2],30),rep(media[1],30)), sigma=Sigma1)
  }
  return(sim1)
}


