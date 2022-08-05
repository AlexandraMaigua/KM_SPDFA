source("functions.R")

#library(npsp) github
library("rnaturalearth")
library("rnaturalearthdata")
library(rgeos)
world <- ne_countries(scale = "medium", returnclass = "sf")



load("./1_Datos_MF/data.RData")
load("./1_Datos_MF/coord.RData")

data=t(data)
#Eliminación de atípicos

library(fdaoutlier)

out=msplot(t(data))

data=data[,-out$outliers]

coord=coord[-out$outliers,]
coord=coord[,1:2]

#Estimación del variograma a ojo
emp.trace.vari<-variog(coords=coord, data=apply(data,2,sum),option="bin",message=FALSE,
                       trend = "2nd",max.dist = 2)
plot(emp.trace.vari)
eyefit(emp.trace.vari)

b=kmeans.fdas(fdataobj = t(data),ncl = 5,Vmdist = T,coord = coord,
              cov.model = "wave",Kappa = NULL,multivgm = FALSE)

load("bout.Rdata")

b=bout
coord$Grupo=as.factor(b$cluster)

#Visualización de grupos

library(ggplot2)
library(sf)
g<-ggplot(data = world) +
  geom_sf() +
  geom_point(data = coord,aes(x,y, color=Grupo),size=1 )+
  coord_sf(xlim = c(-81, -75), ylim = c(-5, 1.5))+ 
  labs(title = "clasificacion k-means")
library(plotly)
ggplotly(g)
ggplot(data = coord)+
  geom_point(aes(x=x,y=y))


