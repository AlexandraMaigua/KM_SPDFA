library(rgdal)
library(broom)
library(plotly)
library("rgdal")
library("rgeos")
library(viridis)
library(viridisLite)
library(readr)
library(npsp)
library(ggmap)
library(ggthemes)
source("functions.R")

shp=readOGR(dsn = "./data/ShapesRellenos/Shapefiles/demarcacion_hidrografica_250k_2014_geo.shp")

spdf <- spTransform(shp, CRS("+proj=longlat +datum=WGS84"))
spdf<- tidy(spdf) 


## con shp hidrografica/data

spdf=spdf %>% filter(long>=-81)
dh_simbol=as.character(shp@data[["dh_simbol"]])
dh_simbol=data.frame(dh_simbol,id=c(0:8))
spdf=merge(spdf,dh_simbol,by="id")

##

## para la presentación:
names(coord)[3]="Cluster" # "Grupo"
names(spdf)[8]="Basins"  # "Cuenca"

##data para nombres de cuencas

## y para anombres  en git_pixeles.R, para spxi

cuenca=c("DHESMERALDAS","DHGUAYAS","DHJUBONES","DHMANABÍ","DHMIRA","DHNAPO",
         "DHPASTAZA","DHPUYANGO-CATAMAYO","DHSANTIAGO")
x=c(-79.3,-79.8,-79.7,-80.3,-78,-77.5,-77,-79.7,-78)
y=c(0.5,-1.8,-3.3,-1,0.7,-0.5,-2,-4.1,-2.7)
dfc=data.frame(cuenca,x,y)

cols=c("#FF0000","#00FFFF","#FF00FF","#00FF00","#FF7F00")

#ggg=
spxg=ggplot() +
  geom_point(data = coord ,aes(x,y,color=Cluster),alpha=1,size=3)+
  scale_color_manual(values=cols)+
  geom_polygon(data = spdf, aes( x = long, y = lat, group = group,fill=Basins),
               alpha=0.14,color="black")+
  scale_fill_viridis_d(option = "turbo",direction = -1,begin = 0,end = 1)+ ## 0.4-0.7 
  labs(x="",y="",title = "Classification",
       subtitle="Hydrographical Demarcation"
  )+theme_bw()+
  theme(axis.text=element_text(size=20),
        panel.grid = element_blank(),
        #plot.subtitle = element_text(face = "bold",size = 15),
        plot.title = element_text(face = "bold",size = 18),
        legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))+
  geom_label(data=dfc,aes(label=cuenca,x=x,y=y,fontface=2))+
  scale_y_continuous(breaks = c(2:-5),
                     labels = c("2N","1°N","0°","1°S","2°S","3°S","4°S","5°S"))+
  scale_x_continuous(breaks = c(-81:-75),
                     labels = c("81°W","80°W","79°W","78°W","77°W","76°W","75°W"))

spxg

ggarrange(spxg1,spxg,ncol = 2,nrow = 1) 

# mapa de calor sobre shape
# ----------------------------------
# ----------------------------------


shp=readOGR(dsn = "./data/cantones/nxprovincias.shp")
spdf <- spTransform(shp, CRS("+proj=longlat +datum=WGS84"))
spdf<- tidy(spdf) 

spdf=spdf %>% filter(long>=-81)

# con data2 de aplderiv.R

datah=data2[,c(2:4)]
names(datah)[3]='NDVI'

cities=c("QUITO","GUAYAQUIL","CUENCA","TENA","RIOBAMBA","LATACUNGA")
x=c(-78.32495,-79.88,-79.15,-77.81286,-78.64712,-78.51554)
y=c(-0.001,-1.97,-2.96,-0.68,-1.87,-0.81)
dfc=data.frame(cities,x,y)

library(ggrepel)
heat<-ggplot() +
  geom_polygon(data = spdf, aes( x = long, y = lat, group = group,fill=id),
               alpha=0.1,color="black")+
  guides(fill = FALSE)+
  geom_point(data = datah ,aes(x,y,color=NDVI),alpha=2,size=2)+
  scale_color_gradientn(colours = rainbow(12,start =0.3,end =1,s = 1,v = 1),values = c(1,0))+
  scale_fill_viridis_d(option = "turbo",direction = -1,begin = 0,end = 1)+ ## 0.4-0.7 
  labs(x="",y="")+#"Clasification: Hydrographical Demarcation"  )+
  theme_bw()+
  theme(axis.text=element_text(size=15),
        panel.grid = element_blank(),
        plot.subtitle = element_text(face = "bold",size = 15),
        plot.title = element_text(face = "bold",size = 18))+
  geom_label_repel(data=dfc,aes(label=cities,x=x,y=y,fontface=2),alpha=0.7)+
  scale_y_continuous(breaks = c(2,1,0,-1,-2,-3,-4),labels = c(2,1,0,-1,-2,-3,-4))

heat


