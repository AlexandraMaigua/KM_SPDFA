library(ggmap)

shp=readOGR(dsn = "./data/ShapesRellenos/Shapefiles/demarcacion_hidrografica_250k_2014_geo.shp")

spdf <- spTransform(shp, CRS("+proj=longlat +datum=WGS84"))
spdf<- tidy(spdf) 

spdf=spdf %>% filter(long>=-81)
dh_simbol=as.character(shp@data[["dh_simbol"]])
dh_simbol=data.frame(dh_simbol,id=c(0:8))
spdf=merge(spdf,dh_simbol,by="id")

##

## para la presentación:
names(coord)[3]="Cluster" # "Grupo"
names(spdf)[8]="Basins"  # "Cuenca"


# Plot it
cols=c("#FF0000","#00FFFF")

spx1=ggplot() +  
  geom_point(data = coord %>% filter(Cluster %in%  c(1:2)),aes(x,y,color=Cluster),size=4,alpha=1)+
  scale_color_manual(values=cols)+
  geom_polygon(data = spdf, aes( x = long, y = lat, group = group,fill=Basins),alpha=0.14,color="black")+#fill=dh_simbo
  scale_fill_viridis_d(option = "turbo",direction = -1,begin = 0,end = 1)+
  labs(title="Classification",subtitle = 'Hydrographical Demarcation',x="",y="")+theme_bw()+
  theme(axis.text=element_text(size=20),
        panel.grid = element_blank(),
        #plot.subtitle = element_text(face = "bold",size = 15),
        plot.title = element_text(face = "bold",size = 18),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=10))+
  geom_label(data=dfc,aes(label=cuenca,x=x,y=y,fontface=2))+
  scale_y_continuous(breaks = c(2:-5),
                     labels = c("2N","1°N","0°","1°S","2°S","3°S","4°S","5°S"))+
  scale_x_continuous(breaks = c(-81:-75),
                     labels = c("81°W","80°W","79°W","78°W","77°W","76°W","75°W"))

cols=c("#FF00FF","#00FF00")

spx2=ggplot() +  
  geom_point(data = coord %>% filter(Cluster %in%  c(3:4)),aes(x,y,color=Cluster),size=4,alpha=1)+
  scale_color_manual(values=cols)+
  geom_polygon(data = spdf, aes( x = long, y = lat, group = group,fill=Basins),alpha=0.14,color="black")+#fill=dh_simbo
  scale_fill_viridis_d(option = "turbo",direction = -1,begin = 0,end = 1)+
  labs(title="Classification",subtitle = 'Hydrographical Demarcation',x="",y="")+theme_bw()+
  theme(axis.text=element_text(size=20),
        panel.grid = element_blank(),
        #plot.subtitle = element_text(face = "bold",size = 15),
        plot.title = element_text(face = "bold",size = 18),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=10))+
  geom_label(data=dfc,aes(label=cuenca,x=x,y=y,fontface=2))+
  scale_y_continuous(breaks = c(2:-5),
                     labels = c("2N","1°N","0°","1°S","2°S","3°S","4°S","5°S"))+
  scale_x_continuous(breaks = c(-81:-75),
                     labels = c("81°W","80°W","79°W","78°W","77°W","76°W","75°W"))

cols=c("#FF7F00")

spx3=ggplot() +  
  geom_point(data = coord %>% filter(Cluster %in%  c(5)),aes(x,y,color=Cluster),size=4,alpha=1)+
  scale_color_manual(values=cols)+
  geom_polygon(data = spdf, aes( x = long, y = lat, group = group,fill=Basins),alpha=0.14,color="black")+#fill=dh_simbo
  scale_fill_viridis_d(option = "turbo",direction = -1,begin = 0,end = 1)+
  labs(title="Classification",subtitle = 'Hydrographical Demarcation',x="",y="")+theme_bw()+
  theme(axis.text=element_text(size=20),
        panel.grid = element_blank(),
        #plot.subtitle = element_text(face = "bold",size = 15),
        plot.title = element_text(face = "bold",size = 18),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=10))+
  geom_label(data=dfc,aes(label=cuenca,x=x,y=y,fontface=2))+
  scale_y_continuous(breaks = c(2:-5),
                     labels = c("2N","1°N","0°","1°S","2°S","3°S","4°S","5°S"))+
  scale_x_continuous(breaks = c(-81:-75),
                     labels = c("81°W","80°W","79°W","78°W","77°W","76°W","75°W"))

