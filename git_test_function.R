library(tidyverse)
library(TSclust)
library(ggpubr)
library(grid)
library(gridExtra)
library(tidyverse)
source("functions.R")

load("bout.Rdata")
b=bout

### para hacer el arrange---------------
a=b
rec=length(a$centers)
#colores=c(2,5,2,5,2)
colores=c("#FF0000","#00FFFF","#FF00FF","#00FF00","#FF7F00")

p=list()
for(i in 1:rec){
  b<-which(a$cluster==i)
  c<-data[,b]
  centros<-a$centers$data 
  c=data.frame(c) 
  c1=gather(data = c,estacion, valor)
  c1$dia=rep(c(1:365),length(c[1,])) 
  
  media=data.frame(centros[i,])
  media$dia=rep(c(1:365),1)
    names(media)[1]="valor"
  
  titulo=paste("Cluster",sep = " ",i)
  
  
  p[[i]]=ggplot(c1)+geom_line(aes(x=dia,y=valor,coluor=estacion),alpha=0.2)+
    geom_line(data=media,aes(x=dia,y=valor),color=colores[i],size=2)+
    theme_bw(base_size = 15)+
    theme(axis.text = element_text(size=20))+
    labs(title = titulo,x="Day",y="NDVI")+
    scale_x_continuous(breaks = c(1,100,200,300),labels = c(1,100,200,300))
}

arrangeGrob(spx2,arrangeGrob(p[[3]],p[[4]],ncol = 1),ncol=2)%>% 
  as_ggplot()
arrangeGrob(spx3,arrangeGrob(p[[5]],ncol = 1),ncol=2)%>% 
  as_ggplot()

## spx sale de git_pixeles.R
lay=rbind(c(1,1,1,2,2),
          c(1,1,1,2,2),
          c(1,1,1,2,2))
g1<-arrangeGrob(spx1,arrangeGrob(p[[1]],p[[2]],ncol = 1),ncol=2,layout_matrix = lay) %>% 
  as_ggplot()
g2<-arrangeGrob(spx2,arrangeGrob(p[[3]],p[[4]],ncol = 1),ncol=2,layout_matrix = lay) %>% 
  as_ggplot()
g3<-arrangeGrob(spx3,arrangeGrob(p[[5]],ncol = 1),ncol=2,layout_matrix = lay) %>% 
  as_ggplot()

#---------------------

#correlaciones



## coors bout
mcc9=c(0.7698621, 0.8172329, 0.8075263, 0.7936363, 0.8694921)
ccorem=c(0.5932437, 0.6689656, 0.6521596, 0.6302183, 0.7561212)

# kmean-mod
cgeary=c(0.0050331671, 0.0028864678, 0.0009766287, 0.0021945645, 0.0028154862)
cmoran=c(-0.64) ## abajo esa el moran incluido ggplot

ccm=data.frame(Grupo=factor(c(1:5)),Cor=round(mcc9,2))
ccm2=data.frame(Grupo=factor(c(1:5)),Cor=round(ccorem,2))
cgry=data.frame(Grupo=factor(c(1:5)),Cor=round(cgeary,5))
cmrn=data.frame(Algo=c("Moran"),Moran=cmoran)

cols=c("#FF0000","#00FFFF","#FF00FF","#00FF00","#FF7F00")

gcorr=ggplot(ccm,aes(x=Grupo,y=Cor,fill=factor(Grupo)))+
  geom_bar(stat = "identity",position = "dodge")+
  scale_fill_manual(values = cols)+
  geom_text(aes(label=Cor), position=position_dodge(width=0.9), vjust=-0.25)+theme_bw()+
  theme(axis.text = element_text(size=15,face = "bold"), legend.position = "none",
        axis.title = element_text(size=15),
        plot.title = element_text(size=15))+
  ylim(c(0,1))+
  labs(y="Correlation",x="Cluster",title = "Correlation: centroid-members")

gcorr1=ggplot(ccm1,aes(x=Grupo,y=Cor,fill=factor(Grupo)))+
  geom_bar(stat = "identity",position = "dodge")+
  scale_fill_manual(values = cols)+
  geom_text(aes(label=Cor), position=position_dodge(width=0.9), vjust=-0.25)+
  theme(axis.text = element_text(size=15,face = "bold"), legend.position = "none",
        axis.title = element_text(size=15),
        plot.title = element_text(size=15))+
  ylim(c(0,1))+
  labs(y="Correlation",title = "Correlation centroid-mean")

gcorr2=ggplot(ccm2,aes(x=Grupo,y=Cor,fill=factor(Grupo)))+
  geom_bar(stat = "identity",position = "dodge")+
  scale_fill_manual(values = cols)+
  geom_text(aes(label=Cor), position=position_dodge(width=0.9), vjust=-0.25)+theme_bw()+
  theme(axis.text = element_text(size=15,face = "bold"), legend.position = "none",
        axis.title = element_text(size=15),
        plot.title = element_text(size=15))+
  ylim(c(0,1))+
  labs(y="Correlation",x="Cluster",title = "Correlation: between members")

ggeary=ggplot(cgry,aes(x=Grupo,y=Cor,fill=factor(Grupo)))+
  geom_bar(stat = "identity",position = "dodge")+
  scale_fill_manual(values = cols)+
  geom_text(aes(label=Cor), position=position_dodge(width=0.9), vjust=-0.25)+theme_bw()+
  theme(axis.text = element_text(size=15,face = "bold"), legend.position = "none",
        axis.title = element_text(size=15),
        plot.title = element_text(size=15))+
  ylim(c(0,max(cgry$Cor)+0.002))+
  labs(y="Geary Index",x="Cluster",title = "Spatial correlation")

ggmoran=ggplot(cmrn,aes(Algo,Moran,fill=Algo))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(values="#EBF203")+
  annotate("text",x=1,y=-0.3,size=7,label="-0.64",fontface="bold")+theme_bw()+
  theme(axis.text = element_text(size=13,face = "bold"), legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_text(size=13,face = "bold"))+
  ylim(c(-0.7,0))+
  labs(x="Global Moran Index",y="")


library(patchwork)
lay=rbind(c(1,1,1,2,2),
          c(1,1,1,2,2),
          c(1,1,1,2,2))
gr=grid.arrange(spxg,arrangeGrob(gcorr,gcorr2,ggeary,ncol = 1),ncol=2)
gr=arrangeGrob(spxg,arrangeGrob(gcorr,gcorr2,ggeary,ncol = 1),ncol=2,layout_matrix = lay)
gr=as_ggplot(gr)
gr+inset_element(ggmoran, left = 0.3, bottom = 0.1, right = 0.5, top = 0.3)


gind1<-ggarrange(gcorr,gcorr2,ggeary,ncol = 1,nrow = 3)

gind<-ggarrange(gcorr,gcorr2,ncol = 1,nrow = 2)

gind<-annotate_figure(gind,top = text_grob("k-medias funcional",
                                           color = "black",
                                           face = "bold", size = 15))
ggarrange(gind,gind1,ncol = 2,nrow = 1)


## moran incluido ggplot
spxg1=spxg+inset_element(ggmoran, left = 0.6, bottom = 0.05, right = 0.9, top = 0.3)
save(spxg1,file="spxg1.Rdata")
load("spxg1.Rdata")


