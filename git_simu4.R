source("functions.R")

set.seed(123)

c1=c(runif(30,-2,0) ,runif(30,-2,0), runif(30,0,2), runif(30,0,2))
c2=c(runif(30,-2,0) ,runif(30,0,2), runif(30,-2,0), runif(30,0,2))

coord=data.frame(cbind(c1,c2))
names(coord)=c("x","y")

Eu.d <-as.matrix(dist(coord,method="euclidian"))

### EXPONENCIAL

set.seed(999)
data=Datsim(Eu.d,cov.model = "Gau",cov.pars = c(2,1.8),
            Kappa =NULL,media = c(5,15))

# grafico de curvas simuladas

matplot((data),type="l",xlab = "Day",
        ylab = "",font.lab=2, cex.lab = 2.5 , cex.axis=2.3, cex.main =2.3,
        col = c(rep("#33FF9F",30),rep("#FF8D33",30),rep("#FF8D33",30),rep("#33FF9F",30)))

coord=data.frame(cbind(c1,c2))
names(coord)=c("x","y")

set.seed(123)

a<-kmeans.fdas(fdataobj = t(data),ncl = 4,Vmdist = TRUE,coord = coord,cov.model = "Gau",
               Kappa =NULL,multivgm = FALSE)
table(a$cluster)


grupo=data.frame(grupo=a$cluster)
coord=data.frame(coord)
coord$Grupo=as.factor(grupo$grupo)

ggplot(coord,aes(x=x,y=y,color=Grupo))+
  geom_point(size=2.3)+
  labs(title = "Modelo gaussiano",subtitle = "u1=5, u2=15, parámetros=(4,2.5)")+
  theme_bw()+
  theme_bw(base_size = 20)+
  theme(legend.position="none")

