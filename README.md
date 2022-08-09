# KM_SPDFA
K-means algorithm  for functional spatially correlated  data 

## Instructivo

### `functions.R`

Contiene todas la mayotia de funciones para la ejecucuíon del algoritmo **K-medias**

+ `Ml2_dist(data, nbasis=65)`
  + Calcula la matriz de distancia entre curvas utilizando la norma $L^2$.
  + `data` matriz de datos con curvas.
  + `nbasis` número de funciones base para la suavización de las curvas (`data`)

+ `kmeans.center.iniL2(fdataobj,Vmdist=FALSE,coord=NULL, ncl = 2, draw = FALSE,
                               cov.model="spherical", Kappa=NULL,multivgm=multivgm,
                               method = "exact", max.iter = 10000, max.comb = 1e+06, 
                               par.metric = NULL, ...)`
                               
  + `fdataobj` objeto funcional que contiene los datos de las curvas.
  + `Vmdist` parámetro que especifíca que el cálculo de la matriz de distancia se realice con ponderación mtv o tvm. Por defecto `FALSE`.
  +  `coord` campo que especifíca que las coordenadas a utilizarse en el cálculo de la matriz de distancia siempre y cuando `Vmdist=TRUE`.





