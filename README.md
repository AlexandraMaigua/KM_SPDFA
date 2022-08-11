# KM_SPDFA
K-means algorithm  for functional spatially correlated  data 

Este algoritmo se basó en la función `kmeans.fd` de la librería [`fda.usc`](https://github.com/moviedo5/fda.usc) .

## Instructivo

------------------------------------------------------------------------------------
### `functions.R`

Contiene todas la mayoria de funciones para la ejecucuíon del algoritmo **K-medias**

+ `Ml2_dist(data, nbasis=65)`
  + Calcula la matriz de distancia entre curvas utilizando la norma $L^2$.
  + `data` matriz de datos con curvas.
  + `nbasis` número de funciones base para la suavización de las curvas (`data`)

+ `kmeans.center.iniL2(fdataobj,Vmdist=FALSE,coord=NULL, ncl = 2, 
                               cov.model="spherical", Kappa=NULL,multivgm=multivgm,
                               method = "sample", max.iter = 10000, max.comb = 1e+06, ...)`
  + Inicializa los centroides.                        
  + `fdataobj` objeto funcional que contiene los datos de las curvas.
  + `Vmdist` parámetro que especifíca que el cálculo de la matriz de distancia se realice con ponderación mtv o tvm. Por defecto `FALSE`.
  +  `coord` campo que especifíca que las coordenadas a utilizarse en el cálculo de la matriz de distancia siempre y cuando `Vmdist=TRUE`.
  +  `cov.model` especifica el modelo de covarianza a utilizarse para la modelización del variograma y la ponderación de la matriz. Trabaja con modelos comunes de `vgm` y `cov.spatial`.
  +  `Kappa` parámetro que trabaja con `cov.model="Mat"`.
  +  `multivgm` parámetro que indica si la ponderación se realiza con el variograma multivariado.
  +  `method` el método utilizado para inicializar los centroides. `sample` toma una muestra de curvas y calcula las `K` máximas distancias para que sean los centroides. `exact` calcula todas las distancias entre las curvas y selecciona `K` centroides. Por defecto `sample`.
  +  `max.iter` número máximo de iteraciones para la elección de los centroides, funciona si `method="sample"` y sirve como criterio de parada.
  +  `max.comb` número máximo de selección de centros, funciona si `method=exact` .
  
  
+ `MVari(data,coord,nugget.fix=NULL, max.dist.variogram=NULL,
                multivgm=multivgm ,nbasis=65,cov.model="Gau",Kappa=NULL)`
  + Calcula la matriz de distancia entre curvas ponderada por el trazo-variograma o multivariograma.
  +  `data` matriz con datos de las curvas.
  +  `coord` matriz de coordenadas asociadas a las curvas de `data`.
  +  `nugget.fix` especifíca si para la estimación del efecto pepita, el `nugget` se fija a algún valor o lo estima automático.
  +  `max.dist.variogram` distancia máxima utilizada en el cálculo del variograma empírico y la estimación del mismo.
  +  `multivgm` especifíca si la ponderación se la realiza mediante multivariograma.
  +  `nbasis` número de funciones base para la suavización de las curvas y su cálculo de distancias con `Ml2_dist`.
  +  `cov.model` modelo de covarianza espacial utilizada en la estimación del variograma.
  +  `Kappa` parámetro utilizado si `cov.model="Mat"`.


+ `kmeans.assig.groups(out,draw=FALSE,...)`
  + Asigna las curvas a los grupos.
  + `out` objeto resultante de aplicar `kmeans.center.iniL2`
  + `draw` especifíca si se quiere visualizar la asignación de las curvas en cada iteración.


+ `kmeans.centers.update(out,group
                               ,dfunc=func.trim.FM,draw=FALSE
                               ,par.dfunc=list(trim=0.05)
                               ,...)`
  + Actualiza los centroides.
  + `out` objeto que contiene una asignación de curvas a centros previa como `kmeans.center.iniL2`.
  + `dfunc` medida de profundiad, retorna el promedio de las (1-`trim`)%  curvas mas profundas, por defecto se utiliza la de Fraiman-Muniz (`func.trim.FM`).
  + `draw` especifíca si se quiere visualizar la actualización de los centroides de cada iteración.
  + `par.dfunc` objeto que almacena los parámetros para el uso de la medidad de profundidad `dfunc`.

+ `kmeans.fdM(fdataobj, ncl = 2, Vmdist=FALSE,coord=NULL,
                      cov.model="spherical", Kappa=NULL, multivgm=multivgm,
                      dfunc = func.trim.FM, 
                      max.iter = 10000, par.dfunc = list(trim = 0.05), 
                      method = "sample", cluster.size = 1, draw = FALSE, ...) `
                      
  + Realiza el proceso completo del algoritmo k-medias para datos funcionales con correlación espacial.
  + `fdataobj` objeto de tipo `fdata` que contiene las curvas.
  + `ncl` número de clústers a formarse.
  + `Vmdist` parámetro que especifíca que el cálculo de la matriz de distancia se realice con ponderación mtv o tvm. Por defecto `FALSE`.
  +  `coord` campo que especifíca que las coordenadas a utilizarse en el cálculo de la matriz de distancia siempre y cuando `Vmdist=TRUE`.
  +  `cov.model` especifica el modelo de covarianza a utilizarse para la modelización del variograma y la ponderación de la matriz. Trabaja con modelos comunes de `vgm` y `cov.spatial`.
  +  `Kappa` parámetro que trabaja con `cov.model="Mat"`.
  +  `multivgm` parámetro que indica si la ponderación se realiza con el variograma multivariado.
  +  `dfunc` medida de profundiad, retorna el promedio de las (1-`trim`)%  curvas mas profundas, por defecto se utiliza la de Fraiman-Muniz (`func.trim.FM`).
  +  `max.iter` número máximo de iteraciones para la elección de los centroides, funciona si `method="sample"` y sirve como criterio de parada.
  + `par.dfunc` objeto que almacena los parámetros para el uso de la medidad de profundidad `dfunc`.
  +  `method` el método utilizado para inicializar los centroides. `sample` toma una muestra de curvas y calcula las `K` máximas distancias para que sean los centroides. `exact` calcula todas las distancias entre las curvas y selecciona `K` centroides. Por defecto `sample`.
  +  `cluster.size` número de elementos mínimo en los clústers.
  + `draw` especifíca si se quiere visualizar la actualización de los centroides de cada iteración.
   
   
+ `kmeans.fdas(fdataobj, ncl = 2,Vmdist=FALSE,coord=NULL,
                      cov.model=cov.model,Kappa=NULL,
                      multivgm=FALSE)`
  + Permite realizar la clusterización utilizando la información espacial o solo la de las curvas.
  + Los parámetros de la función son los mismos que la función `kmeans.fdM`.



+ `SSBindex1(a)`
  + Calcula el índice SSB, entre grupos. Trabaja solo con la información de los centroides.
  + `a` objeto resultante de aplicar la función `kmeans.fdas`.

+ `SSWindex<-function(a,data)`
   + Calcula el índice SSW, intra grupos. Trabaja con los centroides y las curvas asociadas a cada cluster.
   + `a` objeto resultante de aplicar la función `kmeans.fdas`.
   + `data` matriz con curvas.

+ `metod(datat,metodo="tvm")`
  + Calcula el porcentaje de correcta clasificación.
  + `datat` objeto resultante del proceso de simulación (`simu9.R`).
  + `metodo` especifica los resultados obtenidos del método `mtv` o `tvm`.

+ `Datsim(Eu.d,cov.model="exponential",cov.pars=c(0.5,1.5),Kappa=NULL,
                media=c(5,6))`
  + Simula curvas con correlación espacial.
  + `Eu.d` matriz de distancias espacial.
  +  `cov.model` modelo de covarianza espacial para la dependencia espacial.
  +  `cov.pars` parámetros del modelo de covarianza.
  +  `Kappa` parámetro si se utiliza `cov.model="Mat"`.
  +  `media` vector de valores de las medias funcionales constantes.

+ `Datsim1(Eu.d=NULL,cov.model="exponential",cov.pars=c(0.5,1.5),Kappa=NULL,
                media=c(5,15))`
  + Simula curvas sin correlación espacial.
  
-------------------------------------------------------------------------------
### `aplderiv.R`

Script que contiene la aplicación del método anterior (`kmean.fdas`) a los datos del NDVI del Ecuador y también realiza el gráfico de los grupos obtenidos.

-------------------------------------------------------------------------------
### `bout.Rdata`

Objeto que contiene la clasificación final de los datos del NDVI.

-------------------------------------------------------------------------------
### `covariances.R` (revisar)

Script que genera la matriz de covarianza a través de los coeficientes de la funciones base. Usa parámetros del modelo lineal de corregionalización.

-------------------------------------------------------------------------------
### `fgeary.R`

Script que contine la función para el cálculo de correlación espacial utilizando el índice de Geary. La función está diseñada para que se utilice en procesos de paralelos. 

-------------------------------------------------------------------------------
### `fit.lmc.R` (revisar)

Estima un modelo lineal de corregionalización a los coeficientes de la base de funciones.

-------------------------------------------------------------------------------
### `fmoran.R`

Script que contine la función para el cálculo de correlación espacial utilizando el índice de Moran. La función está diseñada para que se utilice en procesos de paralelos. 


-------------------------------------------------------------------------------
### `git_ghsp.R`

Script que se utiliza para la generación de gráficos de la clasificación de grupos sobre shapes de demarcación hidrográfica.

-------------------------------------------------------------------------------
### `git_pixeles.R`

Script utilizado para crear gráficos de grupos obtenidos por pares, sobre shapes de demarcación hidrográfica.

-------------------------------------------------------------------------------
### `git_simu4.R`

Script utilizado para realizar simulaciones específicas y obtener la clasificación de manera gráfica.

-------------------------------------------------------------------------------
### `git_simu9.R`

Script utilizado para realizar las simulaciones en paralelo y obtener los porcentajes de correcta clasificación.

-------------------------------------------------------------------------------
### `git_test_function.R`

Script para la realización de gráfcios de curvas en cada grupo y correlaciones espaciales y temporales.

-------------------------------------------------------------------------------
### `multiv.R` (**revisar teoría**)

Script con la función de estimación del multivariograma, utilizando principios del LMC y una base de funciones otonormales.

-------------------------------------------------------------------------------
### `parameters.lmc.R` (revisar)

Genera la matriz con parámetros del modelo lineal de corregionalización ajustado a las coefcientes de las funciones base.

-------------------------------------------------------------------------------
### `variogrmas.R` (revisar)

Función que genera las matrices de variogramas a través de los coeficientes de las funciones base utilizando parámetros del modelo lineal de corregionalización.






