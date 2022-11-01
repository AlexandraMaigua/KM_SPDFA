# KM_SPDFA
K-means algorithm  for functional spatially correlated  data 

This algorithm was based on the `kmeans.fd` function of library [`fda.usc`](https://github.com/moviedo5/fda.usc) .

## Functions description

------------------------------------------------------------------------------------
### `functions.R`

Contains all major functions for the execution of the **K-means** algorithm.

+ `Ml2_dist(data, nbasis=65)`
  + Calculate the distance matrix between curves using the $L^2$ norm.
  + `data` data matrix with curves.
  + `nbasis` number of basis functions for curve smoothing (`data`)

+ `kmeans.center.iniL2(fdataobj,Vmdist=FALSE,coord=NULL, ncl = 2, 
                               cov.model="spherical", Kappa=NULL,multivgm=multivgm,
                               method = "sample", max.iter = 10000, max.comb = 1e+06, ...)`
  + Initializes the centroids.
  + `fdataobj` functional object containing the curves data.
  + `Vmdist` parameter that specifies that the calculation of the distance matrix is performed with mtv or tvm weighting. Default `FALSE`.
  +  `coord` field specifying that the coordinates to be used in the calculation of the distance matrix as long as `Vmdist=TRUE`.
  +  `cov.model` specifies the covariance model to be used for variogram modeling and matrix weighting. It works with common `vgm` and `cov.spatial` models.
  +  `Kappa` parameter that works with `cov.model="Mat"`.
  +  `multivgm` parameter indicating whether the weighting is performed with the multivariate variogram.
  +  `method` the method used to initialize the centroids. `sample` takes a sample of curves and calculates the `K` maximum distances to be the centroids. `exact` calculates all distances between curves and selects `K` centroids. The default `sample`.
  +  `max.iter` maximum number of iterations for the choice of centroids, works if `method="sample"` and serves as a stopping criterion.
  +  `max.comb` maximum number of site selection, works if `method=exact`.
  
  
+ `MVari(data,coord,nugget.fix=NULL, max.dist.variogram=NULL,
                multivgm=multivgm ,nbasis=65,cov.model="Gau",Kappa=NULL)`
  + Calculates the distance matrix between curves weighted by the plot-variogram or multivariogram.
  +  `data` matrix with curves data.
  +  `coord` matrix of coordinates associated with the `data` curves.
  +  `nugget.fix` specifies whether for the estimation of the nugget effect, the `nugget` is set to some value or it is estimated automatically.
  +  `max.dist.variogram` maximum distance used in the calculation of the empirical variogram and its estimation.
  +  `multivgm` specifies whether the weighting is done by multivariogram.
  +  `nbasis` number of basis functions for the smoothing of the curves and their distance calculation with `Ml2_dist`.
  +  `cov.model` spatial covariance model used in the estimation of the variogram.
  +  `Kappa` parameter used if `cov.model="Mat"`.


+ `kmeans.assig.groups(out,draw=FALSE,...)`
  + Assign the curves to the groups.
  + `out` object resulting from applying `kmeans.center.iniL2`.
  + `draw` specifies whether you want to display the assignment of the curves at each iteration.


+ `kmeans.centers.update(out,group
                               ,dfunc=func.trim.FM,draw=FALSE
                               ,par.dfunc=list(trim=0.05)
                               ,...)`
  + Updates the centroids.
  + `out` object containing a previous mapping of curves to centers as `kmeans.center.iniL2`.
  + `dfunc` depth measure, returns the average of the deepest (1-`trim`)% curves, by default the Fraiman-Muniz (`func.trim.FM`) is used.
  + `draw` specifies whether to display the centroid update for each iteration.
  + `par.dfunc` object that stores the parameters for the use of the depth measurement `dfunc`.

+ `kmeans.fdM(fdataobj, ncl = 2, Vmdist=FALSE,coord=NULL,
                      cov.model="spherical", Kappa=NULL, multivgm=multivgm,
                      dfunc = func.trim.FM, 
                      max.iter = 10000, par.dfunc = list(trim = 0.05), 
                      method = "sample", cluster.size = 1, draw = FALSE, ...) `
                      
  + Performs the full process of the k-means algorithm for spatially correlated functional data.
  + `fdataobj` object of type `fdata` containing the curves.
  + `ncl` number of clusters to be formed.
  + `Vmdist` parameter that specifies that the calculation of the distance matrix is performed with mtv or tvm weighting. Default `FALSE`.
  +  `coord` field specifying that the coordinates to be used in the calculation of the distance matrix as long as `Vmdist=TRUE`.
  +  `cov.model` specifies the covariance model to be used for variogram modeling and matrix weighting. It works with common `vgm` and `cov.spatial` models.
  +  `Kappa` parameter that works with `cov.model="Mat"`.
  +  `multivgm` parameter indicating whether the weighting is performed with the multivariate variogram.
  +  `dfunc` depth measure, returns the average of the deepest (1-`trim`)% curves, by default the Fraiman-Muniz (`func.trim.FM`) is used.
  +  `max.iter` maximum number of iterations for the choice of centroids, works if `method="sample"` and serves as a stopping criterion.
  + `par.dfunc` object that stores the parameters for the use of the depth measurement `dfunc`.
  +  `method` the method used to initialize the centroids. `sample` takes a sample of curves and calculates the `K` maximum distances to be the centroids. `exact` calculates all distances between curves and selects `K` centroids. The default `sample`.
  +  `cluster.size` minimum number of elements in the clusters.
  + `draw` specifies whether to display the centroid update for each iteration.
   
   
+ `kmeans.fdas(fdataobj, ncl = 2,Vmdist=FALSE,coord=NULL,
                      cov.model=cov.model,Kappa=NULL,
                      multivgm=FALSE)`
  + It allows clustering using spatial information or only the information of the curves.
  + The function parameters are the same as the `kmeans.fdM` function.



+ `SSBindex1(a)`
  + Calculates the SSB index, between groups. Works only with centroid information.
  + `a` object resulting from applying the `kmeans.fdas` function.

+ `SSWindex<-function(a,data)`
   + Calculates the SSW index, within clusters. Works with the centroids and curves associated with each cluster.
   + `a` object resulting from applying the `kmeans.fdas` function.
   + `data` matrix with curves.

+ `metod(datat,metodo="tvm")`
  + Calculates the percentage of correct classification.
  + `datat` object resulting from the simulation process (`simu9.R`).
  + `metodo` specifies the results obtained from the `mtv` or `tvm` method.

+ `Datsim(Eu.d,cov.model="exponential",cov.pars=c(0.5,1.5),Kappa=NULL,
                media=c(5,6))`
  + Simulates curves with spatial correlation.
  + `Eu.d` spatial distance matrix.
  +  `cov.model` spatial covariance model for spatial dependence.
  +  `cov.pars` parameters of the covariance model.
  +  `Kappa` parameter if `cov.model="Mat"` is used.
  +  `media` vector of constant functional mean values.

+ `Datsim1(Eu.d=NULL,cov.model="exponential",cov.pars=c(0.5,1.5),Kappa=NULL,
                media=c(5,15))`
  + Simulates curves without spatial correlation.
  
 -------------------------------------------------------------------------------
### `git_moranindex.R`

Contains function for the calculation of the moran index for inter-group and intra-group.

+ `MoranGfda(b,weig.mat,data,coord,tipo="entre")`
  + `b` object resulting from the `kmean.fdas` algorithm.
  + `weig.mat` matrix of weights (check theoretical index). It is possible to use either the ayacency matrix of the groups or the matrix of distances between the curves.
  + `data` matrix of curves.
  + `coord` coordinate matrix.
  + `tipo` parameter that specifies whether the index calculation is performed between groups or within groups, by default `type="between"`.

-------------------------------------------------------------------------------
### `git_gearyindex.R`

Contains function for the calculation of the geary index for inter-group and intra-group.

+ `GearyGfda(b,weig.mat,data,coord,tipo="entre")`
  + `b` object resulting from the `kmean.fdas` algorithm.
  + `weig.mat` matrix of weights (check theoretical index). It is possible to use either the ayacency matrix of the groups or the matrix of distances between the curves.
  + `data` matrix of curves.
  + `coord` coordinate matrix.
  + `tipo` parameter that specifies whether the index calculation is performed between groups or within groups, by default `type="between"`. 
  
-------------------------------------------------------------------------------
### `aplderiv.R`

Script that contains the application of the previous method (`kmean.fdas`) to the NDVI data of Ecuador and also graphs the obtained groups.

-------------------------------------------------------------------------------
### `bout.Rdata`

Object containing the final classification of the NDVI data.

-------------------------------------------------------------------------------
### `covariances.R` (review)

Script that generates the covariance matrix through the coefficients of the basis functions. It uses parameters of the linear co-regionalization model.

-------------------------------------------------------------------------------
### `fgeary.R`

Script containing the function for the calculation of spatial correlation using the Geary index. The function is designed to be used in parallel processes. 

-------------------------------------------------------------------------------
### `fit.lmc.R` (review)

Estimates a linear model of co-regionalization to the coefficients of the function basis.

-------------------------------------------------------------------------------
### `fmoran.R`

Script containing the function for the calculation of spatial correlation using the Moran index. The function is designed to be used in parallel processes. 


-------------------------------------------------------------------------------
### `git_ghsp.R`

Script que se utiliza para la generación de gráficos de la clasificación de grupos sobre shapes de demarcación hidrográfica.

-------------------------------------------------------------------------------
### `git_pixeles.R`

Script used to create group plots obtained by pairs, on hydrographic demarcation shapes.

-------------------------------------------------------------------------------
### `git_simu4.R`

Script used to perform specific simulations and obtain the classification graphically.

-------------------------------------------------------------------------------
### `git_simu9.R`

Script used to perform the simulations in parallel and obtain the percentages of correct classification.

-------------------------------------------------------------------------------
### `git_test_function.R`

Script for the realization of graphs of curves in each group and spatial and temporal correlations.

-------------------------------------------------------------------------------
### `multiv.R` (**review theory**)

Script with the multivariogram estimation function, using LMC principles and a basis of otonormal functions.

-------------------------------------------------------------------------------
### `parameters.lmc.R` (review)

Generates the matrix with parameters of the linear co-regionalization model fitted to the quotients of the basis functions.

-------------------------------------------------------------------------------
### `variogrmas.R` (revisar)

Function that generates the variogram matrices through the coefficients of the basis functions using parameters of the linear co-regionalization model.

-------------------------------------------------------------------------------
### `git_moranindex.R`

Contains function for the calculation of the moran index for inter-group and intra-group.

+ `MoranGfda(b,weig.mat,data,coord,tipo="entre")`
  + `b` object resulting from the `kmean.fdas` algorithm.
  + `weig.mat` matrix of weights (check theoretical index). It is possible to use either the ayacency matrix of the groups or the matrix of distances between the curves.
  + `data` matrix of curves.
  + `coord` coordinate matrix.
  + `tipo` parameter that specifies whether the index calculation is performed between groups or within groups, by default `type="between"`.

-------------------------------------------------------------------------------
### `git_gearyindex.R`

Contains function for the calculation of the geary index for inter-group and intra-group.

+ `GearyGfda(b,weig.mat,data,coord,tipo="entre")`
  + `b` object resulting from the `kmean.fdas` algorithm.
  + `weig.mat` matrix of weights (check theoretical index). It is possible to use either the ayacency matrix of the groups or the matrix of distances between the curves.
  + `data` matrix of curves.
  + `coord` coordinate matrix.
  + `tipo` parameter that specifies whether the index calculation is performed between groups or within groups, by default `type="between"`.





