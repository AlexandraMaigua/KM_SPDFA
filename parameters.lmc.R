###################################################################################################
#  This function generate a matrix with parameters of the linear model of corregionalization
#  fitted to coeffcients of Fourier basis functions
###################################################################################################



parlmc<-function(k,g.fit)
   {
   
   paralmc<-array(0,dim=c(k,k,3))

   m1 <- matrix(0,k,k)
   m1[lower.tri(m1,diag=TRUE)] <- 1:(k*(k+1)/2)
   m2 <- t(m1)
   m1[upper.tri(m1)]<-m2[upper.tri(m2)]

   secuencia<-matrix(m1,nrow=k*k)
   for (j1 in 1:k){
      for (j2 in 1:k)
        {
         j <- (j1-1)*k + j2
         i <- secuencia[j]
         paralmc[j1,j2,] <-c(g.fit$model[[i]][,2],g.fit$model[[i]][2,3])
        }  
     }
   return(paralmc)
}



