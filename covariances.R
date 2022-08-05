
###########################################################################################
# Function for generating covariance matrices among coefficients of Fourier basis functions
# This function uses parameters of the linear model o coregionalization 
###########################################################################################


mcov<-function(para,d)
{
       s <- dim(d)[1]
       k <- dim(para)[1]
       covari<-array(0,dim=c(s,s,k,k))
 
       for (l in 1:s){
             for(i in l:s)
                   {
                    covari[l,i,,] <- para[,,2]*exp(-(d[i,l]/para[,,3]))
                    covari[i,l,,] <- covari[l,i,,]
                    }
             }
      return(covari)
}


