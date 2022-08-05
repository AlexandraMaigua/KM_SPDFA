###########################################################################################
# Function for generating variogram matrices among coefficients of Fourier basis functions
# This function uses parameters of the linear model o coregionalization 
###########################################################################################


mvario<-function(para,d)
{
       s <- dim(d)[1]
       k <- dim(para)[1]
       vario<-array(0,dim=c(s,s,k,k))

       for (l in 1:s){
             for(i in l:s)
                   {
                    vario[l,i,,] <- para[,,1]+para[,,2]*(1-exp(-(d[i,l]/para[,,3])))
                    vario[i,l,,] <- vario[l,i,,]
                    }
             }
      return(vario)
}

