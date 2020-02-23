#general function for phi_3 , ... , phi_128
hat.k=c(k1$par, k2.L12$par[1])

hat.L=c(0,k2.L12$par[2])

#function that give part of the log-likehood 
M <- function(hat.K, hat.L, X ,j ){
  A=0
  A=sum( - log(besselI( hat.K[1] ,0))  + hat.K[1]*cos(X[,1]) ) 
  
  
  
  
  
  for (i in 2:(j-1)) {
    B=0
    for (l in 1:(j-1)) {
      B=B+sum( hat.L[l]*sin(X[,i])  )
    }
    
    A=A+ sum( - log(  besselI(  sqrt(hat.K[i]^2+B^2) ,0))     + sqrt(hat.K[i]^2+B^2)*cos( X[,i] - atan( B/ hat.K[i] )  ))
   
  }
  return(A)
}




frng <- function(parms, x, j, memo){ 

  

   -(   -50*j*log(2*pi)+memo+sum( -log(besselI( sqrt( parms[1]^2 +(sum(parms[2:j]%*%t( rep(x[,j],j-2 ))  ))^2 )  ,0)) +  sqrt( parms[1]^2 + (sum(parms[2:j]%*%t( rep(x[,j],j-2 ))  ))^2 )*cos(x[,j] - atan( sum(parms[2:j]%*%t( rep(x[,j],j-2 )) )/parms[1]) )  )   
   )
  
   
    
}


fr<- function(parms, x, j, memo){ 
  
  
  
  -(   -50*j*log(2*pi)+memo+sum(  -log(besselI( parms[1] ,0)) + parms[1]*cos(x[,j] - parms[2])   )  )
  
  
  
}

kj.Lj <- optim(c(1,0), fr,  x=Phi2, j=3, memo=g, #method = "BFGS")
               method="L-BFGS-B", hessian = 1,
               lower = c(0, rep(-Inf, 2)), 
               upper = c(Inf, rep(Inf, 2)))

#check convergence
kj.Lj$convergence

#estimates
kj.Lj$par

#RSS
diag(solve(kj.Lj$hessian))


#j is refer to parameters, j \in [3,128]
#aggiorno le stime
#hat.k=cbind(hat.k, kj.Lj$par[1])
#hat.L=c(hat.L, kj.Lj$par[-1])


#calcolo la likehood aggiornata
#g=M(hat.K, hat.L, Phi2 ,j )

j=3
g=M(hat.k, hat.L, Phi2 ,j=3 )
g
kj.Lj <- optim(c(1,0,0), frng,  x=Phi2, j=3, memo=g, #method = "BFGS")
                method="L-BFGS-B", hessian = 1,
                lower = c(0, rep(-Inf, j-1)), 
                upper = c(Inf, rep(Inf, j-1)))

#check convergence
kj.Lj$convergence

#estimates
kj.Lj$par

#RSS
diag(solve(kj.Lj$hessian))
