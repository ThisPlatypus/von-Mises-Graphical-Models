count=0
Stime=matrix(0, ncol = 4, nrow = 128)

Stime[1:2,1]=c(k1$par, k2.L12$par[1])
Stime[1:2,2]=c(0,k2.L12$par[2])
for (p in 3:128) {
  
fr<- function(parms, x, j, memo){ 
  
  
  
  -(   -50*j*log(2*pi)+memo+sum(  -log(besselI( parms[1] ,0)) + parms[1]*cos(x[,j] - parms[2])   )  )
  
  
  
}
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

g=M(Stime[1:p,1], Stime[1:p,2], Phi2 ,j=p )
kj.Lj <- optim(c(0.5,-1), fr,  x=Phi2, j=p, memo=g, #method = "BFGS")
               method="L-BFGS-B", hessian = 1,
               lower = c(0, rep(-Inf, 2)), 
               upper = c(Inf, rep(Inf, 2)))

#check convergence
count=count+kj.Lj$convergence

#estimates
Stime[p,1]=kj.Lj$par[1]
Stime[p,2]=kj.Lj$par[2]
Stime[p, 3:4]=diag(solve(kj.Lj$hessian))
}

count
Stime
round(stime1[,1],3)==round(Stime[,1],3)
round(stime1[,2],3)==round(Stime[,2],3)
stime1=Stime


