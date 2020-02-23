#stima della prima f(phi_1)####

frn <- function(parms, X){
  k=parms
  n=length(X)
  x=X
  
  f<-function(k=k, n=n, x=x){
    -n*log(2*pi)+ sum(-log( besselI(k,0) )+ k*cos(x))
    
  }
  #perchè voglio il min della - log-like
  -f(k,n,x)
  
}

k1 <- optim(99, frn,  X=Phi2[,1], method = "L-BFGS-B", hessian = T,lower = 0, upper = 100 )

#stima k1
k1$par

#stima SE(k1)
k1$hessian^(-1)


##stima  f(phi_2 | phi_1)####


frn1 <- function(parms, X){ 
  k2=parms[1]
  l=parms[2]
  x=X[,1]
  y=X[,2]
 # a=sqrt(k2^2+sum(l12^2*sin(x))^2)
  f<-function(parms,x,y){
  -50*log(2*pi)+sum( - log(besselI(  sqrt( parms[1]^2+parms[2]^2*sin(x)^2 ) ,0) ) + parms[1]*cos(y)+parms[2]*sin(x)*sin(y)  )
  }
  #perchè voglio il min della - log-like
  -f(c(k2,l),x,y)
  
}


k2 <- optim(c(4,0), frn1,  X=Phi2,
  method="L-BFGS-B", hessian = 1,
  lower = c(0, -Inf), 
  upper = c(Inf, Inf))

#check convergence
k2$convergence

#estimates
k2$par

#RSS
diag(solve(k2$hessian))
