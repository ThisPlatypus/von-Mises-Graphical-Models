#general function for phi_3 , ... , phi_128
hat.k=c(k1$par, k2.L12$par[1])

hat.L=c(0,k2.L12$par[2])


d=3.265555
for (p in 3:128) {
  

#function that give part of the log-likehood 
#M <- function(hat.K, hat.L, X ,j ){
  #A=0
  #per il primo pezzo p=1
  #A=sum( - log(besselI( hat.K[1] ,0))  + hat.K[1]*cos(X[,1]) ) 
  
   # for (i in 2:(j-1)) {
  #  B=0
   # for (l in 1:(j-1)) {
    #  B=B+sum( hat.L[l]*sin(X[,i])  )
    #}
    
    #A=A+ sum( - log(  besselI(  sqrt(hat.K[i]^2+B^2) ,0))     + sqrt(hat.K[i]^2+B^2)*cos( X[,i] - atan( B/ hat.K[i] )  ))
   
  #}
  #return(A) }


#pezzo nuovo
regress= function(k,l,X,j,g ){
  g+sum( -log( besselI( sqrt(k^2+sum((l)%*%t(sin(X[,j])) )^2) ,0)   ) + sqrt(k^2+sum((l)%*%t(sin(X[,p]))) ^2) *cos(X[,j]- atan(  (sum((l)%*%t(sin(X[,p]))) ) /k))  )
}



#frng <- function(parms, X, j, memo){ 

 # make_B <- function(parms, x,j){
   # B=0
    #for (l in 2:j) {
  #    B=B+sum( parms[l]*sin(x[,j])  )
   # }
    #B
    
  #}

  f<-function(parms,x,j, memo){
    
  # B=make_B(parms,x,j)
    
   -( -50*j*log(2*pi)+memo+sum( -log(besselI( sqrt( parms[1]^2 +(sum( parms[-1]%*%t(sin(x[,j])) ) )^2 )  ,0)) +  sqrt( parms[1]^2 + sum( parms[-1]%*%t(sin(x[,j]) ) )^2 )*cos(x[,j] - atan( sum( parms[-1]%*%t(sin(x[,j]) ) )/parms[1]) )  )   
   )
  
    }
  #perchè voglio il min della - log-like
#  -f(parms,X, j, memo)

  
    
#}


#j is refer to parameters, j \in [3,128]
#aggiorno le stime
#hat.k=cbind(hat.k, kj.Lj$par[1])
#hat.L=c(hat.L, kj.Lj$par[-1])


#calcolo la likehood aggiornata
#g=M(hat.K, hat.L, Phi2 ,j )

#j=3
#g=M(hat.k, hat.L, Phi2 ,j=p )

kj.Lj <- optim(par=c(0.4,0.5,-0.3), f,  x=Phi2, j=p, memo=d, #method = "BFGS")
                method =  "L-BFGS-B", hessian=T,
                lower = c(0.00000000001, rep(-Inf, p-1)), 
                upper = c(Inf, rep(Inf, p-1)))
d=regress( kj.Lj$par[1], kj.Lj$par[-1], Phi2,p, d)

#check convergence
count=count+kj.Lj$convergence

#estimates
stime[p,1]=kj.Lj$par[1]
stime[p,2:length(kj.Lj$par)-1]=kj.Lj$par[-1]

}





kj.Lj
#check convergence
kj.Lj$convergence

#estimates
kj.Lj$par

#RSS
diag(solve(kj.Lj$hessian))
