#optime minimize, so i want minimize the - f(mi, k)


f1 <- function(par, stime, x,p){
  L=par
  k.rest=stime[p,1]
  m.rest=stime[p,2]
      fr<-function(L, k.rest, m.rest,x,p){
        m.rest- atan( (sum(L*x[,p]))/(sqrt( k.rest^2+  sum(L*x[,p])^2  )   )      )
        
      }
  fr(L,k.rest, m.rest,x,p)
}            



kprov <- optim( par=c(1,0.1) , fn=f1, stime=Stime, x=Phi2, p=3, method = "L-BFGS-B",hessian = T,lower =c( -Inf, -Inf), upper = c(Inf, Inf) )
 

kprov$convergence
kprov$par
                  