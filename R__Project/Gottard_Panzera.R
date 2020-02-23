library(circular)
param=c(rep(1,3))
j=3
n=length(Phi2[,1])
Phi2=as.matrix(Phi2)
Phi2=circular(Phi2, units = "radians")

loglike<- function(par, Phi2, j, quali =1:(j-1) ){
  kj=par[1]
  lambda=par[-1]
  cmu=atan((kj)^(-1)*(lambda%*%t(sin(Phi2[,quali]) )))
ckappa=sqrt( kj^2+ (  lambda%*%t(sin(Phi2[,quali]) ) )^2  )
cmu=apply(cmu, 1, circular)
llike=0
for (i  in 1:n) {
  llike=llike+dvonmises(Phi2[i,j], cmu[i], ckappa[i], log=T)
}
return(-llike)
}


Stime=NULL
#da far girare

for (p in 3:9) {
  k.l<- optim(c(rep(0.2,p)), loglike, Phi2=Phi2, j=p, hessian = T, lower = c(0, rep(-Inf, p-1)), upper = c(rep(Inf,p)))
  Stime[[p]] <- list( StimePar=k.l$par, DiagHess=diag(solve(k.l$hessian) ), logL=k.l$value )
  
}

save(Stime, file="fname.RData")
#va bene fino alla 51
#load(fname.RData)

