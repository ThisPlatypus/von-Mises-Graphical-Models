#Mi####
make_mi=function(data){
  mi=circular(0)
  for (i in 1:dim(data)[2]) {
    mi[i]=mean.circular(data[,i], units="degrees")
  }
  return(mi)
}


#Kappa####
make_k=function(data){
  ko=circular(0, units="degrees")
  for (i in 1:dim(data)[2]) {
    ko[i]=var.circular(data[,i])
  }
return(ko)
  
}



#Starter value####
medie=make_mi(Phi1)

kappino=var(Phi1)

#ho deciso arbitrariamente di metterci la matrice di varianza-covarianza
lappo=cov(Phi1)
diag(lappo)=0
