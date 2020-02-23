
###Tutte le derivate sono rispetto alla log L, quindi vanno cambiate di segno nell'aggiornamento


#Per debug
theta=Phi1
k=kappino
mu=medie
L=lappo

  


trust=function(theta, mu, k, L){
  

#Starter
kStore=mStore=matrix(NA, dim(theta)[1], dim(theta)[2] )

for (i in 1:dim(theta)[1]) {
  
  kStore[i,2]=sqrt(k[2]^2+ (L[1,2]*sin(c(theta[i,2]-mu[2])))^2)
  mStore[i,2]=mu[2]+ atan(L[1,2]*sin(c(theta[i,2]-mu[2]))/k[2])
  
  for (j in 3:dim(theta)[2]) {
      #somma parziale del k 
#deve scorrere j e i 
Sk=sum(L[j,1:(j-1)]*rapply(c(theta[i,1:(j-1)]-mu[1:(j-1)]),  sin))^2
#kappa a meno di j nel passo i
kStore[i,j]=sqrt(k[j]^2+Sk)

#somma parziale di mu
Sm=sum( L[j,1:(j-1)]*rapply(c(theta[i, 1:(j-1)]-mu[1:(j-1)]), sin) )
#mu a meno di j nel passo i 
mStore[i,j]=mu[j]+ atan(Sm/k[j])


 }

  
}


#####LAMBDA


part_L=matrix(0, dim(theta)[2], dim(theta)[2])

for (s in 1:(dim(theta)[2])) {
  for (i in 1:(dim(theta)[1])) {
  ##QUESTO PERCHè RAPPLY() NON VA BENE SE IL VETTORE HA DIMENSIONE 1
  #pezzi della derivata della likelihood rispetto ai lambda
  D_k_ij_D_L_rs=sum(L[2,1]*sin(c(theta[i, 2]-mu[2]))*sin(c(theta[i, s]-mu[s])))/kStore[i,2]
  
  D_m_ij_D_L_rs=(sin(c(theta[i,s]-mu[s]))/ (k[2]*(1+  (  sum(L[j, 1:2]*sin(c(theta[i, 2]-mu[2])))/ kStore[i,2]  )^2 )  ))
  
  #Derivata parziale intera
  #Riempio solo la diagonale superiore, quindi j(2-(p-1)) ed s (1-p)
  part_L[1,2]= D_k_ij_D_L_rs * ( -A1FirstDerivative(kStore[i,2])+ cos( c(theta[i,2]- mStore[i,2]) )  )+
    D_m_ij_D_L_rs * kStore[i,2]*sin(  c(theta[i,j]- mStore[i,2]) )   
  }
  
for (j in 3:dim(theta)[2]) {
  for (i in 1:dim(theta)[1]) {
    if(s!=j & s<j ){
     #pezzi della derivata della likelihood rispetto ai lambda
 D_k_ij_D_L_rs=sum(L[j,1:(j-1)]*rapply(c(theta[i, 1:(j-1)]-mu[1:(j-1)]), sin)*sin(c(theta[i, s]-mu[s])))/kStore[i,j]

D_m_ij_D_L_rs=(sin(c(theta[i,s]-mu[s]))/ (k[j]*(1+  (  sum(L[j, 1:(j-1)]*rapply(c(theta[i, 1:(j-1)]-mu[1:(j-1)]), sin))/ kStore[i,j]  )^2 )  ))

#Derivata parziale intera
#Riempio solo la diagonale superiore, quindi j(2-(p-1)) ed s (1-p)
part_L[s,j]= D_k_ij_D_L_rs * ( -A1FirstDerivative(kStore[i,j])+ cos( c(theta[i,j]- mStore[i,j]) )  )+
                            D_m_ij_D_L_rs * kStore[i,j]*sin(  c(theta[i,j]- mStore[i,j]) )   
 
    }  
    }

}
}


########KAPPA al meno di j
##vengono negtivi
#Perchè???? no, va bene, io voglio siano negativi, anche perchè devo cambiare di segno e questo è un bene
#allora perchè mi vengono positivi??


part_k=array(0,dim(theta)[2])

for (i in 1:dim(theta)[1]) {
  part_k[2]=k[2]+sum(cos(c(theta[i,2]-mStore[i,2]))- A1FirstDerivative(kStore[i,2])+
                     sin(c(theta[i,2]-kStore[i,2]))*L[2,1]*sin(c(theta[i,2]-mu[2])))/kStore[i,2]

}

for (j in 3:dim(theta)[2]) {
  for (i in 1:dim(theta)[1]) {
    sum_k=sum(cos(c(theta[i,j]-mStore[i,j]))- A1FirstDerivative(kStore[i,j])+
         sin(c(theta[i,j]-kStore[i,j]))*L[j,1:(j-1)]*rapply(c(theta[i,1:(j-1)]-mu[1:(j-1)]), sin))/kStore[i,j]

  }
  

part_k[j]=k[j]+sum_k
}



########Mu_MU_MU_MU_MU_MU_
part_m=array(0,dim(theta)[2])
for (j in 1:dim(theta)[2]) {
  
for (i in 1:dim(theta)[1]) {
  part_m[j]= - kStore[i,j]*sin(c(theta[i,j]-mStore[i,j]))
}
}

whole=list(part_k=part_k, part_L=part_L, part_m=part_m  )
return(whole)
}

start.time <- Sys.time()
a=trust(theta=Phi1, mu=medie, k=kappino, L=lappo)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
estK=-a[[1]]*t
estL=-a[[2]]*t
estM=-a[[3]]*t

#MANCANO le derivate degli shrink e gli aggiornamenti

#aggiornamenti
start.time <- Sys.time()
t=2.5
alpha=0.5
b=trust(theta=Phi1, mu=a[[1]], k=a[[3]], L=a[[2]])
estK=estK-t*(b[[1]]+(1/(1+exp((alpha)*b[[1]])+(exp(alpha*b[[1]])/(1+exp(alpha*b[[1]])))))     )
estL=estL-t*(b[[2]]+(1/(1+exp((alpha)*b[[2]])+(exp(alpha*b[[2]])/(1+exp(alpha*b[[2]])))))     )
estM=estM-t*(b[[3]]+(1/(1+exp((alpha)*b[[3]])+(exp(alpha*b[[3]])/(1+exp(alpha*b[[3]])))))     )
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken



#mi vengono negativi i cazzo dei k









