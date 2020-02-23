####Test sui parametri####
#sistemo Stime
Stime[[1]]=list( StimePar=k1$par, DiagHess=diag(solve(k1$hessian) ), logL=k1$value )
Stime[[2]]=list( StimePar=k2$par, DiagHess=diag(solve(k2$hessian) ), logL=k2$value )
#calcolo le statistiche test
#la mia hessiana è negativa, dalla 51
Stat=NULL
for (p in 1:9) {
  Stat[[p]]=list(Stime[[p]]$StimePar/(sqrt(Stime[[p]]$DiagHess)))
}

Stat[[2]]
#scelgo alpha
alpha=qnorm(0.05,mean = 0,1)
alpha
#alpha=qnorm(0.025,mean = 0,1)

#creo il vettore che mi dice quali devo prendere
flag=prova=NULL

#faccio il confronto
for (p in 1:9) {
  prova[[p]]=list(dnorm(Stat[[p]][[1]],0,1)  )
 flag[[p]]=list(!c(alpha<= Stat[[p]][[1]] & Stat[[p]][[1]] <=-alpha) )
}


prova
#True=rifiuto h0
(flag)

for (p in 1:9) {
   flag[[p]][[1]][1]=TRUE
}

###LR TEST#####

#Log-likelihood modello intero
Ls=0
for (p in 1:9) {
  Ls=Ls+Stime[[p]]$logL
}
Ls=-Ls
Ls
#modello 2




L2=0
for (p in 1:9) {
  
    b=array(Stime[[p]]$StimePar)
    C=c(array(flag[[p]]))[[1]]
      a=c(c(array(flag[[p]]))[[1]], rep(FALSE, 9-p))
  if(sum(a[-1])==0){
    k=b[1]
      L2=L2+-50*log(2*pi)+ sum(-log( besselI(k,0) )+ k*cos(Phi2[,p]))
  }else if ( sum(a[-1])!=0  ){
    
    cmu=atan((b[1])^(-1)*(b[C]%*%t(sin(Phi2[,a]) )))
    ckappa=sqrt( b[1]^2+ (  b[C]%*%t(sin(Phi2[,a]) ) )^2  )
    cmu=apply(cmu, 1, circular)
    llike=0
    for (i  in 1:n) {
      llike=llike+dvonmises(Phi2[i,p], cmu[i], ckappa[i], log=T)
    }
    L2=L2+llike
  } 
  
    
   }
L2      


#LR test#####
L=-2*(L2-Ls)
L
#p-value=0
#H_1 : il modello più grande è meglio
#H_0 : il modello più semplice è il migliore
#chi^2_{45- 6=39}
sum(rapply(flag, function(x) sum(!x)))
sum(rapply(flag, function(x) sum(x)))
#L2 ha 6 parametri
#Ls ha 45 parametri
p.value=(dchisq(L, 33, ncp = 0, log = FALSE)) 
p.value

# quindi rifiuto H0, ovvero prendo il modello L2



#AIC /BIC TEST
AIC_Ls=2*45-2*(Ls)
AIC_Ls
45*log(50)-2*Ls


AIC_L2=2*sum(rapply(flagp, function(x) sum(!x)))-2*(L2)
AIC_L2
sum(rapply(flagp, function(x) sum(!x)))*log(50)-2*L2

#smaller is best!
