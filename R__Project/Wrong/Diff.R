#Differenziale per un k.r####
diff.k.r = function(dati, mi, K, R, lap) {
  S = circular(0)
  for (i in 1:dim(dati)[1]) {
    a=dati[i, R] - mi.r(dati[i,], lap, K, mi, R)
    S = S + (
      cos(a) -A1FirstDerivative(K[r]) ) / k.r(dati[i,], lap, K, mi, R) +(
        sin(dati[i, R] - mi.r(dati[i,], lap, K, mi, R)) * Sp(lap, dati[i,], mi, R)
      )
    
  }
  
  
  d.k = K[R] * S
  
  return(d.k)
}




#Questa funione mi da la media al netto per ogni passo i
#devi vedere se devi storarlo, così non lo stai facendo
mi.r = function(angle, lap, k, mi, R) {
  P = circular(0)
  for (l in 1:(R-1)) {
    P = P + ((lap[R, l] * sin(angle[l]) - mi[l]) / k[R])
  }
  
  mi.r = mi[R] + atan(P)
  
  return(mi.r)
}

#Questa funione mi da la k al netto per ogni passo i
#devi vedere se devi storarlo, così non lo stai facendo

k.r = function(angle, lap, k, mi, R) {
  P = circular(0)
  for (l in 1:(R - 1)) {
    P = P + lap[R, l] * sin(angle[l] - mi[l])
  }
  
  k.r = sqrt(k ^ 2 + P ^ 2)
  return(k.r)
}



###Funzione che mi da parte del differenziale di k.r

Sp = function(lap, angle, mi, R) {
  Q = circular(0)
  for (l in 1:(R-1)) {
    Q = Q + lap[R, l] * sin(angle[l] - mi[l])
  }
  return(Q)
}

#Differenziale per un l.(r,s)####
diff.l.rs = function(dati, mi, K, R, lap) {
  d.lrs = circular(0)
  for (j in 1:dim(dati)[2]) {
    for (i in 1:dim(dati)[1]) {
      d.lrs = d.lrs + PKL(dati[i,], lap, k, mi, j) * (-1 * A1FirstDerivative(K)) + cos(dati[i, j] - mi.r(dati[i,], lap, K[j], mi, j)) + PML(dati[i,], lap, K, mi, j) *
        (k.r(dati[i,], lap, K[j], mi, j) * sin(dati[i, j] - mi.r(dati[i,], lap, K[j], mi, j)))
    }
    
    
  }
  
  return(d.lrs)
  
}



#dK/dLam
PKL = function(angle, lap, k, mi, R) {
  a = circular(0)
  c = k.r(angle, lap, k, mi, R)
  for (l in 1:(R - 1)) {
    a = a + (lap[R, l] * sin(angle[l] - mi[l]) * sin(angle[R] - mi[R])) / c
  }
  return(a)
}


PML=function(angle, lap, k, mi, R){
  b=circular(0)
  for (l in 1:(R-1)) {
    b=b+((lap[R,l]*sin(angle[l]-mi[l]))/k[R])
  }

  a=(sin(angle[R]-mi[R])/(k[R]+ (1+b^2) ))
  return(a)
  
  }

