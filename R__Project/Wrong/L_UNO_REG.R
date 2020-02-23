L1=function(dati, mi, lap, K, B=2, t=0.3, alpha=0.2,w=2){
  SK=K
  SL=lap
  
  for (q in 1:B) {
    if(w==1){
    for (r in 2:length(K)) {
      SK[r]=SK[r]+t* (diff.k.r(dati, mi, K, r, SL))
    }}
    if(w==0){
    for (r in 1:length(K)) {
      for (s in 1:(r-1)) {
        SL[r,s]=SL[r,s]-t*(diff.l.rs(dati, mi, SK, r, SL)+((1/(1+exp(alpha*SL[r,s])))+(exp(alpha*SL[r,s])/(1+exp(alpha*SL[r,s])))))
      }
          
    }}

    }
  if(w==1){return(SK)}
  if(w==0){return(SL)}
}

#la derivata rispetto a t??? Con quella ci controllo la convergenza?

genio=L1(Phi1, mi=medie,lap=lappo,K=kappino,w=1)


#+((1/(1+exp((alpha)*numeric(SK[r]))))+(exp(alpha*numeric(SK[r]))/(1+exp(alpha*numeric(SK[r]))))))