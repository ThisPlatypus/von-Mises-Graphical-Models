SK=kappino
SL=lappo

for (q in 1:2) {
  i=0
     for (r in 2:length(kappino)) {
      SK[r]=SK[r]-0.5* (diff.k.r(Phi1, medie, kappino, r, lappo))
      i=i+1
     }
  i
  
  }  




start.time <- Sys.time()
genio=L1(Phi1, mi=medie,lap=lappo,K=kappino,w=1)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken













