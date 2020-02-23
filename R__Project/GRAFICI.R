#grafo correlazioni####
dev.new()
ggcorrplot(pp, hc.order = T, method = "square", "lower", lab=T, show.legend = F)
ggcorrplot(pp, hc.order = T, method = "circle", "lower")
pp=cor.circular(Phi2)
diag(pp)=0
pp[lower.tri(pp)]=0
colnames(pp)=row.names(pp)=colnames(Phi2)
colnames(pp)[8]=row.names(pp)[8]="PRO1"

colnames(pp)[9]=row.names(pp)[9]="CYS1"
library(corrplot)
corrplot(pp, type="upper", order="hclust", method = "ellipse",
         tl.col="black", tl.srt=45, #Text label color and rotation
         # hide correlation coefficient on the principal diagonal
         diag=FALSE, col=rainbow(10))

#grafico matrice di adiacenza dei due modelli 
make_A <- function(x, p ){
  A=matrix(0, nrow = 8, ncol = 8)
 
  for (i in 1:8) {
   if(p==1) A[i,]=c(rapply(x[[i]], function(x) as.numeric(x)  ), rep(0,8-i))
   if(p!=1) A[i,]=c(rapply(x[[i]], function(x) as.numeric(x)  ), rep(0,c(0:7)[i]))
  }
  #tolgo i cappi
  diag(A)=0
  colnames(A)=row.names(A)=colnames(Phi2)
  colnames(A)[7]=row.names(A)[7]="PRO1"
  return(A)
}

A_L2=make_A(flag, 1)
A_Ls=matrix(1, 8,8)
diag(A_Ls)=0
A_Ls[lower.tri(A_Ls)]=0
colnames(A_Ls)=row.names(A_Ls)=colnames(Phi2)
colnames(A_Ls)[7]=row.names(A_Ls)[7]="PRO1"

#torno indietro
A_L3=make_A(flag2,0)
A_L3
A_Ls2=matrix(1, 8,8)
diag(A_Ls2)=0
A_Ls2[lower.tri(A_Ls2)]=0
colnames(A_Ls2)=row.names(A_Ls2)=colnames(Phi2)
colnames(A_Ls2)[7]=row.names(A_Ls2)[7]="PRO1"


ggcorrplot(A_Ls, hc.order = T, method = "circle", "lower", show.legend = F, show.diag = F)
ggcorrplot(A_Ls, hc.order = T, method = "circle", "upper", show.legend = F, show.diag = F)
ggcorrplot(A_L2, hc.order = T, method = "circle", show.legend = F, show.diag = F)
ggcorrplot(A_L3, hc.order = T, method = "circle",  show.legend = F, show.diag = F)


A_L2

plot(Phi2)
