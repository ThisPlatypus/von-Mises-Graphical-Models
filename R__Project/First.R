# Import librerie
library(rstudioapi) 
library(circular)
library(bio3d)

# Imposto la cartella di lavoro
setwd(dirname(getActiveDocumentContext()$path))
getwd()
##################

# Carico i dati 1YT6
dati <- read.pdb("D:/PROTEIN/R__Project/Data/1YT6PeptideSD50/1yt6.pdb", multi = TRUE)

names(dati)
dati$atom$resno # indicatore dell'aminoacido 
dati$seqres # sequenza aminoacidi
namin<-length(dati$seqres)
dim(dati$xyz) ###nb le righe sono i modelli, le colonne sono le 3 coordinate per gli  atomi della proteina

dati1<- dati
angoli <- array(NA, c(50, namin, 7)) # model, aminoacido, angolo
for (i in 1:50){
  dati1$xyz <- dati$xyz[i,]
  tor<-torsion.pdb(dati1)
  angoli[i,,] <- tor$tbl
}
angoli[1,,]

dimnames(angoli)[[1]] <- paste("Mod", 1:50, sep="")
dimnames(angoli)[[2]] <- as.vector(dati$seqres)
dimnames(angoli)[[3]] <- colnames(tor$tbl)


angle<- angoli



#tolgo la prima colonna perchè ha tutti NA
Phi=angle[,,"phi"]
Phi=Phi[,-1]
dev.new()
ggpairs(as.data.frame(Phi),mapping = aes(color = "#B3CDE3"), diag=NULL, upper = list(continuous = wrap("density", alpha = 0.5)))


Phi1=as.circular(Phi, type="angles", units="radians", zero=circular(0), rotation="counter", template="none")

#centro i dati
  make_centred=function(data){
    mi=circular(0)
    for (i in 1:dim(data)[2]) {
      data[,i]=data[,i]-mean.circular(data[,i], units="radians")
    }
    return(data)
  }

  
  Phi2 <- make_centred(Phi1)

  dev.new()
  plot(Phi2[,2], stack=TRUE, bins=50, main="Scatter plot of 1yt6", col="#B3CDE3" , tcl.text=1.65, shrink=1.5) 
  
  
  lines(density.circular(Phi2, bw=50), shrink = 1.3, main = "Density of 1yt6 data", col="#FBB4AE" )
  

###############################################################################################################################

#Grafici

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

###############################################################################################################################

#Starter value

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
###############################################################################################################################

# Problema bivariato

count=0
Stime=matrix(0, ncol = 4, nrow = 128)

Stime[1:2,1]=c(k1$par, k2.L12$par[1])
Stime[1:2,2]=c(0,k2.L12$par[2])

for (p in 3:128) {  

# Funzione da ottimizzare
fr<- function(parms, x, j, memo){ 
   -(   -50*j*log(2*pi)+memo+sum(  -log(besselI( parms[1] ,0)) + parms[1]*cos(x[,j] - parms[2])   )  )
   
}

#Funzione gradiente

M <- function(hat.K, hat.L, X ,j ){
  A=0
  A=sum( - log(besselI( hat.K[1] ,0))  + hat.K[1]*cos(X[,1]) ) 
  
    for (i in 2:(j-1)) {
    B=0
    for (l in 1:(j-1)) {
      B=B+sum( hat.L[l]*sin(X[,i])  )
    }
    
    A=A+ sum( - log(  besselI(  sqrt(hat.K[i]^2+B^2) ,0))     + sqrt(hat.K[i]^2+B^2)*cos( X[,i] - atan( B/ hat.K[i] )  ))
    
  }
  return(A)
}

#calcolo gradiente 
g=M(Stime[1:p,1], Stime[1:p,2], Phi2 ,j=p )

# Stime ottime

kj.Lj <- optim(c(0.5,-1), fr,  x=Phi2, j=p, memo=g, #method = "BFGS")
               method="L-BFGS-B", hessian = 1,
               lower = c(0, rep(-Inf, 2)), 
               upper = c(Inf, rep(Inf, 2)))

#check convergence
count=count+kj.Lj$convergence

#estimates
Stime[p,1]=kj.Lj$par[1]
Stime[p,2]=kj.Lj$par[2]
Stime[p, 3:4]=diag(solve(kj.Lj$hessian))
}

#check
count
Stime
round(stime1[,1],3)==round(Stime[,1],3)
round(stime1[,2],3)==round(Stime[,2],3)
stime1=Stime

###################################################################################################

# Problema multivariato  

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

for (p in 3:9) {
  k.l<- optim(c(rep(0.2,p)), loglike, Phi2=Phi2, j=p, hessian = T, lower = c(0, rep(-Inf, p-1)), upper = c(rep(Inf,p)))
  Stime[[p]] <- list( StimePar=k.l$par, DiagHess=diag(solve(k.l$hessian) ), logL=k.l$value )
  
}

save(Stime, file="fname.RData")
























