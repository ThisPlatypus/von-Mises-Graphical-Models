library(rstudioapi) 
library(circular)
setwd(dirname(getActiveDocumentContext()$path))
getwd()
##################
# Carico i dati 1YT6
library(bio3d)
### Questa proteina ha 10 aminoacidi
#dati <- read.pdb("/Users/gottard/Documents/Anna/Artic/Circular/Data/1YT6PeptideSD50/1yt6.pdb", multi = TRUE)
dati <- read.pdb("D:/PROTEIN/R__Project/Data/1YT6PeptideSD50/1yt6.pdb", multi = TRUE)
dati
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

#angle1YT6 <- angoli
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
  
  