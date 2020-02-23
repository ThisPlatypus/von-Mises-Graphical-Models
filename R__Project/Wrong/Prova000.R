#####################################
#  prova gestione
# angoli delle proteine
#####################################
library(rstudioapi) 
setwd(dirname(getActiveDocumentContext()$path))
getwd()
##################
# Carico i dati 1YT6
library(bio3d)
### Questa proteina ha 10 aminoacidi
#dati <- read.pdb("/Users/gottard/Documents/Anna/Artic/Circular/Data/1YT6PeptideSD50/1yt6.pdb", multi = TRUE)
dati <- read.pdb("D:/PROTEIN/R__Project/Data/1E8LHenLysozyme50/1e8l.pdb", multi = TRUE)
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
angle1e8l<- angoli

My.readme <- rbind(c("Ogni arrey riguarda una proteina"), c("La prima dimensione sono i modelli (replicazioni)"), c("La seconda dimensione sono gli aminoacidi: il nome di riga dice quali"),c("La terza dimensione concerne gli angoli"), c("NB: non tutti gli angoli ci sono per tutti gli aminoacidi"), c("Per vedere, prova angle1YT6[1,,] ")) 

save(angle1YT6,angle1e8l,My.readme, file="AngleProtein.RData")
load("AngleProtein.RData")
#############################################
### esempio per un singolo aminoacido
library(ggm)
dati$seqres # nomi aminoacidi
angoli.LEU <- as.matrix(angoli[,4,-c(5:7)], nrow=50)
round(parcor(var(angoli.LEU)),2)

angoli.ASP <- as.matrix(angoli[,7,-c(5:7)], nrow=50)
round(parcor(var(angoli.ASP)),2)

angoli.QT <- as.matrix(angoli[, ,-c(5:7)], nrow=500)#quasitutti, quelli cha hanno 4 angoli
round(parcor(var(angoli.TUTTI)),2)

pippo <- matrix(angoli[,,-c(5:7)],500,4)
head(pippo)
dim(pippo)

