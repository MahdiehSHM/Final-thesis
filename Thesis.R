############################################
# desert endophytes
############################################
library(vegan)
library(mvabund)
library(MASS)
library(ggplot2)
library(ggtree)
library(tidyverse)
library(tidytree)
library(gtable)
library(grid)
library(magrittr)
library(scales)
library(reshape2)
library (lattice)
library(reshape)


data <- read.csv("MSHM.csv",header = T, row.names = 1)
otu <- read.csv(file="OTU.csv",header = T, row.names = 1)
str(data)
str(otu)
summary(otu)
summary(data)
colSums(otu)

######## Temperature as fector:
data$TEMP<-as.factor(data$TEMP)
######## MERG the OTU abundance data by 4
S1<-seq(1,121920,4)
S2<-seq(4,121920,4)
O1<-matrix(0,length(S1),133)
for (i in 1:length(S1)) {
  O1[i,]<-colSums(otu[S1[i]:S2[i],])}

########now convert to data frame and rename the columns
OTUabund<-data.frame(O1)
class(OTUabund)
colnames(OTUabund)=colnames(otu)
name<- row.names(data)
row.names(OTUabund)<- name[1:30480]
summary(OTUabund)
########the abundance data is now saved under object name: OTUabund

######## Merging METADATA
S1<-seq(1,121920,4)
S2<-seq(4,121920,4)
D<-matrix(0,length(S1),7)
for (i in 1:length(S1)) {
  D[i,1]<-noquote(paste(data[S2[i],1]))
  D[i,2]<-noquote(paste(data[S2[i],2]))
  D[i,3]<-noquote(paste(data[S2[i],3]))
  D[i,4]<-noquote(paste(data[S2[i],4]))
  D[i,5]<-noquote(paste(data[S2[i],5]))
  D[i,6]<-noquote(paste(data[S2[i],6]))
  D[i,7]<-noquote(paste(data[S2[i],7]))
}
MetaData<-data.frame(D)
class(MetaData)
colnames(MetaData)=colnames(data)
row.names(MetaData)<- name[1:30480]
summary(MetaData)
######## the merged metadata is now saved under object named: MetaData

######## the merged metadata is now saved under object named: MetaData

######## calculating Isolation Rate for each sample, We need this for diversity estimations
isolationRate= apply(OTUabund,1, sum)
MetaData = cbind(MetaData, IR = isolationRate)
######## now check the MetaData object to see that a new variable is added "IR"
###Exporting as csv files

write.csv(MetaData,file="MataDataMerg.csv")
write.csv(OTUabund, file="OTUabundMerg.csv")
##############################################
##### IR histograms for each variable/factor
##############################################

hist(MetaData $ IR) 
hist(log(MetaData $ IR))
hist(log(MetaData $ IR), prob=TRUE)
hist(log(MetaData $ IR), prob=TRUE, breaks=20)

boxplot(IR ~ SITE, data = MetaData)
boxplot(IR ~ SOIL, data = MetaData)
boxplot(IR ~ TIME, data = MetaData)
boxplot(IR ~ HOST, data = MetaData)
boxplot(IR ~ TISSUE, data = MetaData)
boxplot(IR ~ TEMP, data = MetaData)
boxplot(IR ~ MEDIA, data = MetaData)


aggregate(HAE.wh26.Preussia.sp.     ~ MetaData$SOIL, OTUabund, sum)
aggregate(HAE.wh26.Preussia.sp.      ~ MetaData$SITE, OTUabund, sum)
aggregate(HAE.wh26.Preussia.sp.      ~ MetaData$MEDIA, OTUabund, sum)
aggregate(HAE.wh26.Preussia.sp.      ~ MetaData$TEMP, OTUabund, sum)
aggregate(HAE.wh26.Preussia.sp.      ~ MetaData$TIME, OTUabund, sum)
aggregate(HAE.wh26.Preussia.sp.      ~ MetaData$TISSUE, OTUabund, sum)
aggregate(HAE.wh26.Preussia.sp.      ~ MetaData$HOST, OTUabund, sum)


aggregate(. ~ MetaData$SOIL, OTUabund, sum)



##############################################
##### Diversity Indices
##############################################
#Check and see if it worked?
summary(MetaData)

######## Richness (Species number)#################
######## Richness is the nmber of culture observations in the samples. 
######## Some samples had more than one observed species.
hist(isolationRate)

########  For the evaluation of richness we need to remove the samples with zero observations.
NotZero = isolationRate > 0 
########filter for zero-observation samples
######## Keep only the samples with at least one observed species
AbundNotZero=OTUabund[NotZero,]

######## Richness in the samples
Richness = specnumber(AbundNotZero)
hist(Richness)
hist(log(Richness))

######## Remove the samples with zero observation from the metadata
MetaRich = MetaData[NotZero,]

######## Shannon and Simpson indices####################
######## Keep only samples with at least two OTUs
RichNotOne = Richness > 1
AbundNotOne=AbundNotZero[RichNotOne,]

######## This keeps observations with at least two OTUs
MetaNotOne = MetaRich[RichNotOne,] 

######## Calculate diversity indices
shannon = diversity(AbundNotOne,index = "shannon")
simpson = diversity(AbundNotOne,index = "simpson")
hist(shannon)
hist(simpson)
hist(log(shannon))
hist(log(simpson))

####### Subsetting the data for Gypsumsoil
Gypsumsoil = subset (MetaData, SOIL%in%c("Gypsum Soil"))
GypsumsoilOTU = subset (OTUabund, MetaData$SOIL %in% c("Gypsum Soil"))

OTUsINGypsumsoil<-colSums(GypsumsoilOTU)
GypsumsoilOTU<-GypsumsoilOTU[, colSums(GypsumsoilOTU != 0) > 0]
rownames(GypsumsoilOTU)==rownames(Gypsumsoil)
colSums(GypsumsoilOTU)
colnames(GypsumsoilOTU)

View(Gypsumsoil)
View(GypsumsoilOTU)

levels(Gypsumsoil$SOIL)

####### Subsetting the data for Salinesoil
Salinesoil = subset (MetaData, SOIL%in%c("Saline Soil"))
SalinesoilOTU = subset (OTUabund, MetaData$SOIL %in% c("Saline Soil"))

OTUsINSalinesoil<-colSums(SalinesoilOTU)
SalinesoilOTU<-SalinesoilOTU[, colSums(SalinesoilOTU != 0) > 0]
rownames(SalinesoilOTU)==rownames(Salinesoil)
colSums(SalinesoilOTU)
colnames(SalinesoilOTU)

View(Salinesoil)
View(SalinesoilOTU)


####### Subsetting the data for Aridsoil

Aridsoil = subset (MetaData, SOIL%in%c("Arid soil"))
AridsoilOTU = subset (OTUabund, MetaData$SOIL %in% c("Arid soil"))

rownames(AridsoilOTU)==rownames(Aridsoil)
class(Aridsoil)
class(AridsoilOTU)
View(Aridsoil)
View(AridsoilOTU)

colSums(AridsoilOTU)

AridsoilOTU$GME.ss9..Quambalaria.cyanescens<-NULL
AridsoilOTU$PFE.sh7..Rosellinia.limonispora <-NULL
AridsoilOTU$LREwh64..Neocamarosporium.chichastianum<-NULL
AridsoilOTU$SREwh22.Preussia.minimoides<-NULL
AridsoilOTU$THE.ss3.Fusarium.sp.<-NULL
AridsoilOTU$SKE.ws4.botryotrichum.spirotrichum <-NULL
AridsoilOTU$PSE.wh14.Preussia.sp. <-NULL
AridsoilOTU$PSE.wh40.Dimorphosporicola.tragani<-NULL
AridsoilOTU$PSE.wh66.Comoclathris.italica <-NULL
AridsoilOTU$MUE.ss1.Marasmiellus.candidus<-NULL
AridsoilOTU$MUE.ss9.Eucasphaeria.capensis<-NULL
AridsoilOTU$LRE.wh34.Alternaria.quercicola<-NULL
AridsoilOTU$LREwh63.Coniophora.olivacea<-NULL
AridsoilOTU$HAE.wh26.Preussia.sp.<-NULL
AridsoilOTU$HAE.wh10.Raffaelea.montetyi<-NULL
AridsoilOTU$HAE.wh65.Coniophora.marmorata<-NULL
AridsoilOTU$HBE.wh23.Kalmusia.spartii <-NULL
AridsoilOTU$HBE.wh27.Chaetomium.truncatulum<-NULL
AridsoilOTU$GME.ss5.Aspergillus.frequens<-NULL
AridsoilOTU$GAE.ws1.Tricharina.striispora <-NULL
AridsoilOTU$GAE.ss.6.Fusarium.ensiforme<-NULL
AridsoilOTU$GAE.ws4.Stemphylium.vesicarium<-NULL
AridsoilOTU$ERE.ws2.Fusarium.sp.<-NULL
AridsoilOTU$ERE.ws1.Preussia.intermedia<-NULL
AridsoilOTU$ERE.ws4.Talaromyces.pinophilus<-NULL
AridsoilOTU$DHE.ws22.Meyerozyma.sp.<-NULL
AridsoilOTU$DHE.ss2.Alternaria.arborescens<-NULL
AridsoilOTU$AseE.sh3.Microdochium.bolleyi<-NULL
AridsoilOTU$SKE.ss8.Penicillium.canescens<-NULL
AridsoilOTU$HspE.ws3.Humicola.grisea<-NULL
AridsoilOTU$SREwh19.Neocamarosporium.goegapense<-NULL
AridsoilOTU$AVE.sh2.Aspergillus.conversis<-NULL
AridsoilOTU$APEsh6.Dictyosporium.digitatum<-NULL
AridsoilOTU$APE.sh8.Pestalotiopsis.vismiae<-NULL
AridsoilOTU$APE.sh5.Nigrospora.sphaerica<-NULL
AridsoilOTU$TspE.ws5.Chaetomium.spiculipilium<-NULL
AridsoilOTU$SREwh18...Preussia.sp.<-NULL
AridsoilOTU$RAE.sh8.Aspergillus.protuberus<-NULL
AridsoilOTU$RAE.sh3.Aspergillus.nidulans<-NULL
AridsoilOTU$RAE.sh2.Paraconiothyrium.brasiliense<-NULL
AridsoilOTU$TAEwh43.Chaetomium.longicolleum<-NULL
AridsoilOTU$TAE.sh8.Aspergillus.nomius  <-NULL
AridsoilOTU$TAE.sh9.Fusarium.sp.<-NULL
AridsoilOTU$TAEwh25.Alternaria.tenuissimas<-NULL
AridsoilOTU$TAEwh9.Comoclathris.spartii<-NULL
AridsoilOTU$TAEwh20.Thielavia.subthermophila<-NULL
AridsoilOTU$TPEsh53.Thielavia.arenaria<-NULL
AridsoilOTU$LAE.SH.7.Muriphaeosphaeria.viburni<-NULL
AridsoilOTU$LAE.sh1.Acrocalymma.sp.<-NULL
AridsoilOTU$SRE.sh30..Preussia.grandispora<-NULL
AridsoilOTU$PFE.sh4.Coniolariella.sp.<-NULL
AridsoilOTU$TPE.sh31.Chaetomium.interruptum<-NULL
AridsoilOTU$TPE.ss.11.Chaetomium.subaffine<-NULL
AridsoilOTU$TPE.ss.9.Acremonium.sclerotigenum<-NULL
AridsoilOTU$TPE.ss8.Aspergillus.sydowii<-NULL
AridsoilOTU$TPE.wh6..Coniothyrina.agaves<-NULL
AridsoilOTU$TPEsh44.Cladorrhinum.phialophoroides<-NULL
AridsoilOTU$TPEsh17.Sarocladium.sp.<-NULL
AridsoilOTU$ESE.wh29.Schizothecium.inaequale<-NULL
AridsoilOTU$SRE.ss.4.Neocamarosporium.sp.<-NULL
AridsoilOTU$SRE.ws8.Botryotrichum.murorum<-NULL
AridsoilOTU$SRE.ws10.Sarocladium.kiliense<-NULL
AridsoilOTU$SRE.sh7.Chaetomium.cucumericola<-NULL
AridsoilOTU$SRE.sh5.Fusarium.redolens<-NULL
AridsoilOTU$SRE.sh9.Preussia.intermedia <-NULL
AridsoilOTU$SRE.sh4.Penicillium.vinaceum<-NULL
AridsoilOTU$GME.ws1.Preussia.Africana<-NULL
AridsoilOTU$GAE.ss3.Aspergillus.jensenii<-NULL
AridsoilOTU$GAE.ss7.Aspergillus.conversis<-NULL
AridsoilOTU$GAE.ws8.Sporormiella.septenaria<-NULL
AridsoilOTU$GspE.ss6.Talaromyces.cecidicola<-NULL
AridsoilOTU$GspE.ss5.Talaromyces.leycettanus<-NULL
AridsoilOTU$TAEwh25.Alternaria.tenuissima<-NULL
AridsoilOTU$SRE.sh3.Trichoderma.rifaii<-NULL

# OTU frequencies Article 1
colSums(AridsoilOTU)
#OTU list Article 1
colnames(AridsoilOTU)

View(AridsoilOTU)
View(Aridsoil)

levels(Aridsoil$SOIL)

######## calculating Isolation Rate for each sample, We need this for diversity estimations
isolationRateArid= apply(AridsoilOTU,1, sum)
MetaDataArid = cbind(Aridsoil, IR = isolationRateArid)

#######################################
#IR model
#######################################
hist(Aridsoil$IR)
IR.m<-glm(formula =IR~HOST+TIME+SITE+,data = Aridsoil,
          family=poisson(link = "log"))











