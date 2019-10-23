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


aggregate(APE.se5.Staphylotrichum.coccosporum     ~ MetaData$SOIL, OTUabund, sum)
aggregate(APE.se5.Staphylotrichum.coccosporum     ~ MetaData$SITE, OTUabund, sum)
aggregate(APE.se5.Staphylotrichum.coccosporum     ~ MetaData$MEDIA, OTUabund, sum)
aggregate(APE.se5.Staphylotrichum.coccosporum     ~ MetaData$TEMP, OTUabund, sum)
aggregate(APE.se5.Staphylotrichum.coccosporum     ~ MetaData$TIME, OTUabund, sum)
aggregate(APE.se5.Staphylotrichum.coccosporum     ~ MetaData$TISSUE, OTUabund, sum)
aggregate(APE.se5.Staphylotrichum.coccosporum     ~ MetaData$HOST, OTUabund, sum)


aggregate(TPEsh28.Humicola.fuscoatra      ~ MetaData$SOIL, OTUabund, sum)
aggregate(TPEsh28.Humicola.fuscoatra      ~ MetaData$SITE, OTUabund, sum)
aggregate(TPEsh28.Humicola.fuscoatra      ~ MetaData$MEDIA, OTUabund, sum)
aggregate(TPEsh28.Humicola.fuscoatra      ~ MetaData$TEMP, OTUabund, sum)
aggregate(TPEsh28.Humicola.fuscoatra      ~ MetaData$TIME, OTUabund, sum)
aggregate(TPEsh28.Humicola.fuscoatra      ~ MetaData$TISSUE, OTUabund, sum)
aggregate(TPEsh28.Humicola.fuscoatra      ~ MetaData$HOST, OTUabund, sum)



aggregate(ZEE.se11.Periconia.macrospinosa     ~ MetaData$SOIL, OTUabund, sum)
aggregate(ZEE.se11.Periconia.macrospinosa     ~ MetaData$SITE, OTUabund, sum)
aggregate(ZEE.se11.Periconia.macrospinosa     ~ MetaData$MEDIA, OTUabund, sum)
aggregate(ZEE.se11.Periconia.macrospinosa     ~ MetaData$TEMP, OTUabund, sum)
aggregate(ZEE.se11.Periconia.macrospinosa     ~ MetaData$TIME, OTUabund, sum)
aggregate(ZEE.se11.Periconia.macrospinosa     ~ MetaData$TISSUE, OTUabund, sum)
aggregate(ZEE.se11.Periconia.macrospinosa     ~ MetaData$HOST, OTUabund, sum)



aggregate(GME.ss9..Quambalaria.cyanescens     ~ MetaData$SOIL, OTUabund, sum)
aggregate(GME.ss9..Quambalaria.cyanescens     ~ MetaData$SITE, OTUabund, sum)
aggregate(GME.ss9..Quambalaria.cyanescens     ~ MetaData$MEDIA, OTUabund, sum)
aggregate(GME.ss9..Quambalaria.cyanescens     ~ MetaData$TEMP, OTUabund, sum)
aggregate(GME.ss9..Quambalaria.cyanescens     ~ MetaData$TIME, OTUabund, sum)
aggregate(GME.ss9..Quambalaria.cyanescens     ~ MetaData$TISSUE, OTUabund, sum)
aggregate(GME.ss9..Quambalaria.cyanescens     ~ MetaData$HOST, OTUabund, sum)



aggregate(PFE.sh7..Rosellinia.limonispora     ~ MetaData$SOIL, OTUabund, sum)
aggregate(PFE.sh7..Rosellinia.limonispora     ~ MetaData$SITE, OTUabund, sum)
aggregate(PFE.sh7..Rosellinia.limonispora     ~ MetaData$MEDIA, OTUabund, sum)
aggregate(PFE.sh7..Rosellinia.limonispora     ~ MetaData$TEMP, OTUabund, sum)
aggregate(PFE.sh7..Rosellinia.limonispora     ~ MetaData$TIME, OTUabund, sum)
aggregate(PFE.sh7..Rosellinia.limonispora     ~ MetaData$TISSUE, OTUabund, sum)
aggregate(PFE.sh7..Rosellinia.limonispora     ~ MetaData$HOST, OTUabund, sum)




aggregate(LREwh64..Neocamarosporium.chichastianum     ~ MetaData$SOIL, OTUabund, sum)
aggregate(LREwh64..Neocamarosporium.chichastianum     ~ MetaData$SITE, OTUabund, sum)
aggregate(LREwh64..Neocamarosporium.chichastianum     ~ MetaData$MEDIA, OTUabund, sum)
aggregate(LREwh64..Neocamarosporium.chichastianum     ~ MetaData$TEMP, OTUabund, sum)
aggregate(LREwh64..Neocamarosporium.chichastianum     ~ MetaData$TIME, OTUabund, sum)
aggregate(LREwh64..Neocamarosporium.chichastianum     ~ MetaData$TISSUE, OTUabund, sum)
aggregate(LREwh64..Neocamarosporium.chichastianum     ~ MetaData$HOST, OTUabund, sum)


aggregate(SREwh22.Preussia.minimoides     ~ MetaData$SOIL, OTUabund, sum)
aggregate(SREwh22.Preussia.minimoides     ~ MetaData$SITE, OTUabund, sum)
aggregate(SREwh22.Preussia.minimoides     ~ MetaData$MEDIA, OTUabund, sum)
aggregate(SREwh22.Preussia.minimoides     ~ MetaData$TEMP, OTUabund, sum)
aggregate(SREwh22.Preussia.minimoides     ~ MetaData$TIME, OTUabund, sum)
aggregate(SREwh22.Preussia.minimoides     ~ MetaData$TISSUE, OTUabund, sum)
aggregate(SREwh22.Preussia.minimoides     ~ MetaData$HOST, OTUabund, sum)
 

aggregate(THE.we10..Aporospora.terricola       ~ MetaData$SOIL, OTUabund, sum)
aggregate(THE.we10..Aporospora.terricola       ~ MetaData$SITE, OTUabund, sum)
aggregate(THE.we10..Aporospora.terricola       ~ MetaData$MEDIA, OTUabund, sum)
aggregate(THE.we10..Aporospora.terricola       ~ MetaData$TEMP, OTUabund, sum)
aggregate(THE.we10..Aporospora.terricola       ~ MetaData$TIME, OTUabund, sum)
aggregate(THE.we10..Aporospora.terricola       ~ MetaData$TISSUE, OTUabund, sum)
aggregate(THE.we10..Aporospora.terricola       ~ MetaData$HOST, OTUabund, sum)

aggregate(THE.ss3.Fusarium.sp.          ~ MetaData$SOIL, OTUabund, sum)
aggregate(THE.ss3.Fusarium.sp.          ~ MetaData$SITE, OTUabund, sum)
aggregate(THE.ss3.Fusarium.sp.          ~ MetaData$MEDIA, OTUabund, sum)
aggregate(THE.ss3.Fusarium.sp.          ~ MetaData$TEMP, OTUabund, sum)
aggregate(THE.ss3.Fusarium.sp.          ~ MetaData$TIME, OTUabund, sum)
aggregate(THE.ss3.Fusarium.sp.          ~ MetaData$TISSUE, OTUabund, sum)
aggregate(THE.ss3.Fusarium.sp.          ~ MetaData$HOST, OTUabund, sum)


aggregate(SAE.se4.Coniophora.prasinoides     ~ MetaData$SOIL, OTUabund, sum)
aggregate(SAE.se4.Coniophora.prasinoides     ~ MetaData$SITE, OTUabund, sum)
aggregate(SAE.se4.Coniophora.prasinoides     ~ MetaData$MEDIA, OTUabund, sum)
aggregate(SAE.se4.Coniophora.prasinoides     ~ MetaData$TEMP, OTUabund, sum)
aggregate(SAE.se4.Coniophora.prasinoides     ~ MetaData$TIME, OTUabund, sum)
aggregate(SAE.se4.Coniophora.prasinoides     ~ MetaData$TISSUE, OTUabund, sum)
aggregate(SAE.se4.Coniophora.prasinoides     ~ MetaData$HOST, OTUabund, sum)


aggregate(SKE.ws4.botryotrichum.spirotrichum          ~ MetaData$SOIL, OTUabund, sum)
aggregate(SKE.ws4.botryotrichum.spirotrichum          ~ MetaData$SITE, OTUabund, sum)
aggregate(SKE.ws4.botryotrichum.spirotrichum          ~ MetaData$MEDIA, OTUabund, sum)
aggregate(SKE.ws4.botryotrichum.spirotrichum          ~ MetaData$TEMP, OTUabund, sum)
aggregate(SKE.ws4.botryotrichum.spirotrichum          ~ MetaData$TIME, OTUabund, sum)
aggregate(SKE.ws4.botryotrichum.spirotrichum          ~ MetaData$TISSUE, OTUabund, sum)
aggregate(SKE.ws4.botryotrichum.spirotrichum          ~ MetaData$HOST, OTUabund, sum)


aggregate(SIE.sh1.Briansuttonomyces.eucalypti      ~ MetaData$SOIL, OTUabund, sum)
aggregate(SIE.sh1.Briansuttonomyces.eucalypti      ~ MetaData$SITE, OTUabund, sum)
aggregate(SIE.sh1.Briansuttonomyces.eucalypti      ~ MetaData$MEDIA, OTUabund, sum)
aggregate(SIE.sh1.Briansuttonomyces.eucalypti      ~ MetaData$TEMP, OTUabund, sum)
aggregate(SIE.sh1.Briansuttonomyces.eucalypti      ~ MetaData$TIME, OTUabund, sum)
aggregate(SIE.sh1.Briansuttonomyces.eucalypti      ~ MetaData$TISSUE, OTUabund, sum)
aggregate(SIE.sh1.Briansuttonomyces.eucalypti      ~ MetaData$HOST, OTUabund, sum)



aggregate(RAE.sh12.Acrocalymma.vagum       ~ MetaData$SOIL, OTUabund, sum)
aggregate(RAE.sh12.Acrocalymma.vagum       ~ MetaData$SITE, OTUabund, sum)
aggregate(RAE.sh12.Acrocalymma.vagum       ~ MetaData$MEDIA, OTUabund, sum)
aggregate(RAE.sh12.Acrocalymma.vagum       ~ MetaData$TEMP, OTUabund, sum)
aggregate(RAE.sh12.Acrocalymma.vagum       ~ MetaData$TIME, OTUabund, sum)
aggregate(RAE.sh12.Acrocalymma.vagum       ~ MetaData$TISSUE, OTUabund, sum)
aggregate(RAE.sh12.Acrocalymma.vagum       ~ MetaData$HOST, OTUabund, sum)



aggregate(PSE.wh14.Preussia.sp.       ~ MetaData$SOIL, OTUabund, sum)
aggregate(PSE.wh14.Preussia.sp.       ~ MetaData$SITE, OTUabund, sum)
aggregate(PSE.wh14.Preussia.sp.       ~ MetaData$MEDIA, OTUabund, sum)
aggregate(PSE.wh14.Preussia.sp.       ~ MetaData$TEMP, OTUabund, sum)
aggregate(PSE.wh14.Preussia.sp.       ~ MetaData$TIME, OTUabund, sum)
aggregate(PSE.wh14.Preussia.sp.       ~ MetaData$TISSUE, OTUabund, sum)
aggregate(PSE.wh14.Preussia.sp.       ~ MetaData$HOST, OTUabund, sum)


aggregate(PSE.wh40.Dimorphosporicola.tragani           ~ MetaData$SOIL, OTUabund, sum)
aggregate(PSE.wh40.Dimorphosporicola.tragani           ~ MetaData$SITE, OTUabund, sum)
aggregate(PSE.wh40.Dimorphosporicola.tragani           ~ MetaData$MEDIA, OTUabund, sum)
aggregate(PSE.wh40.Dimorphosporicola.tragani           ~ MetaData$TEMP, OTUabund, sum)
aggregate(PSE.wh40.Dimorphosporicola.tragani           ~ MetaData$TIME, OTUabund, sum)
aggregate(PSE.wh40.Dimorphosporicola.tragani           ~ MetaData$TISSUE, OTUabund, sum)
aggregate(PSE.wh40.Dimorphosporicola.tragani           ~ MetaData$HOST, OTUabund, sum)



aggregate(PSE.wh66.Comoclathris.italica            ~ MetaData$SOIL, OTUabund, sum)
aggregate(PSE.wh66.Comoclathris.italica            ~ MetaData$SITE, OTUabund, sum)
aggregate(PSE.wh66.Comoclathris.italica            ~ MetaData$MEDIA, OTUabund, sum)
aggregate(PSE.wh66.Comoclathris.italica            ~ MetaData$TEMP, OTUabund, sum)
aggregate(PSE.wh66.Comoclathris.italica            ~ MetaData$TIME, OTUabund, sum)
aggregate(PSE.wh66.Comoclathris.italica            ~ MetaData$TISSUE, OTUabund, sum)
aggregate(PSE.wh66.Comoclathris.italica            ~ MetaData$HOST, OTUabund, sum)



aggregate(PSE.ss7.Penicillium.sp.             ~ MetaData$SOIL, OTUabund, sum)
aggregate(PSE.ss7.Penicillium.sp.             ~ MetaData$SITE, OTUabund, sum)
aggregate(PSE.ss7.Penicillium.sp.             ~ MetaData$MEDIA, OTUabund, sum)
aggregate(PSE.ss7.Penicillium.sp.             ~ MetaData$TEMP, OTUabund, sum)
aggregate(PSE.ss7.Penicillium.sp.             ~ MetaData$TIME, OTUabund, sum)
aggregate(PSE.ss7.Penicillium.sp.             ~ MetaData$TISSUE, OTUabund, sum)
aggregate(PSE.ss7.Penicillium.sp.             ~ MetaData$HOST, OTUabund, sum)


aggregate(PSE.we4..Chaetomium.globosum             ~ MetaData$SOIL, OTUabund, sum)
aggregate(PSE.we4..Chaetomium.globosum             ~ MetaData$SITE, OTUabund, sum)
aggregate(PSE.we4..Chaetomium.globosum             ~ MetaData$MEDIA, OTUabund, sum)
aggregate(PSE.we4..Chaetomium.globosum             ~ MetaData$TEMP, OTUabund, sum)
aggregate(PSE.we4..Chaetomium.globosum             ~ MetaData$TIME, OTUabund, sum)
aggregate(PSE.we4..Chaetomium.globosum             ~ MetaData$TISSUE, OTUabund, sum)
aggregate(PSE.we4..Chaetomium.globosum             ~ MetaData$HOST, OTUabund, sum)


aggregate(PSE.se8.Podospora.minicauda             ~ MetaData$SOIL, OTUabund, sum)
aggregate(PSE.se8.Podospora.minicauda             ~ MetaData$SITE, OTUabund, sum)
aggregate(PSE.se8.Podospora.minicauda             ~ MetaData$MEDIA, OTUabund, sum)
aggregate(PSE.se8.Podospora.minicauda             ~ MetaData$TEMP, OTUabund, sum)
aggregate(PSE.se8.Podospora.minicauda             ~ MetaData$TIME, OTUabund, sum)
aggregate(PSE.se8.Podospora.minicauda             ~ MetaData$TISSUE, OTUabund, sum)
aggregate(PSE.se8.Podospora.minicauda             ~ MetaData$HOST, OTUabund, sum)



aggregate(PSE.we8.Alternaria.chlamydospora        ~ MetaData$SOIL, OTUabund, sum)
aggregate(PSE.we8.Alternaria.chlamydospora        ~ MetaData$SITE, OTUabund, sum)
aggregate(PSE.we8.Alternaria.chlamydospora        ~ MetaData$MEDIA, OTUabund, sum)
aggregate(PSE.we8.Alternaria.chlamydospora        ~ MetaData$TEMP, OTUabund, sum)
aggregate(PSE.we8.Alternaria.chlamydospora        ~ MetaData$TIME, OTUabund, sum)
aggregate(PSE.we8.Alternaria.chlamydospora        ~ MetaData$TISSUE, OTUabund, sum)
aggregate(PSE.we8.Alternaria.chlamydospora        ~ MetaData$HOST, OTUabund, sum)


aggregate(MUE.ss1.Marasmiellus.candidus                 ~ MetaData$SOIL, OTUabund, sum)
aggregate(MUE.ss1.Marasmiellus.candidus                 ~ MetaData$SITE, OTUabund, sum)
aggregate(MUE.ss1.Marasmiellus.candidus                 ~ MetaData$MEDIA, OTUabund, sum)
aggregate(MUE.ss1.Marasmiellus.candidus                 ~ MetaData$TEMP, OTUabund, sum)
aggregate(MUE.ss1.Marasmiellus.candidus                 ~ MetaData$TIME, OTUabund, sum)
aggregate(MUE.ss1.Marasmiellus.candidus                 ~ MetaData$TISSUE, OTUabund, sum)
aggregate(MUE.ss1.Marasmiellus.candidus                 ~ MetaData$HOST, OTUabund, sum)



aggregate(MUE.ss9.Eucasphaeria.capensis                 ~ MetaData$SOIL, OTUabund, sum)
aggregate(MUE.ss9.Eucasphaeria.capensis                 ~ MetaData$SITE, OTUabund, sum)
aggregate(MUE.ss9.Eucasphaeria.capensis                 ~ MetaData$MEDIA, OTUabund, sum)
aggregate(MUE.ss9.Eucasphaeria.capensis                 ~ MetaData$TEMP, OTUabund, sum)
aggregate(MUE.ss9.Eucasphaeria.capensis                 ~ MetaData$TIME, OTUabund, sum)
aggregate(MUE.ss9.Eucasphaeria.capensis                 ~ MetaData$TISSUE, OTUabund, sum)
aggregate(MUE.ss9.Eucasphaeria.capensis                 ~ MetaData$HOST, OTUabund, sum)


aggregate(LREwh63.Coniophora.olivacea                    ~ MetaData$SOIL, OTUabund, sum)
aggregate(LREwh63.Coniophora.olivacea                    ~ MetaData$SITE, OTUabund, sum)
aggregate(LREwh63.Coniophora.olivacea                    ~ MetaData$MEDIA, OTUabund, sum)
aggregate(LREwh63.Coniophora.olivacea                    ~ MetaData$TEMP, OTUabund, sum)
aggregate(LREwh63.Coniophora.olivacea                    ~ MetaData$TIME, OTUabund, sum)
aggregate(LREwh63.Coniophora.olivacea                    ~ MetaData$TISSUE, OTUabund, sum)
aggregate(LREwh63.Coniophora.olivacea                    ~ MetaData$HOST, OTUabund, sum)


aggregate(LRE.wh34.Alternaria.quercicola                    ~ MetaData$SOIL, OTUabund, sum)
aggregate(LRE.wh34.Alternaria.quercicola                    ~ MetaData$SITE, OTUabund, sum)
aggregate(LRE.wh34.Alternaria.quercicola                    ~ MetaData$MEDIA, OTUabund, sum)
aggregate(LRE.wh34.Alternaria.quercicola                    ~ MetaData$TEMP, OTUabund, sum)
aggregate(LRE.wh34.Alternaria.quercicola                    ~ MetaData$TIME, OTUabund, sum)
aggregate(LRE.wh34.Alternaria.quercicola                    ~ MetaData$TISSUE, OTUabund, sum)
aggregate(LRE.wh34.Alternaria.quercicola                    ~ MetaData$HOST, OTUabund, sum)


aggregate(LDE.se7.Coniothyrium.aleuritis                    ~ MetaData$SOIL, OTUabund, sum)
aggregate(LDE.se7.Coniothyrium.aleuritis                    ~ MetaData$SITE, OTUabund, sum)
aggregate(LDE.se7.Coniothyrium.aleuritis                    ~ MetaData$MEDIA, OTUabund, sum)
aggregate(LDE.se7.Coniothyrium.aleuritis                    ~ MetaData$TEMP, OTUabund, sum)
aggregate(LDE.se7.Coniothyrium.aleuritis                    ~ MetaData$TIME, OTUabund, sum)
aggregate(LDE.se7.Coniothyrium.aleuritis                    ~ MetaData$TISSUE, OTUabund, sum)
aggregate(LDE.se7.Coniothyrium.aleuritis                    ~ MetaData$HOST, OTUabund, sum)


aggregate(LAE.se5.Sordaria.humana                     ~ MetaData$SOIL, OTUabund, sum)
aggregate(LAE.se5.Sordaria.humana                     ~ MetaData$SITE, OTUabund, sum)
aggregate(LAE.se5.Sordaria.humana                     ~ MetaData$MEDIA, OTUabund, sum)
aggregate(LAE.se5.Sordaria.humana                     ~ MetaData$TEMP, OTUabund, sum)
aggregate(LAE.se5.Sordaria.humana                     ~ MetaData$TIME, OTUabund, sum)
aggregate(LAE.se5.Sordaria.humana                     ~ MetaData$TISSUE, OTUabund, sum)
aggregate(LAE.se5.Sordaria.humana                     ~ MetaData$HOST, OTUabund, sum)

aggregate(HAE.se5.Camarosporomyces.flavigenus             ~ MetaData$SOIL, OTUabund, sum)
aggregate(HAE.se5.Camarosporomyces.flavigenus             ~ MetaData$SITE, OTUabund, sum)
aggregate(HAE.se5.Camarosporomyces.flavigenus             ~ MetaData$MEDIA, OTUabund, sum)
aggregate(HAE.se5.Camarosporomyces.flavigenus             ~ MetaData$TEMP, OTUabund, sum)
aggregate(HAE.se5.Camarosporomyces.flavigenus             ~ MetaData$TIME, OTUabund, sum)
aggregate(HAE.se5.Camarosporomyces.flavigenus             ~ MetaData$TISSUE, OTUabund, sum)
aggregate(HAE.se5.Camarosporomyces.flavigenus             ~ MetaData$HOST, OTUabund, sum)


aggregate(HAE.ss4.Cytospora.chrysosperma                ~ MetaData$SOIL, OTUabund, sum)
aggregate(HAE.ss4.Cytospora.chrysosperma                ~ MetaData$SITE, OTUabund, sum)
aggregate(HAE.ss4.Cytospora.chrysosperma                ~ MetaData$MEDIA, OTUabund, sum)
aggregate(HAE.ss4.Cytospora.chrysosperma                ~ MetaData$TEMP, OTUabund, sum)
aggregate(HAE.ss4.Cytospora.chrysosperma                ~ MetaData$TIME, OTUabund, sum)
aggregate(HAE.ss4.Cytospora.chrysosperma                ~ MetaData$TISSUE, OTUabund, sum)
aggregate(HAE.ss4.Cytospora.chrysosperma                ~ MetaData$HOST, OTUabund, sum)


aggregate(HAE.we5.Coniolariella.sp.                ~ MetaData$SOIL, OTUabund, sum)
aggregate(HAE.we5.Coniolariella.sp.                ~ MetaData$SITE, OTUabund, sum)
aggregate(HAE.we5.Coniolariella.sp.                ~ MetaData$MEDIA, OTUabund, sum)
aggregate(HAE.we5.Coniolariella.sp.                ~ MetaData$TEMP, OTUabund, sum)
aggregate(HAE.we5.Coniolariella.sp.                ~ MetaData$TIME, OTUabund, sum)
aggregate(HAE.we5.Coniolariella.sp.                ~ MetaData$TISSUE, OTUabund, sum)
aggregate(HAE.we5.Coniolariella.sp.                ~ MetaData$HOST, OTUabund, sum)


aggregate(HAE.wh26.Preussia.sp.                ~ MetaData$SOIL, OTUabund, sum)
aggregate(HAE.wh26.Preussia.sp.                ~ MetaData$SITE, OTUabund, sum)
aggregate(HAE.wh26.Preussia.sp.                ~ MetaData$MEDIA, OTUabund, sum)
aggregate(HAE.wh26.Preussia.sp.                ~ MetaData$TEMP, OTUabund, sum)
aggregate(HAE.wh26.Preussia.sp.                ~ MetaData$TIME, OTUabund, sum)
aggregate(HAE.wh26.Preussia.sp.                ~ MetaData$TISSUE, OTUabund, sum)
aggregate(HAE.wh26.Preussia.sp.                ~ MetaData$HOST, OTUabund, sum)



aggregate(HAE.wh10.Raffaelea.montetyi                ~ MetaData$SOIL, OTUabund, sum)
aggregate(HAE.wh10.Raffaelea.montetyi                ~ MetaData$SITE, OTUabund, sum)
aggregate(HAE.wh10.Raffaelea.montetyi                ~ MetaData$MEDIA, OTUabund, sum)
aggregate(HAE.wh10.Raffaelea.montetyi                ~ MetaData$TEMP, OTUabund, sum)
aggregate(HAE.wh10.Raffaelea.montetyi                ~ MetaData$TIME, OTUabund, sum)
aggregate(HAE.wh10.Raffaelea.montetyi                ~ MetaData$TISSUE, OTUabund, sum)
aggregate(HAE.wh10.Raffaelea.montetyi                ~ MetaData$HOST, OTUabund, sum)


aggregate(HAE.wh65.Coniophora.marmorata                ~ MetaData$SOIL, OTUabund, sum)
aggregate(HAE.wh65.Coniophora.marmorata                ~ MetaData$SITE, OTUabund, sum)
aggregate(HAE.wh65.Coniophora.marmorata                ~ MetaData$MEDIA, OTUabund, sum)
aggregate(HAE.wh65.Coniophora.marmorata                ~ MetaData$TEMP, OTUabund, sum)
aggregate(HAE.wh65.Coniophora.marmorata                ~ MetaData$TIME, OTUabund, sum)
aggregate(HAE.wh65.Coniophora.marmorata                ~ MetaData$TISSUE, OTUabund, sum)
aggregate(HAE.wh65.Coniophora.marmorata                ~ MetaData$HOST, OTUabund, sum)


aggregate(HAE.se9.Chaetomium.nigricolor        ~ MetaData$SOIL, OTUabund, sum)
aggregate(HAE.se9.Chaetomium.nigricolor        ~ MetaData$SITE, OTUabund, sum)
aggregate(HAE.se9.Chaetomium.nigricolor        ~ MetaData$MEDIA, OTUabund, sum)
aggregate(HAE.se9.Chaetomium.nigricolor        ~ MetaData$TEMP, OTUabund, sum)
aggregate(HAE.se9.Chaetomium.nigricolor        ~ MetaData$TIME, OTUabund, sum)
aggregate(HAE.se9.Chaetomium.nigricolor        ~ MetaData$TISSUE, OTUabund, sum)
aggregate(HAE.se9.Chaetomium.nigricolor        ~ MetaData$HOST, OTUabund, sum)


aggregate(HAE.se1.Acrocalymma.sp.         ~ MetaData$SOIL, OTUabund, sum)
aggregate(HAE.se1.Acrocalymma.sp.         ~ MetaData$SITE, OTUabund, sum)
aggregate(HAE.se1.Acrocalymma.sp.         ~ MetaData$MEDIA, OTUabund, sum)
aggregate(HAE.se1.Acrocalymma.sp.         ~ MetaData$TEMP, OTUabund, sum)
aggregate(HAE.se1.Acrocalymma.sp.         ~ MetaData$TIME, OTUabund, sum)
aggregate(HAE.se1.Acrocalymma.sp.         ~ MetaData$TISSUE, OTUabund, sum)
aggregate(HAE.se1.Acrocalymma.sp.         ~ MetaData$HOST, OTUabund, sum)

aggregate(HBE.wh23.Kalmusia.spartii         ~ MetaData$SOIL, OTUabund, sum)
aggregate(HBE.wh23.Kalmusia.spartii         ~ MetaData$SITE, OTUabund, sum)
aggregate(HBE.wh23.Kalmusia.spartii         ~ MetaData$MEDIA, OTUabund, sum)
aggregate(HBE.wh23.Kalmusia.spartii         ~ MetaData$TEMP, OTUabund, sum)
aggregate(HBE.wh23.Kalmusia.spartii         ~ MetaData$TIME, OTUabund, sum)
aggregate(HBE.wh23.Kalmusia.spartii         ~ MetaData$TISSUE, OTUabund, sum)
aggregate(HBE.wh23.Kalmusia.spartii         ~ MetaData$HOST, OTUabund, sum)



aggregate(HBE.wh27.Chaetomium.truncatulum         ~ MetaData$SOIL, OTUabund, sum)
aggregate(HBE.wh27.Chaetomium.truncatulum         ~ MetaData$SITE, OTUabund, sum)
aggregate(HBE.wh27.Chaetomium.truncatulum         ~ MetaData$MEDIA, OTUabund, sum)
aggregate(HBE.wh27.Chaetomium.truncatulum         ~ MetaData$TEMP, OTUabund, sum)
aggregate(HBE.wh27.Chaetomium.truncatulum         ~ MetaData$TIME, OTUabund, sum)
aggregate(HBE.wh27.Chaetomium.truncatulum         ~ MetaData$TISSUE, OTUabund, sum)
aggregate(HBE.wh27.Chaetomium.truncatulum         ~ MetaData$HOST, OTUabund, sum)

aggregate(GME.ws1.Preussia.Africana           ~ MetaData$SOIL, OTUabund, sum)
aggregate(GME.ws1.Preussia.Africana           ~ MetaData$SITE, OTUabund, sum)
aggregate(GME.ws1.Preussia.Africana           ~ MetaData$MEDIA, OTUabund, sum)
aggregate(GME.ws1.Preussia.Africana           ~ MetaData$TEMP, OTUabund, sum)
aggregate(GME.ws1.Preussia.Africana           ~ MetaData$TIME, OTUabund, sum)
aggregate(GME.ws1.Preussia.Africana           ~ MetaData$TISSUE, OTUabund, sum)
aggregate(GME.ws1.Preussia.Africana           ~ MetaData$HOST, OTUabund, sum)


aggregate(GME.ss5.Aspergillus.frequens           ~ MetaData$SOIL, OTUabund, sum)
aggregate(GME.ss5.Aspergillus.frequens           ~ MetaData$SITE, OTUabund, sum)
aggregate(GME.ss5.Aspergillus.frequens           ~ MetaData$MEDIA, OTUabund, sum)
aggregate(GME.ss5.Aspergillus.frequens           ~ MetaData$TEMP, OTUabund, sum)
aggregate(GME.ss5.Aspergillus.frequens           ~ MetaData$TIME, OTUabund, sum)
aggregate(GME.ss5.Aspergillus.frequens           ~ MetaData$TISSUE, OTUabund, sum)
aggregate(GME.ss5.Aspergillus.frequens           ~ MetaData$HOST, OTUabund, sum)


aggregate(GAE.ws1.Tricharina.striispora            ~ MetaData$SOIL, OTUabund, sum)
aggregate(GAE.ws1.Tricharina.striispora            ~ MetaData$SITE, OTUabund, sum)
aggregate(GAE.ws1.Tricharina.striispora            ~ MetaData$MEDIA, OTUabund, sum)
aggregate(GAE.ws1.Tricharina.striispora            ~ MetaData$TEMP, OTUabund, sum)
aggregate(GAE.ws1.Tricharina.striispora            ~ MetaData$TIME, OTUabund, sum)
aggregate(GAE.ws1.Tricharina.striispora            ~ MetaData$TISSUE, OTUabund, sum)
aggregate(GAE.ws1.Tricharina.striispora            ~ MetaData$HOST, OTUabund, sum)


aggregate(GAE.ss.6.Fusarium.ensiforme            ~ MetaData$SOIL, OTUabund, sum)
aggregate(GAE.ss.6.Fusarium.ensiforme            ~ MetaData$SITE, OTUabund, sum)
aggregate(GAE.ss.6.Fusarium.ensiforme            ~ MetaData$MEDIA, OTUabund, sum)
aggregate(GAE.ss.6.Fusarium.ensiforme            ~ MetaData$TEMP, OTUabund, sum)
aggregate(GAE.ss.6.Fusarium.ensiforme            ~ MetaData$TIME, OTUabund, sum)
aggregate(GAE.ss.6.Fusarium.ensiforme            ~ MetaData$TISSUE, OTUabund, sum)
aggregate(GAE.ss.6.Fusarium.ensiforme            ~ MetaData$HOST, OTUabund, sum)


aggregate(GAE.ss3.Aspergillus.jensenii              ~ MetaData$SOIL, OTUabund, sum)
aggregate(GAE.ss3.Aspergillus.jensenii              ~ MetaData$SITE, OTUabund, sum)
aggregate(GAE.ss3.Aspergillus.jensenii              ~ MetaData$MEDIA, OTUabund, sum)
aggregate(GAE.ss3.Aspergillus.jensenii              ~ MetaData$TEMP, OTUabund, sum)
aggregate(GAE.ss3.Aspergillus.jensenii              ~ MetaData$TIME, OTUabund, sum)
aggregate(GAE.ss3.Aspergillus.jensenii              ~ MetaData$TISSUE, OTUabund, sum)
aggregate(GAE.ss3.Aspergillus.jensenii              ~ MetaData$HOST, OTUabund, sum)


aggregate(GAE.ss7.Aspergillus.conversis               ~ MetaData$SOIL, OTUabund, sum)
aggregate(GAE.ss7.Aspergillus.conversis               ~ MetaData$SITE, OTUabund, sum)
aggregate(GAE.ss7.Aspergillus.conversis              ~ MetaData$MEDIA, OTUabund, sum)
aggregate(GAE.ss7.Aspergillus.conversis              ~ MetaData$TEMP, OTUabund, sum)
aggregate(GAE.ss7.Aspergillus.conversis              ~ MetaData$TIME, OTUabund, sum)
aggregate(GAE.ss7.Aspergillus.conversis              ~ MetaData$TISSUE, OTUabund, sum)
aggregate(GAE.ss7.Aspergillus.conversis              ~ MetaData$HOST, OTUabund, sum)


aggregate(GAE.ws8.Sporormiella.septenaria            ~ MetaData$SOIL, OTUabund, sum)
aggregate(GAE.ws8.Sporormiella.septenaria            ~ MetaData$SITE, OTUabund, sum)
aggregate(GAE.ws8.Sporormiella.septenaria            ~ MetaData$MEDIA, OTUabund, sum)
aggregate(GAE.ws8.Sporormiella.septenaria            ~ MetaData$TEMP, OTUabund, sum)
aggregate(GAE.ws8.Sporormiella.septenaria            ~ MetaData$TIME, OTUabund, sum)
aggregate(GAE.ws8.Sporormiella.septenaria            ~ MetaData$TISSUE, OTUabund, sum)
aggregate(GAE.ws8.Sporormiella.septenaria            ~ MetaData$HOST, OTUabund, sum)


aggregate(GspE.ss6.Talaromyces.cecidicola             ~ MetaData$SOIL, OTUabund, sum)
aggregate(GspE.ss6.Talaromyces.cecidicola             ~ MetaData$SITE, OTUabund, sum)
aggregate(GspE.ss6.Talaromyces.cecidicola             ~ MetaData$MEDIA, OTUabund, sum)
aggregate(GspE.ss6.Talaromyces.cecidicola             ~ MetaData$TEMP, OTUabund, sum)
aggregate(GspE.ss6.Talaromyces.cecidicola             ~ MetaData$TIME, OTUabund, sum)
aggregate(GspE.ss6.Talaromyces.cecidicola             ~ MetaData$TISSUE, OTUabund, sum)
aggregate(GspE.ss6.Talaromyces.cecidicola             ~ MetaData$HOST, OTUabund, sum)


aggregate(GspE.ss5.Talaromyces.leycettanus             ~ MetaData$SOIL, OTUabund, sum)
aggregate(GspE.ss5.Talaromyces.leycettanus             ~ MetaData$SITE, OTUabund, sum)
aggregate(GspE.ss5.Talaromyces.leycettanus             ~ MetaData$MEDIA, OTUabund, sum)
aggregate(GspE.ss5.Talaromyces.leycettanus             ~ MetaData$TEMP, OTUabund, sum)
aggregate(GspE.ss5.Talaromyces.leycettanus             ~ MetaData$TIME, OTUabund, sum)
aggregate(GspE.ss5.Talaromyces.leycettanus             ~ MetaData$TISSUE, OTUabund, sum)
aggregate(GspE.ss5.Talaromyces.leycettanus             ~ MetaData$HOST, OTUabund, sum)


aggregate(GAE.ws4.Stemphylium.vesicarium               ~ MetaData$SOIL, OTUabund, sum)
aggregate(GAE.ws4.Stemphylium.vesicarium               ~ MetaData$SITE, OTUabund, sum)
aggregate(GAE.ws4.Stemphylium.vesicarium               ~ MetaData$MEDIA, OTUabund, sum)
aggregate(GAE.ws4.Stemphylium.vesicarium               ~ MetaData$TEMP, OTUabund, sum)
aggregate(GAE.ws4.Stemphylium.vesicarium               ~ MetaData$TIME, OTUabund, sum)
aggregate(GAE.ws4.Stemphylium.vesicarium               ~ MetaData$TISSUE, OTUabund, sum)
aggregate(GAE.ws4.Stemphylium.vesicarium               ~ MetaData$HOST, OTUabund, sum)



aggregate(ECE.ws8.Sporormiella.dakotensis              ~ MetaData$SOIL, OTUabund, sum)
aggregate(ECE.ws8.Sporormiella.dakotensis              ~ MetaData$SITE, OTUabund, sum)
aggregate(ECE.ws8.Sporormiella.dakotensis              ~ MetaData$MEDIA, OTUabund, sum)
aggregate(ECE.ws8.Sporormiella.dakotensis              ~ MetaData$TEMP, OTUabund, sum)
aggregate(ECE.ws8.Sporormiella.dakotensis              ~ MetaData$TIME, OTUabund, sum)
aggregate(ECE.ws8.Sporormiella.dakotensis              ~ MetaData$TISSUE, OTUabund, sum)
aggregate(ECE.ws8.Sporormiella.dakotensis              ~ MetaData$HOST, OTUabund, sum)

aggregate(ECE.ws8.Sporormiella.dakotensis              ~ MetaData$SOIL, OTUabund, sum)
aggregate(ECE.ws8.Sporormiella.dakotensis              ~ MetaData$SITE, OTUabund, sum)
aggregate(ECE.ws8.Sporormiella.dakotensis              ~ MetaData$MEDIA, OTUabund, sum)
aggregate(ECE.ws8.Sporormiella.dakotensis              ~ MetaData$TEMP, OTUabund, sum)
aggregate(ECE.ws8.Sporormiella.dakotensis              ~ MetaData$TIME, OTUabund, sum)
aggregate(ECE.ws8.Sporormiella.dakotensis              ~ MetaData$TISSUE, OTUabund, sum)
aggregate(ECE.ws8.Sporormiella.dakotensis              ~ MetaData$HOST, OTUabund, sum)

aggregate(ECEws7.Nigrospora.oryzae              ~ MetaData$SOIL, OTUabund, sum)
aggregate(ECEws7.Nigrospora.oryzae              ~ MetaData$SITE, OTUabund, sum)
aggregate(ECEws7.Nigrospora.oryzae              ~ MetaData$MEDIA, OTUabund, sum)
aggregate(ECEws7.Nigrospora.oryzae              ~ MetaData$TEMP, OTUabund, sum)
aggregate(ECEws7.Nigrospora.oryzae              ~ MetaData$TIME, OTUabund, sum)
aggregate(ECEws7.Nigrospora.oryzae              ~ MetaData$TISSUE, OTUabund, sum)
aggregate(ECEws7.Nigrospora.oryzae              ~ MetaData$HOST, OTUabund, sum)

aggregate(ECE..ss10.Meyerozyma.guilliermondii              ~ MetaData$SOIL, OTUabund, sum)
aggregate(ECE..ss10.Meyerozyma.guilliermondii              ~ MetaData$SITE, OTUabund, sum)
aggregate(ECE..ss10.Meyerozyma.guilliermondii              ~ MetaData$MEDIA, OTUabund, sum)
aggregate(ECE..ss10.Meyerozyma.guilliermondii              ~ MetaData$TEMP, OTUabund, sum)
aggregate(ECE..ss10.Meyerozyma.guilliermondii              ~ MetaData$TIME, OTUabund, sum)
aggregate(ECE..ss10.Meyerozyma.guilliermondii              ~ MetaData$TISSUE, OTUabund, sum)
aggregate(ECE..ss10.Meyerozyma.guilliermondii              ~ MetaData$HOST, OTUabund, sum)

aggregate(ERE.ws2.Fusarium.sp.               ~ MetaData$SOIL, OTUabund, sum)
aggregate(ERE.ws2.Fusarium.sp.               ~ MetaData$SITE, OTUabund, sum)
aggregate(ERE.ws2.Fusarium.sp.               ~ MetaData$MEDIA, OTUabund, sum)
aggregate(ERE.ws2.Fusarium.sp.               ~ MetaData$TEMP, OTUabund, sum)
aggregate(ERE.ws2.Fusarium.sp.               ~ MetaData$TIME, OTUabund, sum)
aggregate(ERE.ws2.Fusarium.sp.               ~ MetaData$TISSUE, OTUabund, sum)
aggregate(ERE.ws2.Fusarium.sp.               ~ MetaData$HOST, OTUabund, sum)

aggregate(ERE.ws1.Preussia.intermedia               ~ MetaData$SOIL, OTUabund, sum)
aggregate(ERE.ws1.Preussia.intermedia               ~ MetaData$SITE, OTUabund, sum)
aggregate(ERE.ws1.Preussia.intermedia               ~ MetaData$MEDIA, OTUabund, sum)
aggregate(ERE.ws1.Preussia.intermedia               ~ MetaData$TEMP, OTUabund, sum)
aggregate(ERE.ws1.Preussia.intermedia               ~ MetaData$TIME, OTUabund, sum)
aggregate(ERE.ws1.Preussia.intermedia               ~ MetaData$TISSUE, OTUabund, sum)
aggregate(ERE.ws1.Preussia.intermedia               ~ MetaData$HOST, OTUabund, sum)

aggregate(ERE.ws4.Talaromyces.pinophilus                 ~ MetaData$SOIL, OTUabund, sum)
aggregate(ERE.ws4.Talaromyces.pinophilus                 ~ MetaData$SITE, OTUabund, sum)
aggregate(ERE.ws4.Talaromyces.pinophilus                 ~ MetaData$MEDIA, OTUabund, sum)
aggregate(ERE.ws4.Talaromyces.pinophilus                 ~ MetaData$TEMP, OTUabund, sum)
aggregate(ERE.ws4.Talaromyces.pinophilus                 ~ MetaData$TIME, OTUabund, sum)
aggregate(ERE.ws4.Talaromyces.pinophilus                 ~ MetaData$TISSUE, OTUabund, sum)
aggregate(ERE.ws4.Talaromyces.pinophilus                 ~ MetaData$HOST, OTUabund, sum)


aggregate(ESE.wh29.Schizothecium.inaequale       ~ MetaData$SOIL, OTUabund, sum)
aggregate(ESE.wh29.Schizothecium.inaequale       ~ MetaData$SITE, OTUabund, sum)
aggregate(ESE.wh29.Schizothecium.inaequale       ~ MetaData$MEDIA, OTUabund, sum)
aggregate(ESE.wh29.Schizothecium.inaequale       ~ MetaData$TEMP, OTUabund, sum)
aggregate(ESE.wh29.Schizothecium.inaequale       ~ MetaData$TIME, OTUabund, sum)
aggregate(ESE.wh29.Schizothecium.inaequale       ~ MetaData$TISSUE, OTUabund, sum)
aggregate(ESE.wh29.Schizothecium.inaequale       ~ MetaData$HOST, OTUabund, sum)




aggregate(. ~ MetaData$SOIL, OTUabund, sum)

OTU.slic<- c(866,352,303,589,434,348,535,877,500,394,301, 4603)  #get the valus from OTU.colsum

OTU.lbls<- c(expression(italic("Rosellinia sp.")),expression(italic("Acrocalymma vagum")),
             expression(italic("Dimorphosporicola tragani")),
             expression(italic("Raffaelea montetyi")),expression(italic("Paracamarosporium hawaiiense")),
             expression(italic("Fusariella sinensis")),expression(italic("Humicola fuscoatra")),
             expression(italic("Neocamarosporium chichastianum")),
             expression(italic("Camarosporomyces flavigenus")),
             expression(italic("Preussia sp.")),expression(italic("Coniophora marmorata")),"Other")


OTU.Percent<-round(OTU.slic/sum (OTU.slic)*100, digits=2)
OTU.lbls <- paste(OTU.lbls, OTU.Percent)
order.lbls<-paste(OTU.lbls,"%",sep="")
pie(OTU.slic,labels =OTU.lbls, col = c("red","skyblue1","magenta",
                                       "deeppink1","mediumblue","royalblue1","orchid1","cyan",
                                       "yellow", "springgreen2", "pink","green" ) , main = "OTU", 
    cex=1,border = NA,cex.main= 0.5, radius = 1)

OTU.Percent<-round(OTU.slic/sum (OTU.slic)*100, digits=2)
OTU.lbls<- c (expression (italic("Rosellinia")*' sp. (8.57%)' ), expression (italic("Acrocalymma vagum")*' (3.48%)'), expression(italic("Dimorphosporicola")*' sp. (3%)'),
              expression(italic("Raffaelea montetyi")*' (5.83%)'), expression(italic("Paracamarosporium")*' sp. (4.3%)'), expression(italic("Fusariella")*' sp. (3.44%)'),
              expression(italic("Humicola fuscoatra")*' (5.3%)'),expression(italic("Neocamarosporium chichastianum")*' (8.68%)'),expression(italic("Chaetosphaeronema")*' sp. (4.95%)'),
              expression(italic("Preussia")*' sp. (3.9%)'),expression (italic("Coniophora")*' sp. (2.98%)'), "Other (45.57%)")

pie (OTU.slic,labels =OTU.lbls, col = c("red","skyblue1","magenta",
                                        "deeppink1","mediumblue","royalblue1","orchid1","cyan",
                                        "yellow", "springgreen2", "pink","green" ) , main = "OTUs frequency", cex=0.8,border = NA,cex.main= 1.1, radius =0.85)


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

#######################################
#IR model
#######################################
hist(Aridsoil$IR)
IR.m<-glm(formula =IR~HOST+TIME+SITE,data = Aridsoil,
          family=poisson(link = "log"))



####### Subsetting the data for Gypsumsoil
Gypsumsoil = subset (MetaData, SOIL%in%c("Gypsum Soil"))
GypsumsoilOTU = subset (OTUabund, MetaData$SOIL %in% c("Gypsum Soil"))

OTUsINGypsumsoil<-colSums(GypsumsoilOTU)
GypsumsoilOTU<-GypsumsoilOTU[, colSums(GypsumsoilOTU != 0) > 0]
rownames(GypsumsoilOTU)==rownames(Gypsumsoil)
colSums(GypsumsoilOTU)
colnames(GypsumsoilOTU)

####### Subsetting the data for Salinesoil
Salinesoil = subset (MetaData, SOIL%in%c("Saline Soil"))
SalinesoilOTU = subset (OTUabund, MetaData$SOIL %in% c("Saline Soil"))

OTUsINSalinesoil<-colSums(SalinesoilOTU)
SalinesoilOTU<-SalinesoilOTU[, colSums(SalinesoilOTU != 0) > 0]
rownames(SalinesoilOTU)==rownames(Salinesoil)
colSums(SalinesoilOTU)
colnames(SalinesoilOTU)









