
#### Download all files associated to this script (REPHYTOX data (CSV files starting with Q2_), regulatory thresholds table (DFreglem2012-2019), fishing season table (Saison_de_peche)  
#### Import data on phycotoxin concentration in shellfish (REPHYTOX filtered data: data used in regulation)
#### Chose your working directory using the function setwd()
DataToxins <- rbind.data.frame(read.csv("Q2_210224_60255817_REPHYTOX-AC-Tox Lip-pour-REB-EM_2010-2015.csv",header = TRUE, sep = ";", quote = "\"", dec = ".", fill = TRUE, comment.char = ""),
                               read.csv("Q2_210224_60255835_REPHYTOX-AC-Tox Lip-pour-REB-EM_2016-2020.csv",header = TRUE, sep = ";", quote = "\"", dec = ".", fill = TRUE, comment.char = ""),
                               read.csv("Q2_210224_60255849_REPHYTOX-ASP-pour-REB-EM_1999-2020.csv",header = TRUE, sep = ";", quote = "\"", dec = ".", fill = TRUE, comment.char = ""),
                               read.csv("Q2_210224_60255858_REPHYTOX-PSP-pour-REB-EM_1988-2020.csv",header = TRUE, sep = ";", quote = "\"", dec = ".", fill = TRUE, comment.char = ""))

head(DataToxins)
dim(DataToxins)
DataToxins <- unique(DataToxins) 
dim(DataToxins)

#### Import the "DFreglem2012-2019" table of regulatory thresholds
DFref <- read.table("DFreglem2012-2019.txt",sep="\t",header=TRUE,quote="")
DFref <- unique(DFref) #DFref is the regulatory thresholds table
dim(DFref)
DFref

#### Import the fishing season table
DFsaisons <- read.table("Saison_de_peche.txt",sep="\t",header=TRUE,quote="")
head(DFsaisons)
dim(DFsaisons)


#### Select data concerning the scallop of the Eastern Channel area
library(dplyr)

DataToxinsCSJ <- subset(DataToxins, Lieu.de.surveillance...Entité.de.classement...Identifiant%in%(3:16) &
                          Passage...Année>=2011 &
                          Echantillon...Libellé.du.taxon.support%in%"Pecten maximus")

#### Create the variable "FishingZone" from the column "Lieu.de.surveillance...Libellé"
DataToxinsCSJ$FishingZone <- sapply(as.character(DataToxinsCSJ$Lieu.de.surveillance...Libellé),function(x) strsplit(x,"- ")[[1]][2])
DataToxinsCSJ$FishingZone <- gsub("Zone ","",DataToxinsCSJ$FishingZone)
DataToxinsCSJ$FishingZone[DataToxinsCSJ$Lieu.de.surveillance...Libellé%in%"Nord Cotentin"] <- "Nord Cotentin"
DataToxinsCSJ$FishingZone[DataToxinsCSJ$Lieu.de.surveillance...Libellé%in%"Extérieur gisement baie de Seine"] <- "I"
table(as.character(DataToxinsCSJ$Lieu.de.surveillance...Libellé))
table(DataToxinsCSJ$FishingZone)


DF <- merge(DataToxinsCSJ,DFref,all.x=TRUE)
dim(DF)
head(DF)

unique(DF[!is.na(DF$Groupe.de.toxines),1:3])
unique(DF[is.na(DF$Groupe.de.toxines),1:3])

# We keep only occurrences included in the reference table
DF <- DF[!is.na(DF$Groupe.de.toxines),]  
dim(DF)

##### Create the variable "Toxins" according to regulations (official data)

DF$Toxins <- ""
table(DF$Toxins)
DF[DF$Résultat...Code.paramétre%in%"ASP" & DF$Résultat...Libellé.fraction%in%"Chair totale égouttée","Toxins"] <- "ASP"
table(DF$Toxins)
DF[DF$Résultat...Code.paramétre%in%"ASP" & DF$Résultat...Libellé.fraction%in%c("Muscle","Muscle + gonade"),"Toxins"] <- "ASPshell"
table(DF$Toxins)
DF[DF$Résultat...Code.paramétre%in%c("AO+DTXs+PTXs-TEFs","AZAs-TEFs","YTXs-TEFs"),"Toxins"] <- "DSP"
table(DF$Toxins)
DF[DF$Résultat...Code.paramétre%in%"TOXPSP","Toxins"] <- "PSP"
table(DF$Toxins)

#### Create the variable "Toxicity": if the concentration value is under the regulatory threshold then Toxicity returns "no" otherwise it returns "yes"
DF$Toxicity <- "no"
DF$Toxicity[DF$Résultat...Valeur.de.la.mesure > DF$Seuil.réglementaire] <- "yes"
table(DF$Toxicity)

#### Create the variable "Alert": binary coding of the variable "Toxicity"
DF$Alert <- NA
table(DF$Toxins,DF$Alert)
DF$Alert[DF$Toxins%in%"DSP" & DF$Toxicity%in%"yes"] <- 1
DF$Alert[DF$Toxins%in%"DSP" & DF$Toxicity%in%"no"] <- 0
table(DF$Toxins,DF$Alert)
DF$Alert[DF$Toxins%in%"PSP" & DF$Toxicity%in%"yes"] <- 1
DF$Alert[DF$Toxins%in%"PSP" & DF$Toxicity%in%"no"] <- 0
table(DF$Toxins,DF$Alert)
DF$Alert[DF$Toxins%in%"ASP" & DF$Résultat...Valeur.de.la.mesure>250] <- 1
table(DF$Toxins,DF$Alert)

#### Definition of the secondary component of "Toxicity" for ASP: "ASPshell"
subASP2 <- subset(DF,Toxins%in%"ASPshell")
rownames(subASP2) <- 1:nrow(subASP2)
DF <- subset(DF,!Toxins%in%"ASPshell")
subASP3 <- tapply(1:nrow(subASP2),list(subASP2$Echantillon...Identifiant.interne),function(x) ifelse(any(subASP2[x,"Toxicity"]%in%"yes"),"yes","no"))
DF$ToxicityBis <- subASP3[as.character(DF$Echantillon...Identifiant.interne)]
 
DF$ToxicityBis[!DF$Toxins%in%"ASP"] <- NA

DF$Alert[DF$Toxins%in%"ASP" & DF$Toxicity%in%"yes" & DF$Résultat...Valeur.de.la.mesure<=250 & DF$ToxicityBis%in%"yes"] <- 1
table(DF$Toxins,DF$Alert)
DF$Alert[DF$Toxins%in%"ASP" & DF$Toxicity%in%"yes" & DF$Résultat...Valeur.de.la.mesure<=250 & DF$ToxicityBis%in%"no"] <- 0
table(DF$Toxins,DF$Alert)
DF$Alert[DF$Toxins%in%"ASP" & DF$Toxicity%in%"yes" & is.na(DF$ToxicityBis)] <- 1
table(DF$Toxins,DF$Alert)
DF$Alert[DF$Toxins%in%"ASP" & DF$Toxicity%in%"no"] <- 0
table(DF$Toxins,DF$Alert)

DFcopy <- DF

#### Compile the data to have only one line for each combination FishingZone-Date and put the "Alert" data per toxin familly in column  
DF$Passage...Date <- as.character(DF$Passage...Date) ; DF$FishingZone <- as.character(DF$FishingZone)
MAT <- with(DF,tapply(Alert,list(paste(FishingZone,Passage...Date,sep=":-:"),Toxins),function(x) ifelse(all(is.na(x)),NA,max(x,na.rm=TRUE))))

nam <- c("Lieu.de.surveillance...Entité.de.classement...Identifiant",
         "Lieu.de.surveillance...Entité.de.classement...Libellé",
         "Passage...Date","Passage...Jour","Passage...Mois","Passage...Année",
         "Echantillon...Libellé.du.taxon.support","FishingZone")

# Verification
dim(unique(DF[,nam]))
dim(MAT)


DFalt <- unique(DF[,nam])
DF <- cbind.data.frame(DFalt,MAT[match(paste(DFalt$FishingZone,DFalt$Passage...Date,sep=":-:"),rownames(MAT)),])

#### Create the variable "AlertTot"
DF$AlertTot <- apply(DF[,c("ASP","DSP","PSP")],1,function(x) as.numeric(any(x[!is.na(x)]%in%1)))

#### Convert the columns "date"
DF$Date_Passage_V2 <- as.Date(as.character(DF$Passage...Date), "%d/%m/%Y")
DF$Date <- as.character(DF$Date_Passage_V2)

#### Develop the intermediate days between sampling dates (the dimensional support is (taxon x zone))
support <- unique(DF[,c("Lieu.de.surveillance...Entité.de.classement...Identifiant","Lieu.de.surveillance...Entité.de.classement...Libellé",
                        "FishingZone","Echantillon...Libellé.du.taxon.support")])

#### Develop the period on a daily basis
allDays <- seq(min(DF$Date_Passage_V2),max(DF$Date_Passage_V2),by=1)
FullDF <- cbind.data.frame(support,Date=rep(allDays,each=nrow(support)))
dim(FullDF)
#### Integrate the table: DF (FullDF is the daily table)
FullDF <- merge(FullDF,DF,all.x=TRUE)
                     
#### Associate the number of week to each date (we start from the reference: Monday 03 January 2011 which is the first Monday in 2011, we add 0.02 to cover changes in hours)
FullDF$Date_V2 <- as.Date(as.character(FullDF$Date))
FullDF$CumulW <- difftime(FullDF$Date_V2,strptime("03-01-2011", format = "%d-%m-%Y"),units="weeks")  
FullDF$CumulW <- as.numeric(as.character(floor(FullDF$CumulW+0.02)))+1     

#### Attach the developed table of fishing season over the considered periods 
DFsaisons$Date1 <- as.Date(as.character(DFsaisons$Début), "%d/%m/%Y")
DFsaisons$Date2 <- as.Date(as.character(DFsaisons$Fin), "%d/%m/%Y")
DFsaisons$FishingZone <- DFsaisons$Zone
DFsaisonsDev <- do.call("rbind.data.frame",lapply(1:nrow(DFsaisons),function(x)
  cbind.data.frame(DFsaisons[x,c("FishingZone","Saison")],Date_V2=seq(DFsaisons$Date1[x],DFsaisons$Date2[x],by=1))))
dim(FullDF)
FullDF <- merge(FullDF,DFsaisonsDev,all.x=TRUE)
dim(FullDF)

#### Determine the variable "FishingSeason"
FullDF$FishingSeason <- "no"
FullDF$FishingSeason[!is.na(FullDF$Saison)] <- "yes"

#### Operate the closure assessment considering independently each Fishing Zone
FullDF_V2 <- do.call("rbind.data.frame",lapply(unique(FullDF$FishingZone), function(x) {
  tmp <- subset(FullDF,FishingZone%in%x)
  tmp <- tmp[order(tmp$Date_V2),]
  # Create the vector "Closure" day after day using 'tmp'
  # We suppose that the area was not closed in the beginning of the analysed period 
  CL <- rep(NA,nrow(tmp)) ; CL[1] <- 0 ; countWeeks <- NA ; histTest <- "" 
  for (i in 2:nrow(tmp)) {
    # 3 options
    if (is.na(tmp$AlertTot[i])) CL[i] <- CL[i-1] #nothing changes
    if (tmp$AlertTot[i]%in%"1") {CL[i] <- 1 ; histTest <- "+" ; countWeeks <- tmp$CumulW[i]} 
    if (tmp$AlertTot[i]%in%"0") {
      t1 <- (histTest%in%"") ; t2 <- (histTest%in%"+") ; t3 <- (histTest%in%"+-")
      if (t1) CL[i] <- 0
      if (t2) {histTest <- "+-" ; CL[i] <- 1}
      if (t3) {histTest <- "" ; CL[i] <- 0}   #the area is reopened ans the "histTest" restarts 
    }
  }
  tmp$Closure <- CL
  # Define the "Closure" events: A closure is defined as a 1 preceded by a 0 in the "Closure" column for the "FishingZone" variable (or a 1 at the beginning of the time series)
  tmp$ClosureEvent <- 0
  tmp$ClosureEvent[grep("01",paste(tmp$Closure[-nrow(tmp)],tmp$Closure[-1],sep=""))+1] <- 1
  # Manage the first value if equal to 1 
  if (tmp$Closure[1]%in%1) tmp$ClosureEvent[1] <- 1
  # Add distinctions by toxin family 
  tmp$ClosureEventDSP <- tmp$ClosureEventASP <- tmp$ClosureEventPSP <- 0
  tmp$ClosureEventDSP[tmp$ClosureEvent%in%1 & tmp$DSP%in%1] <- 1
  tmp$ClosureEventASP[tmp$ClosureEvent%in%1 & tmp$ASP%in%1] <- 1
  tmp$ClosureEventPSP[tmp$ClosureEvent%in%1 & tmp$PSP%in%1] <- 1
  tmp$ClosureDSP <- tmp$ClosureEventDSP ; tmp$ClosureASP <- tmp$ClosureEventASP ; tmp$ClosurePSP <- tmp$ClosureEventPSP
  # Browse according to "Closure"
  for (i in 2:nrow(tmp)) {
    if  (tmp$ClosureDSP[i-1]%in%1 & tmp$Closure[i]%in%1 & !(tmp$DSP[i]%in%0 & any(c(tmp$ASP[i],tmp$PSP[i])%in%1))) tmp$ClosureDSP[i] <- 1     #si un test s'est avéré négatif au méme moment mais que la fermeture est maintenue é cause d'un autre test, il ne faut plus associer cette fermeture é cette toxine
    if  (tmp$ClosureASP[i-1]%in%1 & tmp$Closure[i]%in%1 & !(tmp$ASP[i]%in%0 & any(c(tmp$DSP[i],tmp$PSP[i])%in%1))) tmp$ClosureASP[i] <- 1     #si un test s'est avéré négatif au méme moment mais que la fermeture est maintenue é cause d'un autre test, il ne faut plus associer cette fermeture é cette toxine
    if  (tmp$ClosurePSP[i-1]%in%1 & tmp$Closure[i]%in%1 & !(tmp$PSP[i]%in%0 & any(c(tmp$ASP[i],tmp$DSP[i])%in%1))) tmp$ClosurePSP[i] <- 1     #si un test s'est avéré négatif au méme moment mais que la fermeture est maintenue é cause d'un autre test, il ne faut plus associer cette fermeture é cette toxine
  }
  
  return(tmp)
}))

FullDF_V2$jour <- as.numeric(as.character(strftime(FullDF_V2$Date_V2,format="%d")))
FullDF_V2$mois <- as.numeric(as.character(strftime(FullDF_V2$Date_V2,format="%m")))
FullDF_V2$année <- strftime(FullDF_V2$Date_V2,format="%Y")

#### Export the complete table (closure calendar)
write.table(FullDF_V2,file="ClosureCalendar.txt",sep="\t",quote=FALSE,row.names=FALSE)

#### Synthesis of the developed table, in terms of the number and the duration of closures (in days) per area and per week and then per fishing season
   # Filtering the data on the fishing seasons
FullDF_Open <- subset(FullDF_V2,!is.na(Saison))
FullDF_Open$Zone <- factor(as.character(FullDF_Open$FishingZone),levels=c(as.character(1:15),"I","J","Nord Cotentin"))
   # Theoretical closures per week
TheoClosures <- with(FullDF_Open,aggregate(list(ClosureNbrDSP=ClosureEventDSP,ClosureNbrASP=ClosureEventASP,ClosureNbrPSP=ClosureEventPSP,
                                                ClosureDurationDSP=ClosureDSP,ClosureDurationASP=ClosureASP,ClosureDurationPSP=ClosurePSP,
                                                ClosureNbrTot=ClosureEvent,ClosureDurationTot=Closure),list(FishingZone=Zone,Week=CumulW,Season=Saison),sum))

   # Add the reference "Week" (according to the studied period --> weeks out of the fishing season are not shown)
weeksRef <- merge(with(FullDF_Open,aggregate(list(FirstDay=Date),list(Week=CumulW),min)),
                  with(FullDF_Open,aggregate(list(LastDay=Date),list(Week=CumulW),max)),all=TRUE)

dim(TheoClosures)
TheoClosures <- merge(TheoClosures,weeksRef,all.x=TRUE)
dim(TheoClosures)

   # Theoretical closure per fishing season 
TheoClosuresSeasons <- with(FullDF_Open,aggregate(list(ClosureNbrDSP=ClosureEventDSP,ClosureNbrASP=ClosureEventASP,ClosureNbrPSP=ClosureEventPSP,
                                                       ClosureDurationDSP=ClosureDSP,ClosureDurationASP=ClosureASP,ClosureDurationPSP=ClosurePSP,
                                                       ClosureNbrTot=ClosureEvent,ClosureDurationTot=Closure),list(FishingZone=Zone,Season=Saison),sum))


#### Export the synthesized tables of theoretical closures 
write.table(TheoClosures,file="TheoreticalClosures_Weeks.txt",sep="\t",quote=FALSE,row.names=FALSE)
write.table(TheoClosuresSeasons,file="TheoreticalClosures_Seasons.txt",sep="\t",quote=FALSE,row.names=FALSE)

