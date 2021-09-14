library(lubridate)
library(data.table)
library(plyr)
library(dplyr)
library(readr)     
library(ggplot2)
library(cowplot)
library(reshape2)
library(xts)
library(dygraphs)
library(tidyr)

# Read/format chamber SAMI CO2 data ---------------------------------------------

setwd("C:/Users/matt/IDrive-Sync/Postdoc/Chamber metabolism/Chambers2021")
sco2.data = read.csv("chamber.co2.all.csv",header=T)
sco2.data$date.time<-as.POSIXct(paste(sco2.data$date, sco2.data$time),format="%m-%d-%Y %H:%M:%S",tz="UTC")
sco2.data$date.time<-format(sco2.data$date.time, tz="America/Denver")
sco2.data$date.time<-as.POSIXct(sco2.data$date.time,format="%Y-%m-%d %H:%M:%S")
sco2.data$date.time.round<-round_date(sco2.data$date.time, unit = "minute") #for matching with O2
names(sco2.data)<-c("date", "time","temp_sami","CO2.ppmv","battery","date.time", "date.time.round")


# Correct for temperature difference --------------------------------------

#SAMI is measuring temp of the bath. Must adjust to temp in the chamber (from DO)
o2.data = read.csv("chamber.o2.all.csv",header=T)
names(o2.data)<-c("date", "time","pres.mmhg","odo.sat","odo.mgl","odo.local","temp.f")
o2.data$date.time<-as.POSIXct(paste(o2.data$date, o2.data$time),format="%m-%d-%Y %H:%M:%S",tz="")
o2.data$temp.c<-(o2.data$temp.f-32)*5/9
o2.data$odo.uM<-o2.data$odo.mgl*1000/31.99

#New object with average temp per minute (temp with DO is measured every 30 s)
o2.avg.temp<-NA
o2.avg.temp$date.time<-o2.data$date.time
o2.avg.temp$temp.c<-o2.data$temp.c
o2.avg.temp<-as.data.table(o2.avg.temp)
avg.temp<-o2.avg.temp[, c('Month', 'Day','Year','Hour', 'Minute') := .(data.table::month(o2.avg.temp$date.time),mday(o2.avg.temp$date.time),year(o2.avg.temp$date.time),hour(o2.avg.temp$date.time), minute(o2.avg.temp$date.time))
][, .(Avg = mean(temp.c)), .(Month,Day,Year,Hour, Minute)]

#recreate the date
avg.temp$date.time<-as.POSIXct(ISOdate(month=avg.temp$Month, day=avg.temp$Day, year=avg.temp$Year, hour=avg.temp$Hour,min=avg.temp$Minute, sec=0,tz=""))

#match dates between O2 and CO2 objects and then add the temp from O2 to the CO2 object
sco2.data$o2.temp<-avg.temp$Avg[match(sco2.data$date.time.round,avg.temp$date.time)] 

#finally, correct the co2 data with the correct temperature (equation from Mike and Cory)
sco2.data$CO2_corrected = (sco2.data$CO2) * exp(0.0423*(sco2.data$o2.temp - sco2.data$temp_sami))


# Transform from ppmv to mol ----------------------------------------------

## SOLUBILITY CONSTANT FOR CO2 (KH) 
## Use Henry's Law and parameters from Weiss 1974 to calculate KH in units of mol L-1 atm-1 as a function of temperature (also in Demarty et al 2011; many others)
KH.CO2 <- function(temp_sami){
  tempK <- temp_sami + 273.15
  KH.CO2 <- exp( -58.0931 + (90.5069*(100/tempK)) +22.2940*log((tempK/100), base=exp(1)) )
  KH.CO2
}

## ppmv to mol/L
co2 <- function(temp_sami,CO2.ppmv){
  KH.samp <- KH.CO2(temp_sami)
  co2.mol.L<-  KH.samp*CO2.ppmv/1000000
}

#Run the function
sco2.data$CO2.mol<-co2(sco2.data$o2.temp, sco2.data$CO2.ppmv)


# Calculate DIC -----------------------------------------------------------

#Load and match alkalinity data
at.data <- read.csv("chamber.alk.all.csv",header=T)
at.data.avg <- ddply(at.data, .(Date.collected), summarize, 
                     mean(at.mol.L, na.rm=TRUE))
names(at.data.avg)<-c("date", "at.mol")
sco2.data$at.mol<-at.data.avg$at.mol[match(sco2.data$date,at.data.avg$date)]


## Equilibrium constants
K1calc<- function(temp.c) { 10^( (-3404.71/(273.15+temp.c)) + 14.844 -0.033*(temp.c+273.15) )}
K2calc<- function(temp.c) { 10^( (-2902.39/(273.15+temp.c)) + 6.498 -0.0238*(temp.c+273.15) )}



# Carbonate Chemistry from C & A
#A=alkalinity,  StmpCO2.mol=CO2 (mol_L), B=bicarbonate, H_minus=Hydroxide, Ca=carbonate, D=DIC, A_check=alkalinity check
Carbfrom_C_A <- function(K1, K2, StmpCO2.mol , A){
  H <- (((-K1*StmpCO2.mol ))-sqrt(((K1*StmpCO2.mol )^2)-(4*-1*A*2*K1*K2*StmpCO2.mol )))/(2*-1*A)
  pH <- -1*log10((H))
  B <- (K1*StmpCO2.mol )/H
  Ca <- (K2*B)/H
  D <- StmpCO2.mol  + B + Ca
  A_check <- B + 2*Ca
  Carb1 <- list(H, pH, StmpCO2.mol , B, Ca, D, A, A_check)
  names(Carb1) <- c("H", "pH", "StmpCO2.mol ", "B", "Ca", "D", "A","A check")
  Carb1
}




# Load variables and calculate DIC
K1<-K1calc(sco2.data$o2.temp)
K2<-K2calc(sco2.data$o2.temp)
StmpCO2.mol<-sco2.data$CO2.mol #in mol/L
A<-sco2.data$at.mol*2#in mol/L
carb<-Carbfrom_C_A(K1, K2, StmpCO2.mol , A)
sco2.data$DIC.uM<-carb$D*1000*1000



# QAQC data and calculate slopes by treatment -----------------------------

# DIC ---------------------------------------------------------------------


setwd("C:/Users/matt/IDrive-Sync/Postdoc/Chamber metabolism/Chambers2021")
meta<-read.csv("chamber.notes.all.csv",header=T)
meta$start.date.time<-as.POSIXct(paste(meta$date, meta$start.time),format="%m-%d-%Y %H:%M:%S",tz="")
meta$end.date.time<-as.POSIXct(paste(meta$date, meta$end.time),format="%m-%d-%Y %H:%M:%S",tz="")
##############
## SUBSET
#############
## Subset the main data file based on incubation time and experiment

subdat<-list()

for(i in 1:length(meta$level)){
  
  subdat[[i]]<-subset(sco2.data, date.time>meta$start.date.time[i]  & date.time< meta$end.date.time[i])
  subdat[[i]]$samplenum <-i
  
}
subdat <- bind_rows(subdat)

sco2.timediff<-NA
for (i in 1:length(meta$level)){
sco2.timediff[i]<-as.double(difftime(meta$end.date.time[i], meta$start.date.time[i], units="min"))
}

subdat$time.since<-NA
for(i in 1:length(meta$level)){
time.sin<-as.numeric(seq(from=0, to=round(sco2.timediff[i]),by=5)) 
subdat[subdat$samplenum==i,]$time.since<-time.sin[1:length(subdat[subdat$samplenum==i,]$date)]
}

######################
## QA/QC every sample
######################
## Look at time series traces plotted over original data
## and decide on new Reset time if necessary


ggplot(subdat, aes(time.since, DIC.uM))+
  geom_point()+
  geom_smooth(method='lm', formula= y~x)+
  facet_wrap(~samplenum, ncol=2, scales = "free")

##Calculate the slope
slope<-
  subdat %>% 
  group_by(samplenum) %>% 
  do({
    mod = lm(DIC.uM~ time.since , data = .)
    data.frame(DIC.uM.slope = coef(mod)[2])
  })


meta.final<- full_join(meta,slope)



# DO ----------------------------------------------------------------------

##############
## SUBSET
#############
## Subset the main data file based on incubation time and experiment

subdat<-list()

for(i in 1:length(meta$level)){
  
  subdat[[i]]<-subset(o2.data, date.time>meta$start.date.time[i] & o2.data$date.time< meta$end.date.time[i])
  subdat[[i]]$samplenum <-i
  
}
subdat <- bind_rows(subdat)

o2.timediff<-NA
for (i in 1:length(meta$level)){
  o2.timediff[i]<-as.double(difftime(meta$end.date.time[i], meta$start.date.time[i], units="min"))
}

subdat$time.since<-NA
for(i in 1:length(meta$level)){
  time.sin<-as.numeric(seq(from=0, to=round(o2.timediff[i]),by=0.5)) 
  subdat[subdat$samplenum==i,]$time.since<-time.sin[1:length(subdat[subdat$samplenum==i,]$date)]
}

######################
## QA/QC every sample
######################
## Look at time series traces plotted over original data
## and decide on new Reset time if necessary


ggplot(subdat, aes(time.since, odo.uM))+
  geom_point()+
  geom_smooth(method='lm', formula= y~x)+
  facet_wrap(~samplenum, ncol=2)

##Calculate the slope
slope<-
  subdat %>% 
  group_by(samplenum) %>% 
  do({
    mod = lm(odo.uM~ time.since , data = .)
    data.frame(odo.uM.slope = coef(mod)[2])
  })


meta.final<- full_join(meta.final,slope)

meta.final$DIC.uM.slope/meta.final$odo.uM.slope

#Substract ER from NEP to make GPP
meta.final$level<-as.factor(meta.final$level)
for (i in 1:levels(meta.final$level))

PQ.DIC<-NA
PQ.DIC[1]<-meta.final$DIC.uM.slope[1]+meta.final$DIC.uM.slope[2]
PQ.DIC[2]<-meta.final$DIC.uM.slope[3]+meta.final$DIC.uM.slope[4]
PQ.DIC[3]<-meta.final$DIC.uM.slope[5]+meta.final$DIC.uM.slope[6]

PQ.DO<-NA
PQ.DO[1]<-meta.final$odo.uM.slope[1]+meta.final$odo.uM.slope[2]
PQ.DO[2]<-meta.final$odo.uM.slope[3]+meta.final$odo.uM.slope[4]
PQ.DO[3]<-meta.final$odo.uM.slope[5]+meta.final$odo.uM.slope[6]

PQ<-PQ.DIC/PQ.DO

abs(PQ)

write.csv(sco2.data,"Chamber_test_run_SAMI_DIC_08.18.2020.csv",row.names = FALSE)
