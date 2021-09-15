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
at.data.avg <- ddply(at.data, .(date), summarize, 
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
A<-sco2.data$at.mol#in mol/L
carb<-Carbfrom_C_A(K1, K2, StmpCO2.mol , A)
sco2.data$DIC.uM<-carb$D*1000*1000



# QAQC data and calculate slopes by treatment -----------------------------

# DIC ---------------------------------------------------------------------


setwd("C:/Users/matt/IDrive-Sync/Postdoc/Chamber metabolism/Chambers2021")
meta<-read.csv("chamber.notes.all.csv",header=T)
meta$start.date.time<-as.POSIXct(paste(meta$date, meta$start.time),format="%m-%d-%Y %H:%M:%S",tz="")
meta$end.date.time<-as.POSIXct(paste(meta$date, meta$end.time),format="%m-%d-%Y %H:%M:%S",tz="")
meta$end.date.time[7:36]<-meta$end.date.time[7:36]+(60*6)

## SUBSET
## Subset the main data file based on incubation time and experiment

subdat.DIC<-list()

for(i in 1:length(meta$level)){
  
  subdat.DIC[[i]]<-subset(sco2.data, date.time>meta$start.date.time[i]  & date.time< meta$end.date.time[i])
  subdat.DIC[[i]]$samplenum <-i
  
}
subdat.DIC <- bind_rows(subdat.DIC)

## Calculate the difference in time for each experiment
sco2.timediff<-NA
for (i in 1:length(meta$level)){
sco2.timediff[i]<-as.double(difftime(meta$end.date.time[i], meta$start.date.time[i], units="min"))
}

## Add time since each experiment started in minutes
subdat.DIC$time.since<-NA
for(i in 1:length(meta$level)){
time.sin<-as.numeric(seq(from=0, to=round(sco2.timediff[i]),by=5)) 
subdat.DIC[subdat.DIC$samplenum==i,]$time.since<-time.sin[1:length(subdat.DIC[subdat.DIC$samplenum==i,]$date)]
}


## QA/QC every sample

## Look at time series traces plotted over original data
## and decide on new Reset time if necessary


ggplot(subdat.DIC, aes(time.since, DIC.uM))+
  geom_point()+
  geom_smooth(method='lm', formula= y~x)+
  facet_wrap(~samplenum, ncol=4, scales = "free")

##Remove obviously bad data points
subdat.DIC$DIC.uM[c(45,50,59,60,76,85,92,93,98,111,133:135,156,160,163,250,259)]<-NA

#Replot
ggplot(subdat.DIC, aes(time.since, DIC.uM))+
  geom_point()+
  geom_smooth(method='lm', formula= y~x)+
  facet_wrap(~samplenum, ncol=4, scales = "free")

##Calculate the slope for each experiment
slope<-
  subdat.DIC %>% 
  group_by(samplenum) %>% 
  do({
    mod = lm(DIC.uM~ time.since , data = .)
    data.frame(DIC.uM.slope = coef(mod)[2])
  })

##Join slopes to final data object
meta.final<- full_join(meta,slope)



# DO ----------------------------------------------------------------------


## SUBSET

## Subset the main data file based on incubation time and experiment

subdat.DO<-list()

for(i in 1:length(meta$level)){
  
  subdat.DO[[i]]<-subset(o2.data, date.time>meta$start.date.time[i] & o2.data$date.time< meta$end.date.time[i])
  subdat.DO[[i]]$samplenum <-i
  
}
subdat.DO <- bind_rows(subdat.DO)

## Calculate the difference in time for each experiment
o2.timediff<-NA
for (i in 1:length(meta$level)){
  o2.timediff[i]<-as.double(difftime(meta$end.date.time[i], meta$start.date.time[i], units="min"))
}

## Add time since each experiment started in minutes
subdat.DO$time.since<-NA
for(i in 1:length(meta$level)){
  time.sin<-as.numeric(seq(from=0, to=round(o2.timediff[i]),by=0.5)) 
  subdat.DO[subdat.DO$samplenum==i,]$time.since<-time.sin[1:length(subdat.DO[subdat.DO$samplenum==i,]$date)]
}


## QA/QC every sample

## Look at time series traces plotted over original data
## and decide on new Reset time if necessary


ggplot(subdat.DO, aes(time.since, odo.uM))+
  geom_point()+
  geom_smooth(method='lm', formula= y~x)+
  facet_wrap(~samplenum, ncol=4,scales='free')

##Remove obviously bad data points
subdat.DO$odo.uM[c(1158:1198)]<-NA

##Replot
ggplot(subdat.DO, aes(time.since, odo.uM))+
  geom_point()+
  geom_smooth(method='lm', formula= y~x)+
  facet_wrap(~samplenum, ncol=4,scales='free')

##Calculate the slope for each experiment
slope<-
  subdat.DO %>% 
  group_by(samplenum) %>% 
  do({
    mod = lm(odo.uM~ time.since , data = .)
    data.frame(odo.uM.slope = coef(mod)[2])
  })

##Join to meta file
meta.final<- full_join(meta.final,slope)

#meta.final$DIC.uM.slope[3]<-0.15##average of other two dark delta DIC

#Substract ER from NEP to make GPP
meta.final$level<-as.factor(meta.final$level)
#for (i in 1:levels(meta.final$level))


# Calculate PQ ------------------------------------------------------------


PQ.DIC<-NA

PQ.DIC[1]<-meta.final$DIC.uM.slope[1]+meta.final$DIC.uM.slope[2]
PQ.DIC[2]<-meta.final$DIC.uM.slope[3]+meta.final$DIC.uM.slope[4]
PQ.DIC[3]<-meta.final$DIC.uM.slope[5]+meta.final$DIC.uM.slope[6]
PQ.DIC[4]<-meta.final$DIC.uM.slope[7]+meta.final$DIC.uM.slope[8]
PQ.DIC[5]<-meta.final$DIC.uM.slope[9]+meta.final$DIC.uM.slope[10]#bad
PQ.DIC[6]<-meta.final$DIC.uM.slope[11]+meta.final$DIC.uM.slope[12]#bad
PQ.DIC[7]<-meta.final$DIC.uM.slope[13]+meta.final$DIC.uM.slope[14]#bad
PQ.DIC[8]<-meta.final$DIC.uM.slope[15]+meta.final$DIC.uM.slope[16]#bad
PQ.DIC[9]<-meta.final$DIC.uM.slope[17]+meta.final$DIC.uM.slope[18]
PQ.DIC[10]<-meta.final$DIC.uM.slope[19]+meta.final$DIC.uM.slope[20]
PQ.DIC[11]<-meta.final$DIC.uM.slope[21]+meta.final$DIC.uM.slope[22]
PQ.DIC[12]<-meta.final$DIC.uM.slope[23]+meta.final$DIC.uM.slope[24]
PQ.DIC[13]<-meta.final$DIC.uM.slope[25]+meta.final$DIC.uM.slope[26]
PQ.DIC[14]<-meta.final$DIC.uM.slope[27]+meta.final$DIC.uM.slope[28]#bad
PQ.DIC[15]<-meta.final$DIC.uM.slope[29]+meta.final$DIC.uM.slope[30]#bad
PQ.DIC[16]<-meta.final$DIC.uM.slope[31]+meta.final$DIC.uM.slope[32]
PQ.DIC[17]<-meta.final$DIC.uM.slope[33]+meta.final$DIC.uM.slope[34]#bad
PQ.DIC[18]<-meta.final$DIC.uM.slope[35]+meta.final$DIC.uM.slope[36]#bad

PQ.DO<-NA

PQ.DO[1]<-meta.final$odo.uM.slope[1]+meta.final$odo.uM.slope[2]
PQ.DO[2]<-meta.final$odo.uM.slope[3]+meta.final$odo.uM.slope[4]
PQ.DO[3]<-meta.final$odo.uM.slope[5]+meta.final$odo.uM.slope[6]
PQ.DO[4]<-meta.final$odo.uM.slope[7]+meta.final$odo.uM.slope[8]
PQ.DO[5]<-meta.final$odo.uM.slope[9]+meta.final$odo.uM.slope[10]#bad DIC
PQ.DO[6]<-meta.final$odo.uM.slope[11]+meta.final$odo.uM.slope[12]#bad DIC
PQ.DO[7]<-meta.final$odo.uM.slope[13]+meta.final$odo.uM.slope[14]#bad DIC
PQ.DO[8]<-meta.final$odo.uM.slope[15]+meta.final$odo.uM.slope[16]#bad DIC
PQ.DO[9]<-meta.final$odo.uM.slope[17]+meta.final$odo.uM.slope[18]
PQ.DO[10]<-meta.final$odo.uM.slope[19]+meta.final$odo.uM.slope[20]
PQ.DO[11]<-meta.final$odo.uM.slope[21]+meta.final$odo.uM.slope[22]
PQ.DO[12]<-meta.final$odo.uM.slope[23]+meta.final$odo.uM.slope[24]
PQ.DO[13]<-meta.final$odo.uM.slope[25]+meta.final$odo.uM.slope[26]#bad O2
PQ.DO[14]<-meta.final$odo.uM.slope[27]+meta.final$odo.uM.slope[28]#bad DIC
PQ.DO[15]<-meta.final$odo.uM.slope[29]+meta.final$odo.uM.slope[30]#bad DIC
PQ.DO[16]<-meta.final$odo.uM.slope[31]+meta.final$odo.uM.slope[32]
PQ.DO[17]<-meta.final$odo.uM.slope[33]+meta.final$odo.uM.slope[34]#bad DIC
PQ.DO[18]<-meta.final$odo.uM.slope[35]+meta.final$odo.uM.slope[36]#bad DIC
PQ<-PQ.DIC/PQ.DO

PQ[c(5,6,7,8,13,14,15,17,18)]<-NA
abs(PQ)
samplenum<-seq(from=2, to=36, by=2)
abs.PQ<-as.data.frame(cbind(abs(PQ),samplenum))

meta.final<- full_join(meta.final,abs.PQ)





# Calculate ancillary variables ------------------------------------------

#light
meta<-read.csv("chamber.notes.all.csv",header=T)
meta$start.date.time<-as.POSIXct(paste(meta$date, meta$start.time),format="%m-%d-%Y %H:%M:%S",tz="")
meta$end.date.time<-as.POSIXct(paste(meta$date, meta$end.time),format="%m-%d-%Y %H:%M:%S",tz="")

light<-read.csv("chamber.light.all.csv",header=T)
names(light)<-c("date.time", "temp", "brightness")
light$date.time<-as.POSIXct(light$date.time,format="%m-%d-%Y %H:%M:%S",tz="")

subdat<-list()

for(i in 7:length(meta$level)){
  
  subdat[[i]]<-subset(light, date.time>meta$start.date.time[i]  & date.time< meta$end.date.time[i])
  subdat[[i]]$samplenum <-i
  
}
subdat <- bind_rows(subdat)

## Look at time series traces plotted over original data
## and decide on new Reset time if necessary


ggplot(subdat, aes(date.time, brightness))+
  geom_point()+
  facet_wrap(~samplenum, ncol=4,scales='free')


##Calculate the mean for each experiment
result<-ddply(subdat,.(samplenum), summarize,
              mean.light=mean(brightness,na.rm=TRUE))

##Join to meta file
meta.final<- full_join(meta.final,result)

##Plot
ggplot(meta.final, aes(DIC.uM.slope, light, color=treatment))+
  geom_point()
  facet_grid(~metab)
  
#temperature
  meta<-read.csv("chamber.notes.all.csv",header=T)
  meta$start.date.time<-as.POSIXct(paste(meta$date, meta$start.time),format="%m-%d-%Y %H:%M:%S",tz="")
  meta$end.date.time<-as.POSIXct(paste(meta$date, meta$end.time),format="%m-%d-%Y %H:%M:%S",tz="")
  
  
  subdat<-list()
  
  for(i in 1:length(meta$level)){
    
    subdat[[i]]<-subset(o2.data, date.time>meta$start.date.time[i]  & date.time< meta$end.date.time[i])
    subdat[[i]]$samplenum <-i
    
  }
  subdat <- bind_rows(subdat)
  
  ## Look at time series traces plotted over original data
  ## and decide on new Reset time if necessary
  
  
  ggplot(subdat, aes(date.time, temp.c))+
    geom_point()+
    facet_wrap(~samplenum, ncol=4,scales='free')
  
  
  ##Calculate the mean for each experiment
  result.temp<-ddply(subdat,.(samplenum), summarize,
                mean.temp=mean(temp.c,na.rm=TRUE))
  
  ##Join to meta file
  meta.final<- full_join(meta.final,result.temp)
  
  ##Plot
  ggplot(meta.final, aes(mean.temp, DIC.uM.slope, color=treatment, shape=level))+
    geom_point()+
  facet_grid(~metab)
  
  
##Summary plot
  
good.sampnum<-meta.final$samplenum[which(!is.na(PQ))]
subdat.DIC[subdat.DIC$samplnum %in% good.sampnum,]
good.DIC<-subset(subdat.DIC,samplenum==c(1,2,3))
good.DO<-subset()
  
  ##Plot
  ggplot(meta.final, aes(treatment, V1, color=treatment, shape=level))+
    geom_point()+
    facet_grid(metab~samplenum)


