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
meta$end.date.time<-meta$end.date.time+(60*6)
meta$start.date.time<-meta$start.date.time-(60*6)

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

sub<-subset(subdat.DIC, samplenum==1|samplenum==2 |samplenum==3|samplenum==4|samplenum==5|samplenum==6)
ggplot(subdat.DIC, aes(time.since, CO2.ppmv))+
  geom_point()+
  geom_smooth(method='lm', formula= y~x)+
  facet_wrap(~samplenum, ncol=4, scales = "free")

subdat.DIC %>% 
  filter(treatment=='vel', level==1) %>% 
  ggplot(aes(time.since, DIC.uM)) %.>%
  geom_point() %.>%
  geom_smooth(method='lm', formula= y~x) %.>%
  facet_wrap(metab~date,  scales = "free")





##Remove obviously bad data points
subdat.DIC$DIC.uM[c(45,50,59,60,76,85,92,93,98,111,133:135,156,160,163,250,259)]<-NA

#Replot
ggplot(subdat.DIC, aes(time.since, DIC.uM))+
  geom_point()+
  geom_smooth(method='lm', formula= y~x)+
  facet_wrap(~samplenum, ncol=4, scales = "free")

##Add metadata
subdat.DIC<- full_join(meta,subdat.DIC)


##Calculate the slope for each experiment
slope.DIC<-
  subdat.DIC %>% 
  group_by(samplenum) %>% 
  do({
    mod = lm(DIC.uM~ time.since , data = .)
    data.frame(DIC.uM.slope = coef(mod)[2])
  })

##Calculate the slope for each experiment
slope.CO2<-
  subdat.DIC %>% 
  group_by(samplenum) %>% 
  do({
    mod = lm(CO2.ppmv~ time.since , data = .)
    data.frame(CO2.ppmv.slope = coef(mod)[2])
  })


##Join slopes to final data object
meta.final<- full_join(meta,slope.DIC)
meta.final<- full_join(meta.final,slope.CO2)


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

##Add metadata
subdat.DO<- full_join(meta,subdat.DO)

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




# Calculate PQ ------------------------------------------------------------


PQ.DIC<-NA
#Add ER with NEP to make GPP
PQ.DIC[1]<-meta.final$DIC.uM.slope[2]-meta.final$DIC.uM.slope[1]
PQ.DIC[2]<-meta.final$DIC.uM.slope[4]-meta.final$DIC.uM.slope[3]
PQ.DIC[3]<-meta.final$DIC.uM.slope[6]-meta.final$DIC.uM.slope[5]
PQ.DIC[4]<-meta.final$DIC.uM.slope[8]-meta.final$DIC.uM.slope[7]
PQ.DIC[5]<-meta.final$DIC.uM.slope[10]-meta.final$DIC.uM.slope[9]#bad
PQ.DIC[6]<-meta.final$DIC.uM.slope[12]-meta.final$DIC.uM.slope[11]#bad
PQ.DIC[7]<-meta.final$DIC.uM.slope[14]-meta.final$DIC.uM.slope[13]#bad
PQ.DIC[8]<-meta.final$DIC.uM.slope[16]-meta.final$DIC.uM.slope[15]#bad
PQ.DIC[9]<-meta.final$DIC.uM.slope[18]-meta.final$DIC.uM.slope[17]
PQ.DIC[10]<-meta.final$DIC.uM.slope[20]-meta.final$DIC.uM.slope[19]
PQ.DIC[11]<-meta.final$DIC.uM.slope[22]-meta.final$DIC.uM.slope[21]
PQ.DIC[12]<-meta.final$DIC.uM.slope[24]-meta.final$DIC.uM.slope[23]
PQ.DIC[13]<-meta.final$DIC.uM.slope[26]-meta.final$DIC.uM.slope[25]
PQ.DIC[14]<-meta.final$DIC.uM.slope[28]-meta.final$DIC.uM.slope[27]#bad
PQ.DIC[15]<-meta.final$DIC.uM.slope[30]-meta.final$DIC.uM.slope[29]#bad
PQ.DIC[16]<-meta.final$DIC.uM.slope[32]-meta.final$DIC.uM.slope[31]
PQ.DIC[17]<-meta.final$DIC.uM.slope[34]-meta.final$DIC.uM.slope[33]#bad
PQ.DIC[18]<-meta.final$DIC.uM.slope[36]-meta.final$DIC.uM.slope[35]#bad


PQ.CO2<-NA
#Add ER with NEP to make GPP
PQ.CO2[1]<-meta.final$CO2.ppmv.slope[2]-meta.final$CO2.ppmv.slope[1]
PQ.CO2[2]<-meta.final$CO2.ppmv.slope[4]-meta.final$CO2.ppmv.slope[3]
PQ.CO2[3]<-meta.final$CO2.ppmv.slope[6]-meta.final$CO2.ppmv.slope[5]
PQ.CO2[4]<-meta.final$CO2.ppmv.slope[8]-meta.final$CO2.ppmv.slope[7]
PQ.CO2[5]<-meta.final$CO2.ppmv.slope[10]-meta.final$CO2.ppmv.slope[9]#bad
PQ.CO2[6]<-meta.final$CO2.ppmv.slope[12]-meta.final$CO2.ppmv.slope[11]#bad
PQ.CO2[7]<-meta.final$CO2.ppmv.slope[14]-meta.final$CO2.ppmv.slope[13]#bad
PQ.CO2[8]<-meta.final$CO2.ppmv.slope[16]-meta.final$CO2.ppmv.slope[15]#bad
PQ.CO2[9]<-meta.final$CO2.ppmv.slope[18]-meta.final$CO2.ppmv.slope[17]
PQ.CO2[10]<-meta.final$CO2.ppmv.slope[20]-meta.final$CO2.ppmv.slope[19]
PQ.CO2[11]<-meta.final$CO2.ppmv.slope[22]-meta.final$CO2.ppmv.slope[21]
PQ.CO2[12]<-meta.final$CO2.ppmv.slope[24]-meta.final$CO2.ppmv.slope[23]
PQ.CO2[13]<-meta.final$CO2.ppmv.slope[26]-meta.final$CO2.ppmv.slope[25]
PQ.CO2[14]<-meta.final$CO2.ppmv.slope[28]-meta.final$CO2.ppmv.slope[27]#bad
PQ.CO2[15]<-meta.final$CO2.ppmv.slope[30]-meta.final$CO2.ppmv.slope[29]#bad
PQ.CO2[16]<-meta.final$CO2.ppmv.slope[32]-meta.final$CO2.ppmv.slope[31]
PQ.CO2[17]<-meta.final$CO2.ppmv.slope[34]-meta.final$CO2.ppmv.slope[33]#bad
PQ.CO2[18]<-meta.final$CO2.ppmv.slope[36]-meta.final$CO2.ppmv.slope[35]#bad


PQ.DO<-NA

PQ.DO[1]<-meta.final$odo.uM.slope[2]-meta.final$odo.uM.slope[1]
PQ.DO[2]<-meta.final$odo.uM.slope[4]-meta.final$odo.uM.slope[3]
PQ.DO[3]<-meta.final$odo.uM.slope[6]-meta.final$odo.uM.slope[5]
PQ.DO[4]<-meta.final$odo.uM.slope[8]-meta.final$odo.uM.slope[7]
PQ.DO[5]<-meta.final$odo.uM.slope[10]-meta.final$odo.uM.slope[9]#bad DIC
PQ.DO[6]<-meta.final$odo.uM.slope[12]-meta.final$odo.uM.slope[11]#bad DIC
PQ.DO[7]<-meta.final$odo.uM.slope[14]-meta.final$odo.uM.slope[13]#bad DIC
PQ.DO[8]<-meta.final$odo.uM.slope[16]-meta.final$odo.uM.slope[15]#bad DIC
PQ.DO[9]<-meta.final$odo.uM.slope[18]-meta.final$odo.uM.slope[17]
PQ.DO[10]<-meta.final$odo.uM.slope[20]-meta.final$odo.uM.slope[19]
PQ.DO[11]<-meta.final$odo.uM.slope[22]-meta.final$odo.uM.slope[21]
PQ.DO[12]<-meta.final$odo.uM.slope[24]-meta.final$odo.uM.slope[23]
PQ.DO[13]<-meta.final$odo.uM.slope[26]-meta.final$odo.uM.slope[25]#bad O2
PQ.DO[14]<-meta.final$odo.uM.slope[28]-meta.final$odo.uM.slope[27]#bad DIC
PQ.DO[15]<-meta.final$odo.uM.slope[30]-meta.final$odo.uM.slope[29]#bad DIC
PQ.DO[16]<-meta.final$odo.uM.slope[32]-meta.final$odo.uM.slope[31]
PQ.DO[17]<-meta.final$odo.uM.slope[34]-meta.final$odo.uM.slope[33]#bad DIC
PQ.DO[18]<-meta.final$odo.uM.slope[36]-meta.final$odo.uM.slope[35]#bad DIC
PQ<-PQ.DO/PQ.DIC
PQ.C<-PQ.DO/PQ.CO2
PQ[c(4,5,6,7,8,9,13,14,15,16,17,18)]<-NA
PQ.C[c(4,5,6,7,8,9,13,14,15,16,17,18)]<-NA
abs(PQ)
abs(PQ.C)
samplenum<-seq(from=2, to=36, by=2)
abs.PQ<-as.data.frame(cbind(abs(PQ),samplenum,PQ.DO, PQ.DIC))
names(abs.PQ)<-c('PQ','samplenum', 'PQ.DO','PQ.DIC')
meta.final<- full_join(meta.final,abs.PQ)


# Calculate rates ---------------------------------------------------------

volume<-17 #Liters
area<-0.03 #m^2
o2.mass<-32
c.mass<-12
rate<-NA
meta.final$GPP.do<-meta.final$PQ.DO*volume/area*60/1000
rate$ER.do<-subset(meta.final, metab=="dark")$odo.uM.slope*volume/area*60/1000/1000

meta.final$GPP.dic<-meta.final$PQ.DIC*volume/area*60/1000
rate$ER.dic<-subset(meta.final, metab=="dark")$DIC.uM.slope*volume/area*60/1000/1000
rate$GPP.do[c(4,5,6,7,8,9,13,14,15,16,17,18)]<-NA
rate$GPP.dic[c(4,5,6,7,8,9,13,14,15,16,17,18)]<-NA
rate$ER.do[c(4,5,6,7,8,9,13,14,15,16,17,18)]<-NA
rate$ER.dic[c(4,5,6,7,8,9,13,14,15,16,17,18)]<-NA

meta.final$GPP.dic[c(8,10,12,14,16,18,26,28,30,32,34,36)]<-NA
meta.final$GPP.do[c(8,10,12,14,16,18,26,28,30,32,34,36)]<-NA


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
meta.final$metab<-as.factor(meta.final$metab)
meta.final$treatment<-as.factor(meta.final$treatment)
meta.final$level<-as.factor(meta.final$level)


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
  
  
meta.final.good<-subset(meta.final, !is.na(meta.final$PQ))
  
  
##Summary plot
  # New facet label names for dose variable
  metab.labs <- c("Dark", "Light")
  names(metab.labs ) <- c("dark", "light")
  
  # New facet label names for supp variable
  level.prod.labs <- c("High Prod.", "Med Prod.", "Low Prod.")
  names(level.prod.labs) <- c("1", "2", "3")  
  
  # New facet label names for supp variable
  level.vel.labs <- c("Low Vel.","Med Vel." , "High Vel.")
  names(level.vel.labs) <- c("1", "2", "3")  
  
  
#velocity
#DIC
good.sampnum<-meta.final$samplenum[which(!is.na(PQ))]
good.DIC.vel<-subset(subdat.DIC,samplenum==1|samplenum==2|samplenum==3|samplenum==4|samplenum==5|samplenum==6)
good.DIC.vel$level<-as.factor(good.DIC.vel$level)
good.DIC.vel$level<-factor(good.DIC.vel$level,levels=(c("1","2","3")))
good.DIC.vel$metab<-as.factor(good.DIC.vel$metab)

##Plot

ggplot(good.DIC.vel, aes(time.since, DIC.uM))+
  geom_point()+
  geom_smooth(method='lm', formula= y~x)+
  facet_grid(level~metab, scales='free',
             labeller = labeller(level=level.vel.labs, metab=metab.labs))+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=16,color="black"))+
  theme(axis.title.y=element_text(size=16,color="black"))+
  theme(axis.text.y=element_text(size=16,color="black"))+
  theme(axis.text.x=element_text(size=16,color="black"))+
  xlab("Time (mins)")+
  ylab("DIC (uM)")

#DO
good.sampnum<-meta.final$samplenum[which(!is.na(PQ))]
good.DO.vel<-subset(subdat.DO,samplenum==1|samplenum==2|samplenum==3|samplenum==4|samplenum==5|samplenum==6)
good.DO.vel$level<-as.factor(good.DO.vel$level)
good.DO.vel$level<-factor(good.DO.vel$level,levels=(c("1","2","3")))
good.DO.vel$metab<-as.factor(good.DO.vel$metab)

##Plot

ggplot(good.DO.vel, aes(time.since, odo.uM))+
  geom_point()+
  geom_smooth(method='lm', formula= y~x)+
  facet_grid(level~metab, scales='free',
             labeller = labeller(level=level.vel.labs, metab=metab.labs))+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=16,color="black"))+
  theme(axis.title.y=element_text(size=16,color="black"))+
  theme(axis.text.y=element_text(size=16,color="black"))+
  theme(axis.text.x=element_text(size=16,color="black"))+
  xlab("Time (mins)")+
  ylab("DO(uM)")




#productivity
#DIC
  good.DIC.prod<-subset(subdat.DIC,samplenum==19|samplenum==20|samplenum==21|samplenum==22|samplenum==23|samplenum==24)
  good.DIC.prod$level<-as.factor(good.DIC.prod$level)
  good.DIC.prod$level<-factor(good.DIC.prod$level,levels=(c("3","2","1")))
  good.DIC.prod$metab<-as.factor(good.DIC.prod$metab)

  ##Plot
 
  
  ggplot(good.DIC.prod, aes(time.since, DIC.uM))+
    geom_point()+
    geom_smooth(method='lm', formula= y~x)+
    facet_grid(level~metab, scales='free',
               labeller = labeller(level=level.prod.labs, metab=metab.labs))+
    theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_text(size=16,color="black"))+
    theme(axis.title.y=element_text(size=16,color="black"))+
    theme(axis.text.y=element_text(size=16,color="black"))+
    theme(axis.text.x=element_text(size=16,color="black"))+
    xlab("Time (mins)")+
    ylab("DIC (uM)")

#DO
  good.DO.prod<-subset(subdat.DO,samplenum==19|samplenum==20|samplenum==21|samplenum==22|samplenum==23|samplenum==24)
  good.DO.prod$level<-as.factor(good.DO.prod$level)
  good.DO.prod$level<-factor(good.DO.prod$level,levels=(c("3","2","1")))
  good.DO.prod$metab<-as.factor(good.DO.prod$metab)
  ##Plot
  
  
  ggplot(good.DO.prod, aes(time.since, odo.uM))+
    geom_point()+
    geom_smooth(method='lm', formula= y~x)+
    facet_grid(level~metab, scales='free',
               labeller = labeller(level=level.prod.labs, metab=metab.labs))+
    theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_text(size=16,color="black"))+
    theme(axis.title.y=element_text(size=16,color="black"))+
    theme(axis.text.y=element_text(size=16,color="black"))+
    theme(axis.text.x=element_text(size=16,color="black"))+
    xlab("Time (mins)")+
    ylab("DO (uM)")
  
##Plots of raw data
  
good.DO<-subset(subdat.DO, samplenum==5 |samplenum==6)

ggplot(good.DO, aes(time.since, odo.uM))+
  geom_point()+
  geom_smooth(method='lm', formula= y~x)+
  facet_grid(level~metab, scales='free',
             labeller = labeller( metab=metab.labs))+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=16,color="black"))+
  theme(axis.title.y=element_text(size=16,color="black"))+
  theme(axis.text.y=element_text(size=16,color="black"))+
  theme(axis.text.x=element_text(size=16,color="black"))+
  xlab("Time (min.)")+
  ylab(expression(paste("DO (",mu,"M ",O[2],")")))+
  theme(strip.text.y = element_blank(),strip.text.x = element_text(size = 12))

good.DIC<-subset(subdat.DIC, samplenum==5 |samplenum==6)

ggplot(good.DIC, aes(time.since, DIC.uM))+
  geom_point()+
  geom_smooth(method='lm', formula= y~x)+
  facet_grid(level~metab, scales='free',
             labeller = labeller( metab=metab.labs))+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=16,color="black"))+
  theme(axis.title.y=element_text(size=16,color="black"))+
  theme(axis.text.y=element_text(size=16,color="black"))+
  theme(axis.text.x=element_text(size=16,color="black"))+
  xlab("Time (min.)")+
  ylab(expression(paste("DIC (",mu,"M C)")))+
  theme(strip.text.y = element_blank(),strip.text.x = element_text(size = 12))


## productivity vs PQ
ggplot(meta.final, aes(GPP.dic,PQ))+
  geom_point(size=3)+
  geom_smooth(method='lm', formula= y~x,se = FALSE)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=18,color="black"))+
  theme(axis.title.y=element_text(size=18,color="black"))+
  theme(axis.text.y=element_text(size=18,color="black"))+
  theme(axis.text.x=element_text(size=18,color="black"))+
  xlab("GPP")+
  xlim(c(-30,1))+
  ylim(c(0.5,3.5))+
  xlab(expression(paste("GPP (",mu,"moles DIC ",m^-2," ",hr^-1 ,")")))

ggplot(meta.final, aes(GPP.do,PQ))+
  geom_point(size=3)+
  geom_smooth(method='lm', formula= y~x,se = FALSE)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=18,color="black"))+
  theme(axis.title.y=element_text(size=18,color="black"))+
  theme(axis.text.y=element_text(size=18,color="black"))+
  theme(axis.text.x=element_text(size=18,color="black"))+
  xlab("GPP")+
  ylab("PQ")+
  xlim(c(0,30))+
  ylim(c(0.5,3.5))+
  xlab(expression(paste("GPP (",mu,"moles ", O[2]," ", m^-2," ",hr^-1 ,")")))

ggplot(meta.final, aes(mean.light,PQ))+
  geom_point(size=3)+
  geom_smooth(method='lm', formula= y~x,se = FALSE)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=16,color="black"))+
  theme(axis.title.y=element_text(size=16,color="black"))+
  theme(axis.text.y=element_text(size=16,color="black"))+
  theme(axis.text.x=element_text(size=16,color="black"))+
  xlab("GPP")+
  ylab("PQ")
  #xlim(c(0,60000))
  #ylim(c(0.5,3.5))+
  #xlab(expression(paste("GPP (",mu,"moles ", O[2]," ", m^-2," ",hr^-1 ,")")))

lm(PQ~GPP.do, data=meta.final)

ggplot(meta.final, aes(x=abs(GPP.dic), y=GPP.do), size=PQ)+
  geom_point(aes(size=PQ))+
  geom_abline(intercept = 0, slope = 1, linetype=2, size=1.2)+
  #geom_smooth(method='lm', formula= y~x,se = FALSE)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=18,color="black"))+
  theme(axis.title.y=element_text(size=18,color="black"))+
  theme(axis.text.y=element_text(size=18,color="black"))+
  theme(axis.text.x=element_text(size=18,color="black"))+
  xlim(c(0,30))+
  ylim(c(0,30))+
  xlab(expression(paste("GPP (",mu,"mol ", O[2]," ", m^-2," ",hr^-1 ,")")))+
  ylab(expression(paste("GPP (",mu,"mol DIC " , m^-2," ",hr^-1 ,")")))
  
ggplot(meta.final, aes(y=PQ))+
  geom_boxplot(width=0.5)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=18,color="black"))+
  theme(axis.title.y=element_text(size=18,color="black"))+
  theme(axis.text.y=element_text(size=18,color="black"))+
  theme(axis.text.x=element_text(size=18,color="black"))+
  ylim(c(0.4,3.5))
  
