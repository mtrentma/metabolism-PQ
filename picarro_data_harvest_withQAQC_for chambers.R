# Picarro data handling attempt.  Stealing from https://github.com/bpbond/R-data-picarro/blob/master/scripts/picarro.R
library(plyr)
library(dplyr)
library(readr)     
library(ggplot2)
library(cowplot)
library(reshape2)
library(xts)
library(dygraphs)
library(tidyr)

## Prior to running the rest of the code
## 1) Set the wd to folder location of the raw Picarro dat files and meta data
## 2) Specify meta file name

setwd("C:/Users/matt/IDrive-Sync/Postdoc/Chamber metabolism/Chambers2021")
meta_file_name <- "2021_07_19_co2.txt"
## 3) Specify start time of Picarro run (using Picarro time)
picarro_start_time <- "2021-07-20 17:00:00"
## 4) Specify notes file name
notes_file_name <- "2021_07_19_co2 notes.csv"


###########
## IMPORT
###########

## Import meta data for each day:  In this case metadata consists of a sample identifier, and time stamp, where time stamp is when the sample syring hits about 15mL and it 
## looks like the trace is leveling off for that sample.  This metadata can and should have other stuff in it, such as sample location, collection time,
## pH, volume of sample, volume of equlibration gas (typically 70mL), temp of equilibration
meta <- read.table(meta_file_name, sep=",",header=T)
meta$DATETIME <- as.POSIXct(as.character(meta$ResetTime), format="%Y-%m-%d %H:%M:%S")
meta$samplenum<-seq(1:length(meta$ResetTime))

## Import the raw data
filelist <-list.files( pattern="dat", full.names = TRUE)  ##makes a vector of that day's files
rawdata <- list()  #empty list
#below loops over all the files making a list containg a tibble for each file and stuffing it into a list
for(f in filelist) {
  cat("Reading", f, "\n") #no idea why we need this
  read_table(f) %>% 
    select(DATE, TIME, ALARM_STATUS, MPVPosition, `12CO2_dry`, Delta_Raw_iCO2, HP_12CH4, `HP_Delta_iCH4_Raw`) -> rawdata[[f]] #add whichever things we want here
}
rawdata <- bind_rows(rawdata)

## Fix datetime
rawdata$DATETIME <- lubridate::ymd_hms(paste(rawdata$DATE, rawdata$TIME))
rawdata$DATETIME <- as.POSIXct(as.character(rawdata$DATETIME), format="%Y-%m-%d %H:%M:%S")

# Get rid of unneeded fields
rawdata$DATE <- rawdata$TIME <- NULL

##average data for each second where we have a measure.  Not optional
rawdata<- rawdata %>% group_by(as.numeric(rawdata$DATETIME)) %>% summarise_all(funs(mean))

##get rid of all the early in the day data to make a smaller file.  This step is optional
rawdata <- subset(rawdata, rawdata$DATETIME > picarro_start_time)

##if ok, then name the rawdat file for that day here
dat<-rawdata
colnames(dat)[4:7] <- c("CO2","Delta_iCO2","CH4","Delta_iCH4")

#Plot it. Look ok?
plot(dat$DATETIME,dat$CO2)
plot(dat$DATETIME,rawdata$Delta_Raw_iCO2)
plot(dat$DATETIME,rawdata$HP_12CH4)
plot(dat$DATETIME,rawdata$`HP_Delta_iCH4_Raw`)


#############
## SUBSET
#############
## Subset the main datafile based on sample time + some amount of time, 
## 40 seconds appears correct, but always check this time.  
## To short and we lost data.  Too long and we average off the plateau
meta$ENDTIME<-meta$DATETIME+30

subdat<-list()
for(i in 1:length(meta$samplenum)){
  
  subdat[[i]]<-subset(dat, DATETIME>meta$DATETIME[i]  & DATETIME< meta$ENDTIME[i])
  subdat[[i]]$samplenum <-i
}
subdat <- bind_rows(subdat)

## Plot All Together
subdat_melt <- melt(subdat[,4:9], id.vars=c("DATETIME","samplenum"))
ggplot(subdat_melt, aes(DATETIME, value, color=samplenum))+
  geom_point()+
  facet_wrap(~variable, ncol=2, scales = "free")

######################
## QA/QC every sample
######################
## Look at time series traces plotted over original data
## and decide on new Reset time if necessary

vis <- function(d, subd, v){
  
  d <- d[,c("DATETIME",v)]
  subd <- subd[,c("DATETIME",v,"samplenum")]
  comb <- merge(d, subd, by="DATETIME", all.x=T)
  colnames(comb) <- c("DateTime","All","Sub","samplenum")
  comb.xts <- as.xts(comb, order.by=comb$DateTime)
  
  dygraph(comb.xts[,c("All","Sub","samplenum")]) %>%
    dyAxis("y2", label = "samplenum", independentTicks = TRUE)%>%
    dySeries("Sub", fillGraph = T)%>%
    dySeries("samplenum", axis=('y2'))%>%
    dyRangeSelector()
}

## 13C-CO2
vis(dat, subdat, "Delta_iCO2")

## CO2
vis(dat, subdat, "CO2")

## CH4
vis(dat, subdat, "CH4")

## 13C-CH4
vis(dat, subdat, "Delta_iCH4")


#####################################################
## Adjust overall start time and individual sites
####################################################

## This section is sample run specific
## Save this section of code in text file for future reference

## Oak 2018-09-28

## Samples that need later start times
#meta[which(meta$samplenum %in% c(4,6,9,11,14,16,70)),]$DATETIME <- meta[which(meta$samplenum %in% c(4,6,9,11,14,16,70)),]$DATETIME + 5

## Samples that need earlier start times
#meta[which(meta$samplenum %in% c(49,83)),]$DATETIME <- meta[which(meta$samplenum %in% c(49,83)),]$DATETIME - 5

## Averaging could last 5 seconds longer
#meta$ENDTIME<-meta$DATETIME+35

#subdat_rev<-list()
#for(i in 1:length(meta$samplenum)){
  
 # subdat_rev[[i]]<-subset(dat, DATETIME>meta$DATETIME[i]  & DATETIME< meta$ENDTIME[i])
 #  subdat_rev[[i]]$samplenum <-i
#}
#subdat_rev <- bind_rows(subdat_rev)


## 13C-CO2
#vis(dat, subdat_rev, "Delta_Raw_iCO2")
## CO2
#vis(dat, subdat_rev, "12CO2_dry")
## CH4
#vis(dat, subdat_rev, "CH4")
## 13C-CH4
#vis(dat, subdat, "Delta_iCH4")

## Get rid of samples #10, 50, 51
#subdat_rev <- subdat_rev[-which(subdat_rev$samplenum %in% c(10,50,51)),]

#########################################################################

#Summarize per sample
subdatsum<- subdat %>% group_by(samplenum) %>% summarise(CO2=mean(CO2), delCO2=mean(Delta_iCO2))#, 
                                                         #CH4=mean(CH4), delCH4=mean(Delta_iCH4) )

###########################################
##merge with meta, and give it a new name
#############################################
#Thinking about the metafile, here is what should also go in
#water_temp, equil_temp, pH, ANC, air_vol (it won't always be 70), baro_press (absolute, mm Hg)
#not the  kH_CO2 because that is a derived value that we will calc each time
##Load notes
notes<-read.csv(notes_file_name)
notes$date.time<-as.POSIXct(paste(notes$date, notes$time),format="%m-%d-%Y %H:%M:%S")


##Join metadata and notes to data
sum_file<- full_join(meta,subdatsum)  #there you go, joining by samplenum
sum_file<- full_join(sum_file,notes)  #there you go, joining by Syringe

##Remove standards and no-air
sum_file<-sum_file[complete.cases(sum_file), ]


#View(sum_file)
#sum_file <- sum_file[-which(sum_file$Syringe %in% c('Std_Tank1','Std_Tank2','no_air')),]
#meta <- meta[-which(meta$Syringe %in% c('Std_Tank1','Std_Tank2','no_air')),]

## Plot
sum_file_melt <- melt(sum_file[,c("CO2","delCO2","DATETIME","samplenum")], id.vars=c("DATETIME","samplenum"))
ggplot(sum_file_melt, aes(DATETIME, value, color=samplenum))+
  geom_point()+
  facet_wrap(~variable, ncol=2, scales = "free")




#CORRECT CO2 FOR HEADSPACE EQUILIBRATION


###########################################
##Calculate pCO2 of water from pCO2 headspace
###########################################

##### SOLUBILITY CONSTANT FOR CO2 (KH) #####
# Use Henry's Law and parameters from Weiss 1974 to calculate KH in units of mol L-1 atm-1 as a function of temperature (also in Demarty et al 2011; many others)
KH.CO2 <- function(temp_equil){
  tempK <- temp_equil + 273.15
  KH.CO2 <- exp( -58.0931 + (90.5069*(100/tempK)) +22.2940*log((tempK/100), base=exp(1)) )
  KH.CO2
}

#####ppmv to mol/L
co2 <- function(temp_samp,CO2_post ){
  KH.samp <- KH.CO2(temp_samp)
  co2.mol.L<-  KH.samp*CO2_post/1000000
}


###########################################
##Calculate DIC of water from pCO2 of water
###########################################

## Functions to determine equilibrium constants dependent on temperature of water when sampled
#From Milerno et al 2006 "Dissociation constants of carbonic acid in seawater as a function of salinity and temperature"
K1calc<- function(temp_equil) { 10^-(-126.34048+6320.813/(temp_equil+273.15)+19.568224*log(temp_equil+273.15))}
K2calc<- function(temp_equil) { 10^-(-90.18333+5143.692/(temp_equil+273.15)+14.613358*log(temp_equil+273.15))}


### Function to solve carbonate chemistry using CO2 & Alkalinity
## See R Markdown file (XXXXXXX) for derivation of equations
# H=H+ 
# A=alkalinity in mol/L
# StmpCO2.umol=CO2 in mol/L 
# B=bicarbonate 
# Ca=carbonate 
# D=DIC in umol/L
# A_check=alkalinity check
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

### Function to calculate correction of perturbed sample for DIC before equilibration
## Originally from Dickenson et al. 2007--Chapter 4, highlighted again by Koschorreck et al. for freshwater

DIC_correction<-function(CO2_pre,CO2_post,vol_hs,temp_equil,vol_samp){
  delta_DIC<-(((CO2_post-CO2_pre)/1000000))/(R*(temp_equil + 273.15))*(vol_hs/vol_samp)
  delta_DIC
}

### Apply correction to DIC data from Carb1####

DIC_corr<-Carb1$D+delta_DIC

### Recalculate CO2 using DIC and Alkalinity####

Carbfrom_D_A <- function(K1, K2, DIC_corr, A){
  a <- A
  b <- K1*(A-DIC_corr)
  c <- (A-(2*DIC_corr))*K1*K2
  H_t <- ((-1*b)+sqrt((b^2)-(4*a*c)))/(2*a)
  
  pH_t <- -1*log10(H_t)
  
  B_t <- (DIC_corr*K1*H_t)/((H_t^2)+(K1*H_t)+(K1*K2))
  Ca_t <- (DIC_corr*K1*K2)/((H_t^2)+(K1*H_t)+(K1*K2))
  
  C_t <- (H_t*B_t)/K1
  D2 <- C_t + B_t + Ca_t
  A2_check <- B_t + 2*Ca_t
  D_check <- C_t + B_t + Ca_t
  Carb2 <- list(H_t, pH_t, B_t, Ca_t, C_t, D_check, A2_check)
  names(Carb2) <- c("H", "pH", "B", "Ca", "CO2", "D_check", "A check")
  return(Carb2)
}


###Calculate CO2 in water
R= 0.08205601 # gas constant in L*atm/mol*K
dens=0.9998395 # density of freshwater

#Input variables
Syringe=sum_file$Syringe
temp_equil=sum_file$temp_equil        # temperature of water immediately after equilibriaum in C
temp_samp=sum_file$temp_samp        # temperature of water at the time of sampling in C
press_samp=sum_file$press_samp        # pressure at the time of sampling in atm
vol_hs=sum_file$vol_hs            # volume of headspace in mL
vol_samp=sum_file$vol_samp                    # volume of water sample in mL
CO2_pre=sum_file$CO2_pre          # pCO2 of headspace before equilibrium (zero if using zero air) # uatm or ppmv 
CO2_post=sum_file$CO2          # pCO2 of hs after equilibrium  # uatm or ppmv
A=sum_file$alkalinity/100/1000             # alkalinity in mol/L



Calculate_CO2<-function(temp_equil, temp_samp, press_samp, vol_hs, vol_samp, CO2_pre,CO2_post, A){
  StmpCO2.mol<-co2(temp_samp,CO2_post )
  K1<-K1calc(temp_equil)
  K2<-K2calc(temp_equil)
  Carb1<-Carbfrom_C_A(K1, K2, StmpCO2.mol , A)
  delta_DIC<-DIC_correction(CO2_pre,CO2_post,vol_hs,temp_equil,vol_samp)
  DIC_corr<-Carb1$D+delta_DIC
  Carb2<-Carbfrom_D_A(K1, K2, DIC_corr, A)
  Result <- as.data.frame(cbind(sum_file$Syringe,Carb2$pH, CO2_post,Carb2$CO2,Carb2$CO2/(KH.CO2(temp_equil)/1000000),DIC_corr ))
  names(Result) <- c("Syringe", "pH","CO2.ppmv.pre", "CO2.mol","CO2.ppmv", "DIC.mol")
  return(Result)
}


##Returns a list with the corrected DIC (D) and CO2 
result<-Calculate_CO2(temp_equil, temp_samp, press_samp, vol_hs, vol_samp, CO2_pre,CO2_post,A)

##Reformat
result$Syringe<-as.factor(result$Syringe)
result$pH<-as.numeric(result$pH)
result$CO2.ppmv.pre<-as.numeric(result$CO2.ppmv.pre)
result$CO2.mol<-as.numeric(result$CO2.mol)
result$CO2.ppmv<-as.numeric(result$CO2.ppmv)
result$DIC.mol<-as.numeric(result$DIC.mol)


##Join with previous data and reformat date
sum_file_final<- full_join(sum_file,result) 
sum_file_final$date.time <- as.POSIXct(as.character(sum_file_final$date.time), format="%m/%d/%Y %H:%M")

##Plot data in time
plot(sum_file_final$date.time,sum_file_final$CO2.ppmv)
plot(sum_file_final$date.time,sum_file_final$DIC.mol)
plot(sum_file_final$date.time,sum_file_final$delCO2)



sum_file_export<-sum_file_final[,c(1:3,9:22)]


############
## Export
###########
write.csv(sum_file_export, "2021_07_19_CO2_QAQC.csv",row.names = FALSE)


###################STOP########################################################################################################
###############################################################################################################################


Sample.ID=sum_file$Syringe
HS.mCO2.before=0
HS.mCO2.after=sum_file$CO2 
Temp.insitu=sum_file$temp_samp 
Temp.equil=sum_file$temp_equil
Alkalinity.measured=(sum_file$alkalinity*1000/100*100/50) 
Volume.gas=sum_file$vol_hs 
Volume.water=sum_file$vol_samp
Bar.pressure=(sum_file$press_samp*1.013) 
Constants=1
Salinity=0

dataset<-as.data.frame(cbind(Sample.ID,HS.mCO2.before,HS.mCO2.after,Temp.insitu,Temp.equil,Alkalinity.measured,
                             Volume.gas,Volume.water,Bar.pressure,Constants,Salinity))
dataset$Alkalinity.measured<-as.numeric(dataset$Alkalinity.measured)
dataset$HS.mCO2.before<-as.numeric(dataset$HS.mCO2.before)
dataset$HS.mCO2.after<-as.numeric(dataset$HS.mCO2.after)
dataset$Temp.equil<-as.numeric(dataset$Temp.equil)
dataset$Temp.insitu<-as.numeric(dataset$Temp.insitu)
dataset$Volume.gas<-as.numeric(dataset$Volume.gas)
dataset$Volume.water<-as.numeric(dataset$Volume.water)
dataset$Bar.pressure<-as.numeric(dataset$Bar.pressure)
dataset$Constants<-as.numeric(dataset$Constants)
dataset$Salinity<-as.numeric(dataset$Salinity)


##LOAD FUNCTION

Rheadspace <-  function(...){
  arguments <- list(...)
  
  # test arguments and initialize variables
  if (is.data.frame(arguments[[1]])) {
    input.table=arguments[[1]]
    if (dim(input.table)[2]!=11){
      stop("You should input a data frame with 11 columns. See the readme file or comments in the function", call.=FALSE)
    }else{
      Sample.ID = as.character(input.table$Sample.ID)
      mCO2_headspace = input.table$HS.mCO2.before #the mCO2 (ppmv) of the headspace "before" equilibration
      mCO2_eq = input.table$HS.mCO2.after #the measured mCO2 (ppmv) of the headspace "after" equilibration
      temp_insitu = input.table$Temp.insitu #in situ water temperature in degrees celsius
      temp_eq = input.table$Temp.equil #the water temperature after equilibration in degree celsius
      alk = input.table$Alkalinity.measured #Total alkalinity (micro eq/L) of the water sample
      vol_gas = input.table$Volume.gas #Volume of gas in the headspace vessel (mL)
      vol_water = input.table$Volume.water #Volume of water in the headspace vessel (mL)   
      Bar.pressure = input.table$Bar.pressure #Barometric pressure at field conditions in kPa. 101.325 kPa = 1 atm   
      c_constants = input.table$Constants #Constants for carbonate equilibrium (1=Freshwater; 2=Estuarine; 3=Marine) 
      Salinity = input.table$Salinity #Salinity in PSU. Set to zero if Constants = 1
    } 
  } else if (length(arguments)==11) {
    Sample.ID = as.character(arguments[[1]])
    mCO2_headspace = arguments[[2]] #the mCO2 (ppmv) of the headspace "before" equilibration
    mCO2_eq = arguments[[3]] #the measured mCO2 (ppmv) of the headspace "after" equilibration
    temp_insitu = arguments[[4]] #in situ water temperature in degrees celsius
    temp_eq = arguments[[5]] #the water temperature after equilibration in degree celsius
    alk = arguments[[6]] #Total alkalinity (micro eq/L) of the water sample
    vol_gas = arguments[[7]] #Volume of gas in the headspace vessel (mL)
    vol_water = arguments[[8]] #Volume of water in the headspace vessel (mL)   
    Bar.pressure = arguments[[9]] #Barometric pressure at field conditions in kPa. 101.325 kPa = 1 atm   
    c_constants = arguments[[10]] #Constants for carbonate equilibrium (1=Freshwater; 2=Estuarine; 3=Marine) 
    Salinity = arguments[[11]] #Salinity in PSU. Set to zero if Constants = 1
  } else {
    stop("You should input either a data frame or a vector of 11 values. See the readme file or comments in the function", call.=FALSE)
  }
  
  
  
  #initialization of variables
  pCO2_orig <- data.frame(matrix(NA,length(mCO2_headspace),9))
  names(pCO2_orig) <- c("Sample.ID","mCO2 complete headspace (ppmv)","pCO2 complete headspace (micro-atm)", "CO2 concentration complete headspace (micro-mol/L)", "pH", "mCO2 simple headspace (ppmv)", "pCO2 simple headspace (micro-atm)", "CO2 concentration simple headspace (micro-mol/L)", "% error")
  
  R <- 0.082057338 #L atm K-1 mol-1
  
  #the function uniroot cannot handle vectors, so we need a loop
  for (i in 1:length(mCO2_headspace)){ 
    
    AT = alk[i]*(1e-6) #conversion to mol/L
    
    #Constants of the carbonate ewuilibrium
    # Kw = the dissociation constant of H2O into H+ and OH-
    # Kh = the solubility of CO2 in water - equilibration conditions
    # Kh2 = the solubility of CO2 in water - in situ field conditions
    # K1 = the equilibrium constant between CO2 and HCO3-
    # K2 = the equilibrium constant between HCO3- and CO3 2-
    
    # Solubility coefficients from Weiss (1974) with Sal=0 for freshwater option
    # Dissociation of water from Dickson and Riley (1979)
    
    if (c_constants == 1) {
      
      #Millero, F. (1979). The thermodynamics of the carbonate system in seawater
      #Geochimica et Cosmochimica Acta 43(10), 1651 1661.  
      K1=10^-(-126.34048+6320.813/(temp_eq[i]+273.15)+19.568224*log(temp_eq[i]+273.15))
      K2=10^-(-90.18333+5143.692/(temp_eq[i]+273.15)+14.613358*log(temp_eq[i]+273.15))
      
      Kw = exp(148.9652-13847.26/(temp_eq[i]+273.15)-23.6521*log(273.15+temp_eq[i]))
      Kh = 10^((-60.2409+93.4517*(100/(273.15+temp_eq[i]))+23.3585*log((273.15+temp_eq[i])/100))/log(10)) # mol/L/atm equilibration conditions
      Kh2 = 10^((-60.2409+93.4517*(100/(273.15+temp_insitu[i]))+23.3585*log((273.15+temp_insitu[i])/100))/log(10)) # mol/L/atm original conditions
      
      
    } else if (c_constants == 2) {
      
      #Millero, F. (2010). Carbonate constants for estuarine waters Marine and Freshwater
      #Research 61(2), 139. As amended by Orr et al. 2015.
      pK10=(-126.34048+6320.813/(temp_eq[i]+273.15)+19.568224*log(temp_eq[i]+273.15))
      A1 = 13.4038*Salinity[i]^0.5 + 0.03206*Salinity[i] - 5.242e-5*Salinity[i]^2
      B1 = -530.659*Salinity[i]^0.5 - 5.8210*Salinity[i]
      C1 = -2.0664*Salinity[i]^0.5
      pK1 = pK10 + A1 + B1/(temp_eq[i]+273.15) + C1*log(temp_eq[i]+273.15)
      K1 = 10^-pK1;
      pK20=(-90.18333+5143.692/(temp_eq[i]+273.15)+14.613358*log(temp_eq[i]+273.15))
      A2 = 21.3728*Salinity[i]^0.5 + 0.1218*Salinity[i] - 3.688e-4*Salinity[i]^2
      B2 = -788.289*Salinity[i]^0.5 - 19.189*Salinity[i]
      C2 = -3.374*Salinity[i]^0.5
      pK2 = pK20 + A2 + B2/(temp_eq[i]+273.15) + C2*log(temp_eq[i]+273.15)
      K2 = 10^-pK2;
      
      Kw=exp(148.9652-13847.26/(temp_eq[i]+273.15)-23.6521*log(273.15+temp_eq[i])+sqrt(Salinity[i])*(118.67/(temp_eq[i]+273.15)-5.977+1.0495*log(273.15+temp_eq[i]))-0.01615*Salinity[i])
      Kh = exp(-60.2409 + 93.4517 * (100 / (273.15 + temp_eq[i])) + 23.3585 * log( (273.15 + temp_eq[i]) / 100 )+Salinity[i]*(0.023517-0.023656*(273.15+temp_eq[i])/100+0.0047036*((273.15+temp_eq[i])/100)^2))
      Kh2 = exp(-60.2409 + 93.4517 * (100 / (273.15 + temp_insitu[i])) + 23.3585 * log( (273.15 + temp_insitu[i]) / 100 )+Salinity[i]*(0.023517-0.023656*(273.15+temp_insitu[i])/100+0.0047036*((273.15+temp_insitu[i])/100)^2))
      
    } else if (c_constants == 3) {
      
      #Dickson, A. G., Sabine, C. L., and Christian, J. R. (2007): Guide to best practices for
      #ocean CO2 measurements, PICES Special Publication 3, 191 pp.
      K1=10^(-3633.86/ (temp_eq[i] + 273.15)+61.2172-9.67770*log(temp_eq[i]+273.15)+0.011555*Salinity[i]-0.0001152*Salinity[i]^2)
      K2 = 10^(-417.78/ (temp_eq[i] + 273.15) - 25.9290 + 3.16967*log(temp_eq[i]+273.15)+0.01781*Salinity[i]-0.0001112*Salinity[i]^2)
      
      Kw=exp(148.9652-13847.26/(temp_eq[i]+273.15)-23.6521*log(273.15+temp_eq[i])+sqrt(Salinity[i])*(118.67/(temp_eq[i]+273.15)-5.977+1.0495*log(273.15+temp_eq[i]))-0.01615*Salinity[i])
      Kh = exp(-60.2409 + 93.4517 * (100 / (273.15 + temp_eq[i])) + 23.3585 * log( (273.15 + temp_eq[i]) / 100 )+Salinity[i]*(0.023517-0.023656*(273.15+temp_eq[i])/100+0.0047036*((273.15+temp_eq[i])/100)^2))
      Kh2 = exp(-60.2409 + 93.4517 * (100 / (273.15 + temp_insitu[i])) + 23.3585 * log( (273.15 + temp_insitu[i]) / 100 )+Salinity[i]*(0.023517-0.023656*(273.15+temp_insitu[i])/100+0.0047036*((273.15+temp_insitu[i])/100)^2))
      
    } else {
      print(i)
      stop("Option for carbonate equilibrium constants should be a number between 1 and 3", call.=FALSE)
      
    }
    
    HS.ratio <- vol_gas[i]/vol_water[i] #Headspace ratio (=vol of gas/vol of water)
    
    #The following calculations assume 1 atm, this is corrected later for measured pressure in the field.
    
    #DIC at equilibrium
    co2 <- Kh * mCO2_eq[i]/1000000
    h_all <- polyroot(c(-(2*K1*K2*co2),-(co2*K1+Kw),AT,1))
    real<-Re(h_all)
    h <-real[which(real>0)]
    DIC_eq <- co2 * (1 + K1/h + K1 * K2/(h * h))
    
    #DIC in the original sample
    DIC_ori <- DIC_eq + (mCO2_eq[i] - mCO2_headspace[i])/1000000/(R*(temp_eq[i]+273.15))*HS.ratio
    
    #pCO2 in the original sample
    h_all <- polyroot(c(-(K1*K2*Kw),K1*K2*AT-K1*Kw-2*DIC_ori*K1*K2,AT*K1-Kw+K1*K2-DIC_ori*K1,AT+K1,1))
    real<-Re(h_all)
    h <-real[which(real>0)]
    
    co2 <- h* (DIC_ori * h * K1/(h * h + K1 * h + K1 * K2)) / K1
    
    pCO2_orig[i,1] <- as.character(Sample.ID[i])
    pCO2_orig[i,2] <- co2/Kh2*1000000
    pCO2_orig[i,3] <- pCO2_orig[i,2]*Bar.pressure[i]/101.325
    pCO2_orig[i,4] <- co2*1000000
    pCO2_orig[i,5] <- -log10( h )
    
    
    #Calculation not accounting for alkalinity effects and associated error
    
    #concentration and total mass in the water sample assuming ideal gas from the pCO2 measured at the headspace
    CO2_solution <- mCO2_eq[i]/1000000*Kh #mol/L
    CO2_solution_mass <- CO2_solution * vol_water[i]/1000 #mol
    
    #mass of CO2 in the measured headspace
    final_C_headspace_mass <- mCO2_eq[i]/1000000*(vol_gas[i]/1000) / (R * (temp_eq[i]+273.15)) #mol
    
    mols_headspace <- mCO2_headspace[i]/1000000*(vol_gas[i]/1000)/(R * (temp_eq[i]+273.15)) #mol PV / RT = n
    
    #implication: mass, concentration, and partial pressure of CO2 in the original sample (aount in sample and headspace after equilibration minus original mass in the headspace)
    Sample_CO2_mass <- CO2_solution_mass + final_C_headspace_mass - mols_headspace #mol
    Sample_CO2_conc <- Sample_CO2_mass/(vol_water[i]/1000) #mol/L
    pCO2_orig[i,6] <- Sample_CO2_conc/Kh2*1000000 #ppmv
    pCO2_orig[i,7] <- pCO2_orig[i,6]*Bar.pressure[i]/101.325 # micro-atm
    pCO2_orig[i,8] <- Sample_CO2_conc*1000000 # micro-mol/L
    #calculation of the error
    pCO2_orig[i,9] <- (pCO2_orig[i,6]-pCO2_orig[i,2])/pCO2_orig[i,2] *100  #%
  }
  
  
  return(pCO2_orig) #Output data frame
  
}


##Calculate DIC
pCO2 <- Rheadspace(dataset)

