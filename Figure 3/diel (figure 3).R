

# Load Libraries ----------------------------------------------------------

library(chron) 
library(ggplot2)
library(patchwork)
library(reshape2); library(dplyr); library(tidyr)




# Simulate time -----------------------------------------------------------


Sys.setenv(TZ='GMT')

time<-seq(from = as.numeric(chron(dates="05/21/22", times="02:00:00")),
          to=as.numeric(chron(dates="05/22/22", times="04:00:00")), by=10/1440 )
dtime<- chron(time)



# Simulate light data -----------------------------------------------------

## From Yard et al. (1995) Ecological Modelling.  Remember your trig?  
## calculate light as umol photon m-2 s-1.
## Arguments are:  
## time = a date and time input (i.e. a chron object) 
## lat = latitude of field site
## longobs = longitude of field site
## longstd = standard longitude of the field site (NOTE: watch daylight savings time!!!). For PST, longstd is be 120 degrees. But during PDT it is 105 degrees. MST is 105 deg. MDT is 90. 
## year = the year for which you collected data and is entered as "2013-01-01"


## convert degrees to radians
radi<-function(degrees){(degrees*pi/180)}

## function to estimate light
lightest<- function (time, lat, longobs, longstd, year ) {
  
  jday<-as.numeric(trunc(time)-as.numeric(as.Date(year)))
  E<- 9.87*sin(radi((720*(jday-81))/365)) - 7.53*cos(radi((360*(jday-81))/365)) - 1.5*sin(radi((360*(jday-81))/365))
  LST<-as.numeric (time-trunc(time))
  ST<-LST+(3.989/1440)*(longstd-longobs)+E/1440
  solardel<- 23.439*sin(radi(360*((283+jday)/365)))
  hourangle<-(0.5-ST)*360
  theta<- acos( sin(radi(solardel)) * sin(radi(lat)) + cos(radi(solardel)) * cos(radi(lat)) * cos(radi(hourangle)) )
  suncos<-ifelse(cos(theta)<0, 0, cos(theta))
  GI<- suncos*2326
  GI	
  
}
#end of function

## Use function to simulate light
light<- lightest(time=dtime, lat=47,  longobs=114, longstd= 105, year="2022-01-01")


# Simulate temperature ----------------------------------------------------
temp<- rep(15, length(light))  #too cold to swim in



# Function to simulate gas-exchange ---------------------------------------------------


Kcor<-function (temp,K600) {
  K600/(600/(1800.6-(temp*120.1)+(3.7818*temp^2)-(0.047608*temp^3)))^-0.5
}
#end of function

# Function to simulate oxygen concentration

osat<- function(temp, bp) {
  
  tstd<-log((298.15-temp) / (273.15 + temp))
  
  a0<-2.00907
  a1<-3.22014
  a2<-4.0501
  a3<-4.94457
  a4<- -0.256847
  a5<- 3.88767
  
  u<-10^(8.10765-(1750.286/(235+temp)))
  
  sato<-(exp(a0 + a1*tstd + a2*tstd^2 + a3*tstd^3 + a4*tstd^4+a5*tstd^5))*((bp-u)/(760-u))*1.42905
  sato
}

##end of function


# Simulate oxygen ------------------------------------------

# Function simulate low GPP oxygen
onestationsim<-function(GPP, ER, z, temp, K, light, bp, ts) {
  
  oxy.mod<-numeric(length(light))
  oxy.mod[1]<-9.97  #trial end error this first point based on the O2 data at the end of the time series
  
  # this is the metabolism equation as in Van de Bogert et al (2007) L&OMethods
  for (i in 2:length(light)) { oxy.mod[i]<-oxy.mod[i-1]+((GPP/z)*(light[i]/sum(light)))+ ER*ts/z+(Kcor(temp[i],K))*ts*(osat(temp[i],bp)-oxy.mod[i-1]) }
  oxy.mod
}
# end of function

# SCENARIO 1: Low GPP and low NO3
## Known (simulated) parameters. These values were picked to create a scenario where % O2sat stayed at or below %100
## All non-metabolic variables are held constant except light(depth, gas-exchange, pressure, temperature).
lowgpp<-0.6
lowER<--1.5# scales with GPP
depth=1 #constant
K=15 #constant
bp=760 #constant

sim_oxy<- onestationsim(GPP=lowgpp, ER=lowER,z=depth,temp=temp,light=light, K=K, bp=bp, ts=10/1440) #Run simulation function with given parameters
sim_oxy_sat<-100*sim_oxy/osat(temp, 760)# transform from concentration to O2 sat
data_lowgpp<-data.frame(dtime, sim_oxy_sat) #format
data_lowgpp$dtime<-as.POSIXct(data_lowgpp$dtime) #format
data_lowgpp$time<-as.POSIXct(data_lowgpp$dtime, format="%H:%M:%S") #format


## Plot low GPP and low NO3 % oxygen data
p1<-
  ggplot(data_lowgpp, aes(x=time, y=sim_oxy_sat))+
  geom_rect(aes( xmin=time[20], xmax=time[111], ymin=-Inf, ymax=Inf, fill="Daytime"), color=NA, alpha=0.2)+
  geom_line(size=1)+
  theme_classic()+
  ylab(expression(paste(O[2], " saturation (%)")))+
  xlab("Time of Day")+
  theme(axis.title.x=element_text(size=10,colour = "black"))+
  theme(axis.title.y=element_text(size=10,colour = "black"))+
  theme(axis.text.y=element_text(size=10,colour = "black"))+
  theme(axis.text.x=element_text(size=8,colour = "black"))+
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "6 hours", limits=c(as.POSIXct("2022-05-21 02:00:00"), as.POSIXct("2022-05-22 02:00:00")))+
  scale_y_continuous(limits=c(93,110), breaks=seq(from=80, to=130, by=5))+
  geom_hline(yintercept=100,linetype="dashed", size=1)+
  theme(axis.title.x=element_blank())+
  theme(plot.title = element_text(size = 12))+
  scale_fill_manual('',
                    values = 'lightblue',  
                    guide = guide_legend(override.aes = list(alpha = 0.2)))+
  theme(legend.text=element_text(size=12),legend.position="none")
  
  
  

#SCENARIO 2: High GPP and low NO3

# Function simulate high GPP oxygen
onestationsim<-function(GPP, ER, z, temp, K, light, bp, ts) {
  
  oxy.mod<-numeric(length(light))
  oxy.mod[1]<-9.51  #trial end error this first point based on the O2 data at the end of the time series
  
  # this is the metabolism equation as in Van de Bogert et al (2007) L&OMethods
  for (i in 2:length(light)) { oxy.mod[i]<-oxy.mod[i-1]+((GPP/z)*(light[i]/sum(light)))+ ER*ts/z+(Kcor(temp[i],K))*ts*(osat(temp[i],bp)-oxy.mod[i-1]) }
  oxy.mod
}
# end of function

##parameters. These values were picked to have an O2 saturation curve that goes above 110. Only GPP and ER are different from the first scenario.
highgpp<-8
highER<--8# scales with GPP

sim_oxy<- onestationsim(GPP=highgpp, ER=highER,z=depth,temp=temp,light=light, K=K, bp=bp, ts=10/1440)
sim_oxy_sat<-100*sim_oxy/osat(temp, 760)# transform from concentration to O2 sat
data_highgpp<-data.frame(dtime, sim_oxy_sat)# format
data_highgpp$dtime<-as.POSIXct(data_lowgpp$dtime) # format
data_highgpp$time<-as.POSIXct(data_lowgpp$dtime, format="%H:%M:%S")# format


## Plot high GPP and low NO3 % oxygen data
p2<-ggplot(data_highgpp, aes(x=time, y=sim_oxy_sat))+
  geom_rect(aes( xmin=time[20], xmax=time[111], ymin=-Inf, ymax=Inf, fill="Daytime"), color=NA, alpha=0.2)+
  geom_line(size=1)+
  theme_classic()+
  ylab(expression(paste(O[2], " saturation (%)")))+
  xlab("Time of Day")+
  theme(axis.title.x=element_text(size=10,colour = "black"))+
  theme(axis.title.y=element_text(size=10,colour = "black"))+
  theme(axis.text.y=element_text(size=10,colour = "black"))+
  theme(axis.text.x=element_text(size=8,colour = "black"))+
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "6 hours", limits=c(as.POSIXct("2022-05-21 02:00:00"), as.POSIXct("2022-05-22 02:00:00")))+
  scale_y_continuous(limits=c(93,110), breaks=seq(from=80, to=130, by=5))+
  geom_hline(yintercept=100,linetype="dashed", size=1)+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(plot.title = element_text(size = 12))+
  scale_fill_manual('',
                    values = 'lightblue',  
                    guide = guide_legend(override.aes = list(alpha = 0.2)))+
  theme(legend.text=element_text(size=12),legend.position="none")


#SCENARIO 3: High GPP and high NO3
## Plot high GPP and high NO3 % oxygen data. Same as plot 2, as in this scenario, N availability does not affect production.
p3<-ggplot(data_highgpp, aes(x=time, y=sim_oxy_sat))+
  geom_rect(aes( xmin=time[20], xmax=time[111], ymin=-Inf, ymax=Inf, fill="Daytime"), color=NA, alpha=0.2)+
  geom_line(size=1)+
  theme_classic()+
  ylab(expression(paste(O[2], " saturation (%)")))+
  xlab("Time of Day")+
  theme(axis.title.x=element_text(size=10,colour = "black"))+
  theme(axis.title.y=element_text(size=10,colour = "black"))+
  theme(axis.text.y=element_text(size=10,colour = "black"))+
  theme(axis.text.x=element_text(size=8,colour = "black"))+
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "6 hours", limits=c(as.POSIXct("2022-05-21 02:00:00"), as.POSIXct("2022-05-22 02:00:00")))+
  scale_y_continuous(limits=c(93,110), breaks=seq(from=80, to=130, by=5))+
  geom_hline(yintercept=100,linetype="dashed", size=1)+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(plot.title = element_text(size = 12))+
  scale_fill_manual('',
                    values = 'lightblue',  
                    guide = guide_legend(override.aes = list(alpha = 1)))+
  theme(legend.position="right",legend.text=element_text(size=9), legend.title=element_text(size=10), legend.text.align = 0)
  

# Simulate PQ data based on simulated oxygen data -------------------------


# SCENARIO 1: Low GPP and low NO3

#PQ_C is the carbon portion that describes the type of organic molecule synthesized. 
#For phytoplakton this value is believed to be 1.1. 
pq_c<-1.1

#PQ_N is the NO3 portion. NO3 is low in scenario 1, so this is near zero
pq_n<-0

#PQ_pr is the photorespiration portion. Low GPP below 100% sat means no PR
pq_pr<-0

#PQ_net is all combined
pq_net<-pq_c+pq_n+pq_pr

#Format for plotting
pq_lowgpp<-melt(data.frame(as.POSIXct(dtime,format="%H:%M:%S"),pq_c, pq_n, pq_pr,pq_net ),id="as.POSIXct.dtime..format.....H..M..S..")
names(pq_lowgpp)<-c("time","factor","pq_x" )
labs<-c(expression(phantom("--")*PQ['Bio']),
        expression(PQ['C']*phantom("-")),
        expression(PQ['NO_3^-']*phantom("-")),
        expression(PQ['Pr']))
pq_lowgpp$factor<-factor(pq_lowgpp$factor, levels= c("pq_net","pq_c", "pq_n", "pq_pr"))


#Plot PQ for SCENARIO 1
p4<-ggplot(pq_lowgpp, aes(x=time, y=pq_x, color=factor(factor)))+
  geom_line(aes(linetype=factor),size=1)+
  theme_classic()+
  ylab(expression(paste(PQ[x])))+
  xlab("Time of Day")+
  theme(axis.title.x=element_text(size=10,colour = "black"))+
  theme(axis.title.y=element_text(size=10,colour = "black"))+
  theme(axis.text.y=element_text(size=10,colour = "black"))+
  theme(axis.text.x=element_text(size=8,colour = "black"))+
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "4 hours", limits=c(as.POSIXct(data_lowgpp$time[20]), as.POSIXct(data_lowgpp$time[111])))+
  scale_y_continuous(limits=c(-1.3,1.8), breaks=seq(from=-1.2, to=1.8, by=0.6))+
  scale_color_manual(labels=labs, name = "PQ type",values = c("#000000","#D55E00","#0072B2","#009E73"),
                     guide=guide_legend(override.aes=list(linetype=c(1,5,3,4))))+
  scale_linetype_manual(values=c("solid", "dashed","dotted","dotdash"))+
  theme(legend.text=element_text(size=12),legend.position="none")
  
  
#SCENARIO 2: High GPP and low NO3

#PQ_C is the carbon portion that describes the type of organic molecule synthesized. 
#For phytoplakton this value is believed to be 1.1. 
pq_c<-1.1

#PQ_N is the NO3 portion. NO3 is low so this is near zero
pq_n<-0

#PQ_pr is the photorespiration portion. GPP above 100% sat means PR is occurring. 
# We scale the effect by light and a coefficient (0.00061) that represents the literature. In this case a PQ_pr that averages ~ -0.5
pq_pr<- -light*0.00061

#PQ_net is all combined
pq_net<-pq_c+pq_n+pq_pr

#Format for plotting
pq_highgpp_lown<-melt(data.frame(as.POSIXct(dtime,format="%H:%M:%S"),pq_c, pq_n, pq_pr,pq_net ),id="as.POSIXct.dtime..format.....H..M..S..")
names(pq_highgpp_lown)<-c("time","factor","pq_x" )
pq_highgpp_lown$factor<-factor(pq_highgpp_lown$factor, levels= c("pq_net","pq_c", "pq_n", "pq_pr"))


# Plot PQ for SCENARIO 2
p5<-ggplot(pq_highgpp_lown, aes(x=time, y=pq_x, color=factor(factor)))+
  geom_line(aes(linetype=factor),size=1)+
  theme_classic()+
  ylab(expression(paste(PQ[x])))+
  xlab("Time of Day")+
  theme(axis.title.x=element_text(size=10,colour = "black"))+
  theme(axis.title.y=element_text(size=10,colour = "black"))+
  theme(axis.text.y=element_text(size=10,colour = "black"))+
  theme(axis.text.x=element_text(size=8,colour = "black"))+
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "4 hours", limits=c(as.POSIXct(data_lowgpp$time[20]), as.POSIXct(data_lowgpp$time[111])))+
  scale_y_continuous(limits=c(-1.3,1.8), breaks=seq(from=-1.2, to=1.8, by=0.6))+
  scale_color_manual(labels=labs, name = "PQ type",values = c("#000000","#D55E00","#0072B2","#009E73"),
                     guide=guide_legend(override.aes=list(linetype=c(1,5,3,4))))+
  scale_linetype_manual(values=c("solid", "dashed","dotted","dotdash"))+
  theme(legend.position="none")+
  theme(axis.title.y=element_blank())


#SCENARIO 3: High GPP and high NO3

#PQ_C is the carbon portion that describes the type of organic molecule synthesized. 
#For phytoplakton this value is believed to be 1.1. 
pq_c<-1.1

#PQ_N is the NO3 portion. NO3 is high so we scale the effect by light and a coefficient (0.0008) that represents the literature. In this case a PQ_NO3 that averages ~ +0.6
pq_n<- light*0.0008

#PQ_pr is the photorespiration portion. # We scale the effect by light and a coefficient (0.00061) that represents the literature. In this case a PQ_pr that averages ~ -0.5
pq_pr<- -light*0.00061

#PQ_net is all combined
pq_net<-pq_c+pq_n+pq_pr

#Format for plot
pq_highgpp_highn<-melt(data.frame(as.POSIXct(dtime,format="%H:%M:%S"),pq_c, pq_n, pq_pr,pq_net ),id="as.POSIXct.dtime..format.....H..M..S..")
names(pq_highgpp_highn)<-c("time","factor","pq_x" )
labs<-c(expression(PQ['Bio']),
        expression(PQ['C']),
        expression(PQ[NO[3]]),
        expression(PQ['Pr']))
pq_highgpp_highn$factor<-factor(pq_highgpp_highn$factor, levels= c("pq_net","pq_c", "pq_n", "pq_pr"))


# Plot PQ for SCENARIO 3
p6<-ggplot(pq_highgpp_highn, aes(x=time, y=pq_x, color=factor))+
  geom_line(aes(linetype=factor),size=1)+
  theme_classic()+
  ylab(expression(paste(PQ[x])))+
  xlab("Time of Day")+
  theme(axis.title.x=element_text(size=10,colour = "black"))+
  theme(axis.title.y=element_text(size=10,colour = "black"))+
  theme(axis.text.y=element_text(size=10,colour = "black"))+
  theme(axis.text.x=element_text(size=8,colour = "black"))+
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "4 hours", limits=c(as.POSIXct(data_lowgpp$time[20]), as.POSIXct(data_lowgpp$time[111])))+
  scale_y_continuous(limits=c(-1.3,1.8), breaks=seq(from=-1.2, to=1.8, by=0.6))+
  #scale_color_discrete(labels=labs, name = "PQ type")+
  scale_color_manual(labels=labs, name = "PQ type",values = c("#000000","#D55E00","#0072B2","#009E73"),
                     guide=guide_legend(override.aes=list(linetype=c(1,5,3,4))))+
  scale_linetype_manual(values=c("solid", "dashed","dotted","dotdash"))+
  theme(legend.position="right",legend.text=element_text(size=9), legend.title=element_text(size=10), legend.text.align = 0)+
  theme(axis.title.y=element_blank(),legend.key.width = unit(3.5,"line"))+
  guides(linetype = "none")
  


# Calculate carbon-based GPP using simulated oxygen-based GPP --------------


# Common used and literature PQ values
pq_fw1<-1 #Common value used in the freshwater sciences
pq_fw2<-1.2 #Common value used in the freshwater sciences
pq_mar<-1.4 #Common value used in the marine sciences
pq_min<-0.5 #Lower measured value from the literature. This is a conservative value.
pq_max<-3.5 #Higher measured value from the literature. This is a conservative value. 

# SCENARIO 1: Low GPP and low NO3
##Take the mean of PQ_net from Scenario 1
pq_net<-mean(pq_lowgpp[pq_lowgpp$factor == "pq_net",]$pq_x )

##Calculate 
gpp_o2<-19 ##Fixed low oxygen based GPP (in mmmol per meter-cubed per hour)
gpp_c_fw1<-gpp_o2/pq_fw1
gpp_c_fw2<-gpp_o2/pq_fw2
gpp_c_mar<-gpp_o2/pq_mar
gpp_c_min<-gpp_o2/pq_min
gpp_c_max<-gpp_o2/pq_max
gpp_c_true<-gpp_o2/pq_net

#Format
values<-c(gpp_o2, gpp_c_fw1,gpp_c_fw2,gpp_c_mar, gpp_c_true )
var<-c(1,1.5,1.5,1.5,1.5)
factor<-c("", "fw1", "fw2", "mar", "true")
data<-data.frame(values, var, factor)
data$factor<-factor(data$factor, levels= c("true","min", "fw1", "fw2","mar", ""))

# Plot O and C-based metabolism with different PQ values for scenario 1
p7<-ggplot(data, aes(x=var, y=values, group=factor(factor)))+
  geom_jitter(aes(shape=factor(factor)),size=2, width=0.05)+
  theme_classic()+
  ylab(expression(atop("GPP", paste("(mmol ", O[2]," or C ", m^-3," ",hr^-1,")"))))+
  xlab("Metabolism Units")+
  theme(axis.title.x=element_text(size=10,colour = "black"))+
  theme(axis.title.y=element_text(size=10,colour = "black",hjust = 0.5))+
  theme(axis.text.y=element_text(size=10,colour = "black"))+
  theme(axis.text.x=element_text(size=10,colour = "black"))+
  scale_y_continuous(limits=c(0,43), breaks=seq(from=0, to=40, by=10))+
  scale_x_continuous(limits=c(0.8, 1.7), breaks=c(1, 1.5), labels=c(expression(paste(, O[2],)), "C"))+
  scale_shape_manual(values=c(13,15, 16, 17,19),labels=labs, name = "PQ value")+ 
  theme(legend.position="none")+
  theme(axis.title.x=element_blank())+
  annotate("rect", xmin=1.4, xmax=1.6, ymin=gpp_c_min, ymax=gpp_c_max, fill = "black", alpha=.2)


# SCENARIO 2: High GPP and low NO3
pq_net<-mean(pq_highgpp_lown[pq_highgpp_lown$factor == "pq_net",]$pq_x )


gpp_o2<-250 ##Fixed high oxygen based GPP (in mmmol per meter-cubed per hour)
gpp_c_fw1<-gpp_o2/pq_fw1
gpp_c_fw2<-gpp_o2/pq_fw2
gpp_c_mar<-gpp_o2/pq_mar
gpp_c_min<-gpp_o2/pq_min
gpp_c_max<-gpp_o2/pq_max
gpp_c_true<-gpp_o2/pq_net

# format
values<-c(gpp_o2, gpp_c_fw1,gpp_c_fw2,gpp_c_mar, gpp_c_true )
var<-c(1,1.5,1.5,1.5,1.5)
factor<-c("", "fw1", "fw2", "mar", "true")
data<-data.frame(values, var, factor)
data$factor<-factor(data$factor, levels= c("true","min", "fw1", "fw2","mar", ""))

# Plot O and C-based metabolism with different PQ values for scenario 2
p8<-ggplot(data, aes(x=var, y=values, group=factor(factor)))+
  geom_jitter(aes(shape=factor(factor)),size=2, width=0.05)+
  theme_classic()+
  ylab(expression(paste("GPP (mmol ", O[2]," ", m^-3," ",hr^-1 ,")")))+
  xlab("Metabolism Units")+
  theme(axis.title.x=element_text(size=10,colour = "black"))+
  theme(axis.title.y=element_text(size=10,colour = "black"))+
  theme(axis.text.y=element_text(size=10,colour = "black"))+
  theme(axis.text.x=element_text(size=10,colour = "black"))+
  scale_y_continuous(limits=c(50,500), breaks=seq(from=75, to=475, by=100))+
  scale_x_continuous(limits=c(0.8, 1.7), breaks=c(1, 1.5), labels=c(expression(paste(, O[2],)), "C"))+
  scale_shape_manual(values=c(13,15, 16, 17,19),labels=labs, name = "PQ value")+ 
  theme(axis.title.y=element_blank())+
  theme(legend.position="none")+
  annotate("rect", xmin=1.4, xmax=1.6, ymin=gpp_c_min, ymax=gpp_c_max, fill = "black", alpha=.2)


# SCENARIO 2: High GPP and high NO3
pq_net<-mean(pq_highgpp_highn[pq_highgpp_highn$factor == "pq_net",]$pq_x )

gpp_o2<-250 ##Fixed high oxygen based GPP (in mmmol per meter-cubed per hour)
gpp_c_fw1<-gpp_o2/pq_fw1
gpp_c_fw2<-gpp_o2/pq_fw2
gpp_c_mar<-gpp_o2/pq_mar
gpp_c_min<-gpp_o2/pq_min
gpp_c_max<-gpp_o2/pq_max
gpp_c_true<-gpp_o2/pq_net


values<-c(gpp_o2, gpp_c_fw1,gpp_c_fw2,gpp_c_mar, gpp_c_true )
var<-c(1,1.5,1.5,1.5,1.5)
factor<-c("", "fw1", "fw2", "mar", "true")
data<-data.frame(values, var, factor)
data$factor<-factor(data$factor, levels= c("true","min", "fw1", "fw2","mar", ""))
labs<-c(expression(Average~PQ['Bio']),
        "PQ=1 (FW)",
        "PQ=1.2 (FW)",
        "PQ=1.4 (Marine)",
        "NA (true value)")


# Plot O and C-based metabolism with different PQ values for scenario 3
p9<-ggplot(data, aes(x=var, y=values))+
  geom_rect(aes(xmin=1.4, xmax=1.6, ymin=gpp_c_min, ymax=gpp_c_max, fill="Range of GPP\nderived from literature\nPQ range"), color=NA, alpha=0.2)+
  geom_jitter(aes(shape=factor(factor)),size=2, width=0.05)+
  theme_classic()+
  ylab(expression(paste("GPP (mmol ", O[2]," ", m^-3," ",hr^-1 ,")")))+
  xlab("Metabolism Units")+
  theme(axis.title.x=element_text(size=10,colour = "black"))+
  theme(axis.title.y=element_text(size=10,colour = "black"))+
  theme(axis.text.y=element_text(size=10,colour = "black"))+
  theme(axis.text.x=element_text(size=10,colour = "black"))+
  scale_y_continuous(limits=c(50,500), breaks=seq(from=75, to=475, by=100))+
  scale_x_continuous(limits=c(0.8, 1.7), breaks=c(1, 1.5), labels=c(expression(paste(, O[2],)), "C"))+
  scale_shape_manual(values=c(13,15, 16, 17, 19), 
                     breaks= c("true", "fw1", "fw2", "mar", ""), 
                     labels=labs,
                     name = "PQ value")+  
  theme(legend.position="right",legend.text=element_text(size=9), legend.title=element_text(size=10), legend.text.align = 0)+
  theme(axis.title.y=element_blank(),axis.title.x=element_blank())+
  scale_fill_manual('',
                    values = 'gray',  
                    guide = guide_legend(override.aes = list(alpha = 1)))+
  theme(legend.spacing.y=unit(-0.02, "cm"))+
  guides(rect = guide_legend(order=2),
         shape = guide_legend(order=1))

pdf(file="~/GitHub/metabolism-PQ/diel.pdf",
    width=8, height=5)   

#Patch plots together and add labels
patchwork<-p1/p4/p7 | p2/p5/p8 | p3/p6/p9

patchwork + plot_annotation(tag_levels = list(c('A)', 'D)','G)','B)','E)','H)','C)','F)','I)')))
dev.off()
                  
