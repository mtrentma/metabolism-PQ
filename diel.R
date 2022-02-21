library(chron) #this is how we used to do time.  Unwilling to change for this exercise
library(ggplot2)
library(patchwork)
library(reshape2); library(dplyr); library(tidyr)
Sys.setenv(TZ='GMT')
#make up some time data

time<-seq(from = as.numeric(chron(dates="05/21/22", times="02:00:00")),
          to=as.numeric(chron(dates="05/22/22", times="04:00:00")), by=10/1440 )
dtime<- chron(time)


######simulate O2 saturation######

##Make up light data
## now make up light data if you don't have it

## From Yard et al. (1995) Ecological Modelling.  Remember your trig?  
## calculate light as umol photon m-2 s-1.
## Arguments are:  
## time = a date and time input (i.e. a chron object) 
## lat = latitude of field site
## longobs = longitude of field site
## longstd = standard longitude of the field site (NOTE: watch daylight savings time!!!). For PST, longstd is be 120 degrees. But during PDT it is 105 degrees. MST is 105 deg. MDT is 90. 
## year = the year for which you collected data and is entered as "2013-01-01"



# convert degrees to radians
radi<-function(degrees){(degrees*pi/180)}

# function to estimate light
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

## onespecific time just to see if function works
light<- lightest(time=dtime, lat=47,  longobs=114, longstd= 105, year="2022-01-01")

#fake some temp
temp<- rep(15, length(light))  #too cold to swim in

Kcor<-function (temp,K600) {
  K600/(600/(1800.6-(temp*120.1)+(3.7818*temp^2)-(0.047608*temp^3)))^-0.5
}
#end of function

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

#end of function







#Low GPP and low NO3- simulate oxy data
onestationsim<-function(GPP, ER, z, temp, K, light, bp, ts) {
  
  oxy.mod<-numeric(length(light))
  oxy.mod[1]<-9.97  #trial end error this first point based on the O2 data at the end of the time series
  
  # this is the metabolism equation as in Van de Bogert et al (2007) L&OMethods
  for (i in 2:length(light)) { oxy.mod[i]<-oxy.mod[i-1]+((GPP/z)*(light[i]/sum(light)))+ ER*ts/z+(Kcor(temp[i],K))*ts*(osat(temp[i],bp)-oxy.mod[i-1]) }
  oxy.mod
}
# end of function

#tw more functions you need



sim_oxy<- onestationsim(GPP=0.6, ER=-1.5,z=1,temp=temp,light=light, K=15, bp=760, ts=10/1440)
sim_oxy_sat<-100*sim_oxy/osat(temp, 760)
data_lowgpp<-data.frame(dtime, sim_oxy_sat)
data_lowgpp$dtime<-as.POSIXct(data_lowgpp$dtime)
data_lowgpp$time<-as.POSIXct(data_lowgpp$dtime, format="%H:%M:%S")



p1<-
  ggplot(data_lowgpp, aes(x=time, y=sim_oxy_sat))+
  geom_rect(aes( xmin=time[20], xmax=time[111], ymin=-Inf, ymax=Inf, fill="Daytime"), color=NA, alpha=0.2)+
  geom_line(size=1)+
  theme_classic()+
  ylab(expression(paste(O[2], " saturation (%)")))+
  xlab("Time of Day")+
  theme(axis.title.x=element_text(size=14,colour = "black"))+
  theme(axis.title.y=element_text(size=14,colour = "black"))+
  theme(axis.text.y=element_text(size=14,colour = "black"))+
  theme(axis.text.x=element_text(size=14,colour = "black"))+
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "6 hours", limits=c(as.POSIXct("2022-05-21 02:00:00"), as.POSIXct("2022-05-22 02:00:00")))+
  scale_y_continuous(limits=c(90,110), breaks=seq(from=80, to=130, by=10))+
  geom_hline(yintercept=100,linetype="dashed", size=1)+
  #annotate("rect", xmin=data_lowgpp$time[20], xmax=data_lowgpp$time[111], ymin=-Inf, ymax=Inf, fill = "gray", alpha=0.2)+
  ggtitle(expression(paste("GPP= 19 mmol ", O[2]," ", m^-3," ",d^-1 ,"")))+
  theme(axis.title.x=element_blank())+
  theme(plot.title = element_text(size = 12))+
  scale_fill_manual('',
                    values = 'aliceblue',  
                    guide = guide_legend(override.aes = list(alpha = 0.2)))+
  theme(legend.text=element_text(size=12),legend.position="none")
  
  

#High GPP- simulate oxy data
onestationsim<-function(GPP, ER, z, temp, K, light, bp, ts) {
  
  oxy.mod<-numeric(length(light))
  oxy.mod[1]<-9.51  #trial end error this first point based on the O2 data at the end of the time series
  
  # this is the metabolism equation as in Van de Bogert et al (2007) L&OMethods
  for (i in 2:length(light)) { oxy.mod[i]<-oxy.mod[i-1]+((GPP/z)*(light[i]/sum(light)))+ ER*ts/z+(Kcor(temp[i],K))*ts*(osat(temp[i],bp)-oxy.mod[i-1]) }
  oxy.mod
}
# end of function

#tw more functions you need



sim_oxy<- onestationsim(GPP=8, ER=-8,z=1,temp=temp,light=light, K=15, bp=760, ts=10/1440)
sim_oxy_sat<-100*sim_oxy/osat(temp, 760)
data_highgpp<-data.frame(dtime, sim_oxy_sat)
data_highgpp$dtime<-as.POSIXct(data_lowgpp$dtime)
data_highgpp$time<-as.POSIXct(data_lowgpp$dtime, format="%H:%M:%S")



p2<-ggplot(data_highgpp, aes(x=time, y=sim_oxy_sat))+
  geom_rect(aes( xmin=time[20], xmax=time[111], ymin=-Inf, ymax=Inf, fill="Daytime"), color=NA, alpha=0.2)+
  geom_line(size=1)+
  theme_classic()+
  ylab(expression(paste(O[2], " saturation (%)")))+
  xlab("Time of Day")+
  theme(axis.title.x=element_text(size=14,colour = "black"))+
  theme(axis.title.y=element_text(size=14,colour = "black"))+
  theme(axis.text.y=element_text(size=14,colour = "black"))+
  theme(axis.text.x=element_text(size=14,colour = "black"))+
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "6 hours", limits=c(as.POSIXct("2022-05-21 02:00:00"), as.POSIXct("2022-05-22 02:00:00")))+
  scale_y_continuous(limits=c(90,110), breaks=seq(from=80, to=130, by=10))+
  geom_hline(yintercept=100,linetype="dashed", size=1)+
  ggtitle(expression(paste("GPP= 250 mmol ", O[2]," ", m^-3," ",d^-1 ,"; low ", NO[3])))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(plot.title = element_text(size = 12))+
  scale_fill_manual('',
                    values = 'aliceblue',  
                    guide = guide_legend(override.aes = list(alpha = 0.2)))+
  theme(legend.text=element_text(size=12),legend.position="none")

p3<-ggplot(data_highgpp, aes(x=time, y=sim_oxy_sat))+
  geom_rect(aes( xmin=time[20], xmax=time[111], ymin=-Inf, ymax=Inf, fill="Daytime"), color=NA, alpha=0.2)+
  geom_line(size=1)+
  theme_classic()+
  ylab(expression(paste(O[2], " saturation (%)")))+
  xlab("Time of Day")+
  theme(axis.title.x=element_text(size=14,colour = "black"))+
  theme(axis.title.y=element_text(size=14,colour = "black"))+
  theme(axis.text.y=element_text(size=14,colour = "black"))+
  theme(axis.text.x=element_text(size=14,colour = "black"))+
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "6 hours", limits=c(as.POSIXct("2022-05-21 02:00:00"), as.POSIXct("2022-05-22 02:00:00")))+
  scale_y_continuous(limits=c(90,110), breaks=seq(from=80, to=130, by=10))+
  geom_hline(yintercept=100,linetype="dashed", size=1)+
  ggtitle(expression(paste("GPP= 250 mmol ", O[2]," ", m^-3," ",d^-1 ,"; high ", NO[3])))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(plot.title = element_text(size = 12))+
  scale_fill_manual('',
                    values = 'aliceblue',  
                    guide = guide_legend(override.aes = list(alpha = 1)))+
  theme(legend.position="right",legend.text=element_text(size=12), legend.title=element_text(size=15))


######simulate PQ data######

##Low GPP and low NO3

#PQ_C is the carbon portion that describes the type of organic molecule synthesized. 
#For phytoplakton this value is believed to be 1.1. Although, lipid production will likely increase this value
#during the photic period. 
pq_c_start<-rep(1.1,times=length(dtime))
pq_c<-pq_c_start+light*0.000001

#PQ_N is the NO3 portion. NO3 is low so this is near zero
pq_n_start<-rep(0.01,times=length(dtime))
pq_n<-pq_n_start+light*0.000001

#PQ_pr is the photorespiration portion. Low GPP below 100% sat means no PR
pq_pr<-0.001

#PQ_net is all combined
pq_net<-pq_c+pq_n+pq_pr

pq_lowgpp<-melt(data.frame(as.POSIXct(dtime,format="%H:%M:%S"),pq_c, pq_n, pq_pr,pq_net ),id="as.POSIXct.dtime..format.....H..M..S..")
names(pq_lowgpp)<-c("time","factor","pq_x" )
labs<-c(expression(phantom("--")*PQ['net']),
        expression(PQ['c']*phantom("-")),
        expression(PQ['N']*phantom("-")),
        expression(PQ['Pr']))
pq_lowgpp$factor<-factor(pq_lowgpp$factor, levels= c("pq_net","pq_c", "pq_n", "pq_pr"))

p4<-ggplot(pq_lowgpp, aes(x=time, y=pq_x, color=factor(factor)))+
  annotate("rect", xmin=data_lowgpp$time[20], xmax=data_lowgpp$time[111], ymin=-Inf, ymax=Inf, fill = "aliceblue", alpha=1)+
  geom_line(aes(linetype=factor),size=1.5)+
  theme_classic()+
  ylab(expression(paste(PQ[x])))+
  xlab("Time of Day")+
  theme(axis.title.x=element_text(size=14,colour = "black"))+
  theme(axis.title.y=element_text(size=14,colour = "black"))+
  theme(axis.text.y=element_text(size=14,colour = "black"))+
  theme(axis.text.x=element_text(size=14,colour = "black"))+
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "4 hours", limits=c(as.POSIXct(data_lowgpp$time[20]), as.POSIXct(data_lowgpp$time[111])))+
  scale_y_continuous(limits=c(-1.3,1.8), breaks=seq(from=-1.2, to=1.8, by=0.3))+
  scale_color_manual(labels=labs, name = "PQ type",values = c("#000000","#D55E00","#0072B2","#009E73"))+
  scale_linetype_manual(values=c("solid", "dashed","dashed","dashed"))+
  theme(legend.text=element_text(size=12),legend.position="none")
  
  
##High GPP and low NO3

#PQ_C is the carbon portion that describes the type of organic molecule synthesized. 
#For phytoplakton this value is believed to be 1.1. Although, lipid production will likely increase this value
#during the photic period. 
pq_c_start<-rep(1.1,times=length(dtime))
pq_c<-pq_c_start+light*0.000055

#PQ_N is the NO3 portion. NO3 is low so this is near zero
pq_n<-light*0.00001

#PQ_pr is the photorespiration portion. Low GPP below 100% sat means no PR
pq_pr<- -light*0.00061

#PQ_net is all combined
pq_net<-pq_c+pq_n+pq_pr

pq_highgpp_lown<-melt(data.frame(as.POSIXct(dtime,format="%H:%M:%S"),pq_c, pq_n, pq_pr,pq_net ),id="as.POSIXct.dtime..format.....H..M..S..")
names(pq_highgpp_lown)<-c("time","factor","pq_x" )

pq_highgpp_lown$factor<-factor(pq_highgpp_lown$factor, levels= c("pq_net","pq_c", "pq_n", "pq_pr"))



p5<-ggplot(pq_highgpp_lown, aes(x=time, y=pq_x, color=factor(factor)))+
  annotate("rect", xmin=data_lowgpp$time[20], xmax=data_lowgpp$time[111], ymin=-Inf, ymax=Inf, fill = "aliceblue", alpha=1)+
  geom_line(aes(linetype=factor),size=1.5)+
  theme_classic()+
  ylab(expression(paste(PQ[x])))+
  xlab("Time of Day")+
  theme(axis.title.x=element_text(size=14,colour = "black"))+
  theme(axis.title.y=element_text(size=14,colour = "black"))+
  theme(axis.text.y=element_text(size=14,colour = "black"))+
  theme(axis.text.x=element_text(size=14,colour = "black"))+
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "4 hours", limits=c(as.POSIXct(data_lowgpp$time[20]), as.POSIXct(data_lowgpp$time[111])))+
  scale_y_continuous(limits=c(-1.3,1.8), breaks=seq(from=-1.2, to=1.8, by=0.3))+
  scale_color_manual(labels=labs, name = "PQ type",values = c("#000000","#D55E00","#0072B2","#009E73"))+
  scale_linetype_manual(values=c("solid", "dashed","dashed","dashed"))+
  theme(legend.position="none")+
  theme(axis.title.y=element_blank())


##High GPP and low NO3

#PQ_C is the carbon portion that describes the type of organic molecule synthesized. 
#For phytoplakton this value is believed to be 1.1. Although, lipid production will likely increase this value
#during the photic period. 
pq_c_start<-rep(1.1,times=length(dtime))
pq_c<-pq_c_start+light*0.0001

#PQ_N is the NO3 portion. NO3 is low so this is near zero
pq_n<- light*0.0008

#PQ_pr is the photorespiration portion. Low GPP below 100% sat means no PR
pq_pr<- -light*0.00061

#PQ_net is all combined
pq_net<-pq_c+pq_n+pq_pr

pq_highgpp_highn<-melt(data.frame(as.POSIXct(dtime,format="%H:%M:%S"),pq_c, pq_n, pq_pr,pq_net ),id="as.POSIXct.dtime..format.....H..M..S..")
names(pq_highgpp_highn)<-c("time","factor","pq_x" )
labs<-c(expression(phantom("--")*PQ['net']),
        expression(PQ['c']*phantom("-")),
        expression(PQ['N']*phantom("-")),
        expression(PQ['Pr']))

pq_highgpp_highn$factor<-factor(pq_highgpp_highn$factor, levels= c("pq_net","pq_c", "pq_n", "pq_pr"))


p6<-ggplot(pq_highgpp_highn, aes(x=time, y=pq_x, color=factor))+
  annotate("rect", xmin=data_lowgpp$time[20], xmax=data_lowgpp$time[111], ymin=-Inf, ymax=Inf, fill = "aliceblue", alpha=1)+
  geom_line(aes(linetype=factor),size=1.5)+
  theme_classic()+
  ylab(expression(paste(PQ[x])))+
  xlab("Time of Day")+
  theme(axis.title.x=element_text(size=14,colour = "black"))+
  theme(axis.title.y=element_text(size=14,colour = "black"))+
  theme(axis.text.y=element_text(size=14,colour = "black"))+
  theme(axis.text.x=element_text(size=14,colour = "black"))+
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "4 hours", limits=c(as.POSIXct(data_lowgpp$time[20]), as.POSIXct(data_lowgpp$time[111])))+
  scale_y_continuous(limits=c(-1.3,1.8), breaks=seq(from=-1.2, to=1.8, by=0.3))+
  #scale_color_discrete(labels=labs, name = "PQ type")+
  scale_color_manual(labels=labs, name = "PQ type",values = c("#000000","#D55E00","#0072B2","#009E73"))+
  scale_linetype_manual(values=c("solid", "dashed","dashed","dashed"))+
  theme(legend.position="right",legend.text=element_text(size=12),legend.title=element_text(size=15))+
  theme(axis.title.y=element_blank())+
  guides(linetype = "none")   


##GPP estimates with different PQs
pq_fw1<-1
pq_fw2<-1.2
pq_mar<-1.4
pq_min<-0.5
pq_max<-3.5

#low GPP and N
pq_net<-mean(pq_lowgpp[pq_lowgpp$factor == "pq_net",]$pq_x )

gpp_o2<-20
gpp_c_fw1<-gpp_o2/pq_fw1
gpp_c_fw2<-gpp_o2/pq_fw2
gpp_c_mar<-gpp_o2/pq_mar
gpp_c_min<-gpp_o2/pq_min
gpp_c_max<-gpp_o2/pq_max
gpp_c_true<-gpp_o2/pq_net


values<-c(gpp_o2, gpp_c_fw1,gpp_c_fw2,gpp_c_mar, gpp_c_true )
#values<-c(gpp_o2, gpp_c_fw1,gpp_c_fw2,gpp_c_mar,gpp_c_min, gpp_c_max, gpp_c_true )
#var<-c("O2", "C", "C", "C", "C", "C", "C")
var<-c(1,1.5,1.5,1.5,1.5)
#factor<-c("", "fw1", "fw2", "mar","min", "max", "true")
factor<-c("", "fw1", "fw2", "mar", "true")
data<-data.frame(values, var, factor)
#data$var<-factor(data$var, levels= c("O2", "C"))
#data$factor<-factor(data$factor, levels= c("true","min", "fw1", "fw2", "mar", "max", ""))
data$factor<-factor(data$factor, levels= c("true","min", "fw1", "fw2","mar", ""))

p7<-ggplot(data, aes(x=var, y=values, group=factor(factor)))+
  geom_jitter(aes(shape=factor(factor)),size=4, width=0.05)+
  theme_classic()+
  ylab(expression(paste("GPP (mmol ", O[2]," or C ", m^-3," ",hr^-1 ,")")))+
  xlab("Metabolism Units")+
  theme(axis.title.x=element_text(size=14,colour = "black"))+
  theme(axis.title.y=element_text(size=14,colour = "black"))+
  theme(axis.text.y=element_text(size=14,colour = "black"))+
  theme(axis.text.x=element_text(size=14,colour = "black"))+
  scale_y_continuous(limits=c(0,45), breaks=seq(from=0, to=45, by=5))+
  scale_x_continuous(limits=c(0.8, 1.7), breaks=c(1, 1.5), labels=c(expression(paste(, O[2],)), "C"))+
  #scale_x_discrete(labels=c(expression(paste(, O[2],)), "C"))+
  #scale_shape_manual(values=c(13, 0,15, 16, 17, 2,19),labels=labs, name = "PQ value")+
  scale_shape_manual(values=c(13,15, 16, 17,19),labels=labs, name = "PQ value")+ 
  theme(legend.position="none")+
  theme(axis.title.x=element_blank())+
  annotate("rect", xmin=1.4, xmax=1.6, ymin=gpp_c_min, ymax=gpp_c_max, fill = "black", alpha=.2)


#high GPP and low N
pq_net<-mean(pq_highgpp_lown[pq_highgpp_lown$factor == "pq_net",]$pq_x )


gpp_o2<-250
gpp_c_fw1<-gpp_o2/pq_fw1
gpp_c_fw2<-gpp_o2/pq_fw2
gpp_c_mar<-gpp_o2/pq_mar
gpp_c_min<-gpp_o2/pq_min
gpp_c_max<-gpp_o2/pq_max
gpp_c_true<-gpp_o2/pq_net

values<-c(gpp_o2, gpp_c_fw1,gpp_c_fw2,gpp_c_mar, gpp_c_true )
#values<-c(gpp_o2, gpp_c_fw1,gpp_c_fw2,gpp_c_mar,gpp_c_min, gpp_c_max, gpp_c_true )
#var<-c("O2", "C", "C", "C", "C", "C", "C")
var<-c(1,1.5,1.5,1.5,1.5)
factor<-c("", "fw1", "fw2", "mar", "true")
#factor<-c("", "fw1", "fw2", "mar","min", "max", "true")
data<-data.frame(values, var, factor)
#data$var<-factor(data$var, levels= c("O2", "C"))
#data$factor<-factor(data$factor, levels= c("true","min", "fw1", "fw2", "mar", "max", ""))
data$factor<-factor(data$factor, levels= c("true","min", "fw1", "fw2","mar", ""))

p8<-ggplot(data, aes(x=var, y=values, group=factor(factor)))+
  geom_jitter(aes(shape=factor(factor)),size=4, width=0.05)+
  theme_classic()+
  ylab(expression(paste("GPP (mmol ", O[2]," ", m^-3," ",hr^-1 ,")")))+
  xlab("Metabolism Units")+
  theme(axis.title.x=element_text(size=14,colour = "black"))+
  theme(axis.title.y=element_text(size=14,colour = "black"))+
  theme(axis.text.y=element_text(size=14,colour = "black"))+
  theme(axis.text.x=element_text(size=14,colour = "black"))+
  scale_y_continuous(limits=c(70,520), breaks=seq(from=70, to=520, by=50))+
  scale_x_continuous(limits=c(0.8, 1.7), breaks=c(1, 1.5), labels=c(expression(paste(, O[2],)), "C"))+
  #scale_x_discrete(labels=c(expression(paste(, O[2],)), "C"))+
  #scale_shape_manual(values=c(13, 0,15, 16, 17, 2,19),labels=labs, name = "PQ value")+  theme(legend.position="none")+
  scale_shape_manual(values=c(13,15, 16, 17,19),labels=labs, name = "PQ value")+ 
  theme(axis.title.y=element_blank())+
  theme(legend.position="none")+
  annotate("rect", xmin=1.4, xmax=1.6, ymin=gpp_c_min, ymax=gpp_c_max, fill = "black", alpha=.2)


#high GPP and low N
pq_net<-mean(pq_highgpp_highn[pq_highgpp_highn$factor == "pq_net",]$pq_x )

gpp_o2<-250
gpp_c_fw1<-gpp_o2/pq_fw1
gpp_c_fw2<-gpp_o2/pq_fw2
gpp_c_mar<-gpp_o2/pq_mar
gpp_c_min<-gpp_o2/pq_min
gpp_c_max<-gpp_o2/pq_max
gpp_c_true<-gpp_o2/pq_net


values<-c(gpp_o2, gpp_c_fw1,gpp_c_fw2,gpp_c_mar, gpp_c_true )
#var<-c("O2", "C", "C", "C", "C", "C", "C")
var<-c(1,1.5,1.5,1.5,1.5)
factor<-c("", "fw1", "fw2", "mar", "true")
data<-data.frame(values, var, factor)
#data$var<-factor(data$var, levels= c("O2", "C"))
data$factor<-factor(data$factor, levels= c("true","min", "fw1", "fw2","mar", ""))

labs<-c("Net PQ"," PQ=1 (FW)", "PQ=1.2 (FW)", "PQ=1.4 (Marine)","")


p9<-ggplot(data, aes(x=var, y=values))+
  geom_rect(aes(xmin=1.4, xmax=1.6, ymin=gpp_c_min, ymax=gpp_c_max, fill="GPP of Literature \nPQ Range"), color=NA, alpha=0.2)+
  geom_jitter(aes(shape=factor(factor)),size=4, width=0.05)+
  theme_classic()+
  ylab(expression(paste("GPP (mmol ", O[2]," ", m^-3," ",hr^-1 ,")")))+
  xlab("Metabolism Units")+
  theme(axis.title.x=element_text(size=14,colour = "black"))+
  theme(axis.title.y=element_text(size=14,colour = "black"))+
  theme(axis.text.y=element_text(size=14,colour = "black"))+
  theme(axis.text.x=element_text(size=14,colour = "black"))+
  scale_y_continuous(limits=c(70,520), breaks=seq(from=70, to=520, by=50))+
  scale_x_continuous(limits=c(0.8, 1.7), breaks=c(1, 1.5), labels=c(expression(paste(, O[2],)), "C"))+
  #scale_x_discrete(labels=c(expression(paste(, O[2],)), "C"))+
  scale_shape_manual(values=c(13,15, 16, 17,19),labels=labs, name = "PQ value")+  
  theme(legend.position="right",legend.text=element_text(size=12), legend.title=element_text(size=15))+
  theme(axis.title.y=element_blank(),axis.title.x=element_blank())+
   scale_fill_manual('',
                    values = 'gray',  
                    guide = guide_legend(override.aes = list(alpha = 1)))+
  theme(legend.spacing.y=unit(-0.03, "cm"))




p1/p4/p7 | p2/p5/p8 | p3/p6/p9

