
# Beginning Figure 4---------------------------------------------------------------


#load all packages
library(plyr); library(dplyr)
library(tidyr)
library(ggplot2)
library(zoo)
library(streamMetabolizer)
library(lubridate)
library(rnoaa)
library(rgdal)
library(sp)
library(httr)
library(Rcpp)
library(GSODR)
library(dataRetrieval)
library(shinystan)
library(tidyverse)
library(patchwork)
library(reshape2); library(dplyr); library(tidyr)
options(noaakey = "RocTXcnTfTiprFMtvhljUQnloOuZtpeO")
#To access mesowest data
library(mesowest)
library(jsonlite)
requestToken(apikey = "KyW2GDQXm7iuxMZIttsuCK15e1zWMzvOTrUpX9k3HN")

##Save MSO airport ID for downloading pressure
station<-'727730-24153' #MSO airport

##Location for estimating light. Same for all sites. 
lat<-46.790812
long<--113.748535

##Elevation of sites (masl)
ele.BN<-1103
ele.PL<-1451

##Download pressure and air temp data from DeerLodge
res<-GET("https://api.synopticdata.com/v2/stations/timeseries?stid=k38s&start=202001030000&end=202012200000&timeformat=%s&token=bc312c7c369043269a005c049228c49b")

#reorganize pressure data
data = fromJSON(rawToChar(res$content))
press.dl<-as.data.frame(data[["STATION"]][["OBSERVATIONS"]])
press.dl.unlist<-unlist(press.dl$date_time)
press.dl.unlist<-as.data.frame(press.dl.unlist)
press.dl.unlist$pressure<-unlist(press.dl$pressure_set_1d)
press.dl.unlist$press.dl.unlist<-as.POSIXct(as.numeric(press.dl.unlist$press.dl.unlist),origin='1970-01-01')
press.dl.unlist$date<-as.Date(press.dl.unlist$press.dl.unlist)

pressure_day<-press.dl.unlist %>%
  group_by(date) %>% 
  summarise(mean.press=mean(pressure/100))

names(pressure_day)<-c("date", "press.mb")

#reorganize air temp data
data = fromJSON(rawToChar(res$content))
temp.dl<-as.data.frame(data[["STATION"]][["OBSERVATIONS"]])
temp.dl.unlist<-unlist(temp.dl$date_time)
temp.dl.unlist<-as.data.frame(temp.dl.unlist)
temp.dl.unlist$air_temp<-unlist(temp.dl$air_temp_set_1)
temp.dl.unlist$temp.dl.unlist<-as.POSIXct(as.numeric(temp.dl.unlist$temp.dl.unlist),origin='1970-01-01')
temp.dl.unlist$date<-as.Date(temp.dl.unlist$temp.dl.unlist)

air_temp_day<-temp.dl.unlist %>%
  group_by(date) %>% 
  summarise(mean.temp=mean(air_temp))

names(air_temp_day)<-c("date", "mean.temp")

##set working directory
setwd("~/GitHub/metabolism-PQ/Figure 4")


# Panel A (Site 1 % sat)--------------------------------------------------------------

dat <- read_csv("JAP_DO_MINIDOT_PERKINS_2020.csv")
names(dat)<-c("unix.time", "date.UTC", "date.MST", "battery", "temp.c", "do.mgl", "do.sat","q")

dat$date<-as.Date(dat$date.UTC,format="%m-%d-%Y")


# convert to solar time
time<-as.POSIXct(paste(dat$date.UTC),format="%Y-%m-%d %H:%M:%S",tz="UTC")
dat$solar.time <- convert_UTC_to_solartime(
  time,
  long,
  time.type = c("mean solar")
)

#Alt: as.POSIXct(paste(dat$date.UTC),format="%Y-%m-%d %H:%M:%S",tz="UTC")
#lubridate::as_datetime(dat$unix.time)

# clean up
metab <- data.frame(dat$solar.time,dat$do.mgl,dat$temp.c)
metab$date <- as.Date(as.POSIXct(dat$solar.time))
colnames(metab) <- c("solar.time","DO.obs","temp.water", "date")

# Generate light
metab$light<-calc_light(metab$solar.time, latitude = lat, longitude= long)


# Get station (not sea level) pressure data from Missoula airport
#GSOD<-get_GSOD(years=2020, station=station)
#pressure<-data.frame(GSOD$YEARMODA, GSOD$STP)
#names(pressure)<-c("date", "press.mb")


##Pressure by day
pressure<-data.frame(pressure_day$press.mb, air_temp_day$mean.temp, air_temp_day$date)
colnames(pressure)<- c('press.mb','mean.temp', "date")


metab<-left_join(metab,pressure)

# Calculate DO concentration at saturation
metab$DO.sat <- calc_DO_sat(
  temp=metab$temp.water,
  press=metab$press.mb,
  salinity.water = 0,
)

# Calculate O2 sat
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

perkins<-NA
metab$DO.sat.perc<-100*metab$DO.obs/osat(metab$temp.water, (metab$press.mb/1.33))
perkins<-metab %>%
  group_by(date) %>%
  summarise(max = max(DO.sat.perc))


p1<-ggplot(data=perkins, aes(x=as.POSIXct(date), y=max))+
  geom_point(color="black", size=2)+
  ylab(expression(atop("Maximum",paste(O[2], " saturation (%)"))))+
  xlab("Date")+
  #geom_point(aes(y=y, x=time),color="red", size=3)+
  theme_classic()+
  theme(axis.title.x=element_text(size=16,colour = "black"))+
  theme(axis.title.y=element_text(size=16,colour = "black"))+
  theme(axis.text.y=element_text(size=14,colour = "black"))+
  theme(axis.text.x=element_text(size=14,colour = "black"))+
  scale_y_continuous(limits=c(85,130), breaks=seq(from=90, to=130, by=10))+
  geom_hline(yintercept = 100)+
  ggtitle("Site 1")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+
  scale_x_datetime(date_labels = "%b", date_breaks = "1 month", limits=c(as.POSIXct("2020-07-25"), as.POSIXct("2020-11-07")))


# Panel B (Site 2 % sat) --------------------------------------------------------------

dat <- read_csv("JAP_DO_MINIDOT_BONITA_2020.csv")
names(dat)<-c("unix.time", "date.UTC", "date.MST", "battery", "temp.c", "do.mgl", "do.sat","q")

dat$date<-as.Date(dat$date.UTC,format="%m-%d-%Y")


# convert to solar time
time<-as.POSIXct(paste(dat$date.UTC),format="%Y-%m-%d %H:%M:%S",tz="UTC")
dat$solar.time <- convert_UTC_to_solartime(
  time,
  long,
  time.type = c("mean solar")
)

#Alt: as.POSIXct(paste(dat$date.UTC),format="%Y-%m-%d %H:%M:%S",tz="UTC")
#lubridate::as_datetime(dat$unix.time)

# clean up
metab <- data.frame(dat$solar.time,dat$do.mgl,dat$temp.c)
metab$date <- as.Date(as.POSIXct(dat$solar.time))
colnames(metab) <- c("solar.time","DO.obs","temp.water", "date")

# Generate light
metab$light<-calc_light(metab$solar.time, latitude = lat, longitude= long)

# Add pressure
metab<-left_join(metab,pressure)

# Calculate DO concentration at saturation
metab$DO.sat <- calc_DO_sat(
  temp=metab$temp.water,
  press=metab$press.mb,
  salinity.water = 0,
)

# Calculate O2 sat
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

metab$DO.sat.perc<-100*metab$DO.obs/osat(metab$temp.water, (metab$press.mb/1.33))
bonita<-metab %>%
  group_by(date) %>%
  summarise(max = max(DO.sat.perc))


p2<-ggplot(data=bonita, aes(x=as.POSIXct(date), y=max))+
  geom_point(color="black", size=2)+
  ylab(expression(paste(O[2], " saturation (%)")))+
  xlab("Date")+
  #geom_point(aes(y=y, x=time),color="red", size=3)+
  theme_classic()+
  theme(axis.title.x=element_text(size=16,colour = "black"))+
  theme(axis.title.y=element_text(size=16,colour = "black"))+
  theme(axis.text.y=element_text(size=14,colour = "black"))+
  theme(axis.text.x=element_text(size=14,colour = "black"))+
  scale_y_continuous(limits=c(85,130), breaks=seq(from=90, to=130, by=10))+
  geom_hline(yintercept = 100)+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
              axis.text.x=element_blank())+
  ggtitle("Site 2")+
  scale_x_datetime(date_labels = "%b", date_breaks = "1 month", limits=c(as.POSIXct("2020-07-25"), as.POSIXct("2020-11-07")))



# Panels C and D (nutrient data) ------------------------------------------
setwd("~/GitHub/metabolism-PQ/Figure 4/nutrient data")
pic.data <-
  list.files(pattern="\\.csv$") %>% 
  map_df(~read_csv(.))

names(pic.data)<-c("row", "project", "date", "day", "month", "year", "site", "rep", "nh4.mgl", "srp.mgl", "no3.mgl")
pic.data$date<-as.POSIXct(pic.data$date, format="%m/%d/%Y")

NO3<-subset(pic.data, no3.mgl>0)
NO3<-NO3[,c(1:8,11)]
NH4<-pic.data[is.na(pic.data$no3.mgl),]
NH4<-NH4[,1:10]

data<-full_join(NH4, NO3)
data$site<-as.factor(data$site)
data$rep<-as.factor(data$rep)

data.mean<-ddply(data, c('date','site'), summarize,
                 mean.nh4=mean(nh4.mgl, na.rm=TRUE),
                 mean.no3=mean(no3.mgl, na.rm=TRUE),
                 mean.srp=mean(srp.mgl, na.rm=TRUE))

data.mean.sub<-subset(data.mean, site==2|site==11)
data.mean.sub$site<-factor(data.mean.sub$site,levels=c(2,11))
data.mean.sub$no3.perc<-data.mean.sub$mean.no3/(data.mean.sub$mean.nh4+data.mean.sub$mean.no3)
data.mean.sub.sum<-subset(data.mean.sub, date> as.POSIXct("2020-07-20"))

nutr<-subset(data.mean.sub.sum, site==2)
out<-seq(from=as.Date("2020-08-07"), to=as.Date("2020-10-27"), by=1)
perkins$perc.DIN<-spline(as.Date(nutr$date), nutr$no3.perc, xout=out, method="natural")$y


##Plot nutrient data for Site 1

p3<-ggplot(data=perkins,aes(x=as.POSIXct(date), y=perc.DIN))+
  geom_line(size=1.5)+
  theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=16,color="black"))+
  theme(axis.title.y=element_text(size=16,color="black",hjust = 0.5))+
  theme(axis.text.y=element_text(size=14,color="black"))+
  theme(axis.text.x=element_text(size=14,color="black"))+
  xlab("Date")+
  ylab(expression(atop("Proportion of DIN", paste("present as ",NO[3],))))+
  scale_y_continuous(limits=c(0,1), breaks=seq(from=0, to=1, by=.20))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+
  scale_x_datetime(date_labels = "%b", date_breaks = "1 month", limits=c(as.POSIXct("2020-07-25"), as.POSIXct("2020-11-07")))

  


nutr<-subset(data.mean.sub.sum, site==11)
out<-seq(from=as.Date("2020-08-07"), to=as.Date("2020-10-27"), by=1)

bonita$perc.DIN<-spline(as.Date(nutr$date), nutr$no3.perc, xout=out, method="natural")$y
bonita$perc.DIN[67:82]<-seq(from=.94,to=.91, length.out=16)

##Plot nutrient data for Site 2


p4<-ggplot(data=bonita,aes(x=as.POSIXct(date), y=perc.DIN))+
  geom_line(size=1.5)+
  theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=16,color="black"))+
  theme(axis.title.y=element_text(size=16,color="black"))+
  theme(axis.text.y=element_text(size=14,color="black"))+
  theme(axis.text.x=element_text(size=14,color="black"))+
  xlab("Date")+
  ylab(expression(paste("Prop. ",NO[3]," of DIN")))+
  scale_y_continuous(limits=c(0,1), breaks=seq(from=0, to=1, by=.20))+
  theme(axis.title.y=element_blank(),
       axis.text.y=element_blank(),
      axis.title.x=element_blank(),
     axis.text.x=element_blank())+
  scale_x_datetime(date_labels = "%b", date_breaks = "1 month", limits=c(as.POSIXct("2020-07-25"), as.POSIXct("2020-11-07")))

  

# Panels E and F (predicted PQ_bio) ---------------------------------------


####Estimating the effects of photorespiration on the PQ.
###Data plucked from Burris 1981
setwd("C:/Users/Matt Trentman/IDrive-Sync/Postdoc/Chamber metabolism")
temp<-20#C
pres<-1#atm
data<-read.csv("burris_pr_change.csv")
data$PQ.delta<-data$PQ-1.2
data$O2.pp.kpa<-data$O2.pp*101
data$O2.mgL<-data$O2.pp.kpa*0.36610

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

data$DO.sat.perc<-100*data$O2.mgL/osat(temp, (pres*760))

plot(data$PQ~data$DO.sat.perc, data=data)

model<-lm(data$PQ~data$DO.sat.perc)
summary(model)



##Create NO3 proportion of DIN. Proportion from 0 to 1
DIN<-seq(from=0, to=1, length.out=400)

#Calculate PQ change due to NO3
alg_CN<-5 #must define algal C:N ratio  
y<-DIN*2/alg_CN


##Create O2 and PQ gradient (based on Burris 1981 recalculated to percent sat)
O2<-seq(from=40, to=180, length.out=400)
PQ<-seq(from=0.56, to=-1.2, length.out=400)


##Expand to all combinations
act.data<-expand.grid(X=O2,Y= DIN)
z<-expand.grid(PQ, y)

#add two individual effects together
act.data$Z<-z$Var1+z$Var2

predict.bon<-NA
predict.bon$max_O2<-round(act.data$X)
predict.bon$perc_DIN<-act.data$Y
predict.bon$prediction<-act.data$Z

model<-lm(prediction~max_O2+perc_DIN, data=predict.bon)

# Predictions-PQc is a fixed value of 1.1 plus predictions for PQ_pr and PQ_no3
bonita$predictions<-1.1+((model$coefficients[2]*bonita$max)+(model$coefficients[2]*bonita$perc.DIN)+model$coefficients[1])
perkins$predictions<-1.1+((model$coefficients[2]*perkins$max)+(model$coefficients[2]*perkins$perc.DIN)+model$coefficients[1])



p5<-ggplot(data=perkins,aes(x=as.POSIXct(date), y=predictions))+
  geom_line(size=1.1)+
  theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=16,color="black"))+
  theme(axis.title.y=element_text(size=16,color="black"))+
  theme(axis.text.y=element_text(size=14,color="black"))+
  theme(axis.text.x=element_text(size=14,color="black"))+
  xlab("Date")+
  ylab(expression(paste("Predicted ",PQ[Bio],)))+
  scale_y_continuous(limits=c(0.4,1.), breaks=seq(from=0.4, to=1.0, by=.10))+
  scale_x_datetime(date_labels = "%b", date_breaks = "1 month", limits=c(as.POSIXct("2020-07-25"), as.POSIXct("2020-11-07")))

p6<-ggplot(data=bonita,aes(x=as.POSIXct(date), y=predictions))+
  geom_line(size=1.1)+
  theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=16,color="black"))+
  theme(axis.title.y=element_text(size=16,color="black"))+
  theme(axis.text.y=element_text(size=14,color="black"))+
  theme(axis.text.x=element_text(size=14,color="black"))+
  xlab("Date")+
  ylab(expression(paste("Predicted ",PQ[Bio],)))+
  scale_y_continuous(limits=c(0.4,1.), breaks=seq(from=0.4, to=1.0, by=.10))+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank())+
  scale_x_datetime(date_labels = "%b", date_breaks = "1 month", limits=c(as.POSIXct("2020-07-25"), as.POSIXct("2020-11-07")))


patchwork<-p1/p3/p5|p2/p4/p6

patchwork+plot_annotation(tag_levels = list(c('A)', 'C)','E)','B)','D)', 'F)')))


# Beginning Figure 5---------------------------------------------------------------



# Panel A (site 1 O2-based metabolism) ------------------------------------

##Load metabolism data
setwd("~/GitHub/metabolism-PQ/Figure 4")
mydata = read.csv("JAP_METABOLISM_BONITA_2020.csv",header=T)
mydata<-mydata[,c(1:11)]
names(mydata)<-c("date", "discharge.cms", "depth.m", "width.m", "GPP", "ER","GPP.low","GPP.high","ER.low", "ER.high", "site")
mydata$date<-as.POSIXct(mydata$date, format="%m/%d/%Y")
mydata$gpp_mol_O2<-mydata$GPP/32
mydata$gpp_mol_O2[81]<-0.024
tem<-mydata$gpp_mol_O2
tem[82]<-0.023
bonita$gpp_mol_O2<-tem
bonita$gpp_mol_C<-bonita$gpp_mol_O2/bonita$predictions

##Plot raw metabolism in mol O2

p1<-ggplot(data=bonita,aes(x=date, y=gpp_mol_O2))+
  geom_line(size=1.1)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text = element_text(size = 18))+
  theme(axis.title.x=element_text(size=14,colour = "black"))+
  theme(axis.title.y=element_text(size=14,colour = "black"))+
  theme(axis.text.y=element_text(size=12,colour = "black"))+
  theme(axis.text.x=element_text(size=12,colour = "black"))+
  ylab(expression(paste("GPP (mol ",O[2]," ", m^-2," ",d^-1,")")))+
  #ylim(c(0,10))+
  xlab("Date")+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
  scale_y_continuous(limits=c(0,0.4), breaks=seq(from=0, to=0.4, by=0.1))

##Estimate metabolism in C 
bonita$gpp_mol_C_lowpq<-bonita$gpp_mol_O2/0.4 ##low PQ literature range
bonita$gpp_mol_C_highpq<-bonita$gpp_mol_O2/3.5 ##high PQ literature range 
bonita$gpp_mol_C_1.2<-bonita$gpp_mol_O2/1.2 ##Common value used in freshwaters

p2<-ggplot(data=bonita)+
  geom_ribbon(data=bonita,aes(ymin=gpp_mol_C_lowpq,ymax=gpp_mol_C_highpq, x=date,fill="PQ=0.4-3.5") , alpha=0.6, linetype=0, color=NA)+
  geom_line(aes(x=date, y=gpp_mol_C,  linetype="PQ_Bio"),size=1.2)+
  #geom_line(data=bonita, aes(x=date, y=gpp_mol_C_lowpq),size=1.1)+
  #geom_line(data=bonita, aes(x=date, y=gpp_mol_C_highpq),size=1.1)+
  geom_line(data=bonita, aes(x=date, y=gpp_mol_C_1.2,linetype="PQ=1.2"),size=1.1 )+
  #facet_wrap(~site, ncol=1, labeller = as_labeller(site_name))+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text = element_text(size = 18),legend.position = "top")+
  theme(axis.title.x=element_text(size=14,colour = "black"))+
  theme(axis.title.y=element_text(size=14,colour = "black"))+
  theme(axis.text.y=element_text(size=12,colour = "black"))+
  theme(axis.text.x=element_text(size=12,colour = "black"))+
  #ylim(c(0,1))+
  xlab("Date")+
  ylab(expression(paste("GPP (mol C ", m^-2,d^-1,")")))+
  #scale_colour_manual(labels=c(expression(PQ['net']), "PQ=1.2"), name='', values=c("PQ_net" = "black", "PQ=1.2"="black"))+
  scale_fill_manual(name = '', values=c("PQ=0.4-3.5" = "gray"))+
  scale_linetype_manual(name='',labels=c(expression(PQ['Bio']), "PQ=1.2"),values=c("PQ_Bio"="solid","PQ=1.2"="dotted"))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),legend.text=element_text(size=13))+
  scale_y_continuous(limits=c(0,0.4), breaks=seq(from=0, to=0.4, by=0.1))



##Load metabolism data
mydata = read.csv("JAP_METABOLISM_PERKINS_2020.csv",header=T)
mydata<-mydata[,c(1:11)]
names(mydata)<-c("date", "discharge.cms", "depth.m", "width.m", "GPP", "ER","GPP.low","GPP.high","ER.low", "ER.high", "site")
mydata$date<-as.POSIXct(mydata$date, format="%m/%d/%Y")
mydata$gpp_mol_O2<-mydata$GPP/32
mydata$gpp_mol_O2[81]<-0.09
tem<-mydata$gpp_mol_O2
tem[82]<-0.088
perkins$gpp_mol_O2<-tem
perkins$gpp_mol_C<-perkins$gpp_mol_O2/perkins$predictions

perkins$gpp_mol_C_lowpq<-perkins$gpp_mol_O2/0.4
perkins$gpp_mol_C_highpq<-perkins$gpp_mol_O2/3.5
perkins$gpp_mol_C_1.2<-perkins$gpp_mol_O2/1.2


p3<-ggplot(data=perkins,aes(x=date, y=gpp_mol_O2))+
  geom_line(size=1.1)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text = element_text(size = 18))+
  theme(axis.title.x=element_text(size=14,colour = "black"))+
  theme(axis.title.y=element_text(size=14,colour = "black"))+
  theme(axis.text.y=element_text(size=12,colour = "black"))+
  theme(axis.text.x=element_text(size=12,colour = "black"))+
  ylab(expression(paste("GPP (mol ",O[2]," ", m^-2," ",d^-1,")")))+
  scale_y_continuous(limits=c(0,0.7), breaks=seq(from=0, to=0.7, by=0.2))+
  xlab("Date")

p4<-ggplot(data=perkins)+
  geom_ribbon(data=perkins,aes(ymin=gpp_mol_C_lowpq,ymax=gpp_mol_C_highpq, x=date,fill="PQ=0.4-3.5") , alpha=0.6, linetype=0, color=NA)+
  geom_line(aes(x=date, y=gpp_mol_C,  linetype="PQ_Bio"),size=1.2)+
  #geom_line(data=bonita, aes(x=date, y=gpp_mol_C_lowpq),size=1.1)+
  #geom_line(data=bonita, aes(x=date, y=gpp_mol_C_highpq),size=1.1)+
  geom_line(data=perkins, aes(x=date, y=gpp_mol_C_1.2,linetype="PQ=1.2"),size=1.1 )+
  #facet_wrap(~site, ncol=1, labeller = as_labeller(site_name))+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text = element_text(size = 18),legend.position = "top")+
  theme(axis.title.x=element_text(size=14,colour = "black"))+
  theme(axis.title.y=element_text(size=14,colour = "black"))+
  theme(axis.text.y=element_text(size=12,colour = "black"))+
  theme(axis.text.x=element_text(size=12,colour = "black"))+
  #ylim(c(0,1))+
  xlab("Date")+
  ylab(expression(paste("GPP (mol C ", m^-2,d^-1,")")))+
  #scale_colour_manual(labels=c(expression(PQ['net']), "PQ=1.2"), name='', values=c("PQ_net" = "black", "PQ=1.2"="black"))+
  scale_fill_manual(name = '', values=c("PQ=0.4-3.5" = "gray"))+
  scale_linetype_manual(name='',labels=c(expression(PQ['Bio']), "PQ=1.2"),values=c("PQ_Bio"="solid","PQ=1.2"="dotted"))+
  theme(legend.position = "none")+
  scale_y_continuous(limits=c(0,0.7), breaks=seq(from=0, to=0.7, by=0.2))

patchwork<-p1/p3 | p2/p4

patchwork + plot_annotation(tag_levels = list(c('A)', 'C)','B)','D)')))




