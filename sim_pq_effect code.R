library(ggplot2)
library(chron) #this is how we used to do time.  Unwilling to change for this exercise
library(ggplot2)
library(patchwork)
library(reshape2); library(dplyr); library(tidyr)


PQ_n<-seq(from=0.05, to=5, by=.01)
PQ_d<-rep(1, times=length(PQ_n))
PQ<-(PQ_n/PQ_d)
#PQ<- seq(from=0.1, to=3.5, length.out=length(GPP.O.d))
#GPP.C<-GPP.O*mr/PQ
data<-data.frame(PQ,PQ_n )

p1<-ggplot(data=data,aes(x=PQ_n, y=PQ))+
  geom_line(size=1.5)+
  theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=14,color="black"))+
  theme(axis.title.y=element_text(size=14,color="black"))+
  theme(axis.text.y=element_text(size=11,color="black"))+
  theme(axis.text.x=element_text(size=11,color="black"))+
  xlab("Change in Oxygen")+
  ylab("PQ")
  #scale_y_continuous(limits=c(-250,800), breaks=seq(from=-250, to=800, by=150))+
  #scale_x_continuous(limits=c(0,3.5), breaks=seq(from=0, to=3.6, by=.4))+
  #geom_hline(yintercept=0,linetype="dashed")


PQ_d<-seq(from=0.05, to=5, by=.01)
PQ_n<-rep(1, times=length(PQ_d))
PQ<-(PQ_n/PQ_d)
#PQ<- seq(from=0.1, to=3.5, length.out=length(GPP.O.d))
#GPP.C<-GPP.O*mr/PQ
#GPP.C<-GPP.O*mr/PQ
data<-data.frame(PQ,PQ_d )

p2<-ggplot(data=data,aes(x=PQ_d, y=PQ))+
  geom_line(size=1.5)+
  theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=14,color="black"))+
  theme(axis.title.y=element_text(size=14,color="black"))+
  theme(axis.text.y=element_text(size=11,color="black"))+
  theme(axis.text.x=element_text(size=11,color="black"))+
  xlab("Change in DIC")+
  ylab("PQ")
#scale_y_continuous(limits=c(-250,800), breaks=seq(from=-250, to=800, by=150))+
#scale_x_continuous(limits=c(0,3.5), breaks=seq(from=0, to=3.6, by=.4))+
#geom_hline(yintercept=0,linetype="dashed")

mr<-14/32
PQ<- seq(from=0.1, to=3.5, by=0.01)
GPP.O<-rep(250, times=length(PQ))
GPP.O.2<-rep(19, times=length(PQ))
GPP.C<-GPP.O/PQ
GPP.C.2<-GPP.O.2/PQ
data<-data.frame(PQ, GPP.C, GPP.O, GPP.C.2 )

p3<-ggplot(data=data,aes(x=PQ, y=(GPP.C-GPP.O), color="High rate"))+
  geom_line(size=1.5)+
  geom_line(aes(x=PQ, y=GPP.C.2-GPP.O.2, color="Low rate"),size=1.5)+
  theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=14,color="black"))+
  theme(axis.title.y=element_text(size=14,color="black"))+
  theme(axis.text.y=element_text(size=11,color="black"))+
  theme(axis.text.x=element_text(size=11,color="black"))+
  xlab("PQ")+
  ylab(expression(paste(,GPP[C],"-",GPP[O],)))+
  scale_y_continuous(limits=c(-200,2250), breaks=seq(from=-200, to=2250, by=200))+
  scale_x_continuous(limits=c(0,3.5), breaks=seq(from=0, to=3.6, by=.4))+
  geom_hline(yintercept=0,linetype="dashed")+
  scale_color_manual(values=c("High rate"="#000000","Low rate"="#D55E00"))+
  guides(color=guide_legend(title=element_blank()))+
  theme(legend.position='top')



p1/p2/p3



PQ_n<-seq(from=5, to=0.05, by=-.1)
PQ_d<-seq(from=0.05, to=5, by=.1)
PQ<-PQ_n/PQ_d

data<-data.frame(PQ_d,PQ_n, PQ )
#p4<-
ggplot(data=data,aes(x=PQ_d, y=PQ_n))+
  geom_point(aes(size=PQ))+
  theme_classic()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_text(size=14,color="black"))+
  theme(axis.title.y=element_text(size=14,color="black"))+
  theme(axis.text.y=element_text(size=11,color="black"))+
  theme(axis.text.x=element_text(size=11,color="black"))+
  xlab("Change in DIC")+
  ylab("Change in Oxygen")
