
# Load libraries ----------------------------------------------------------

library("lattice")
library("viridis")
library(RColorBrewer)
library(patchwork)
library(gridExtra)
library(reshape2); library(dplyr); library(tidyr)
library(cowplot)


# Estimating the effects of photorespiration on the PQ --------------------

###Data plucked from Burris 1981
setwd("~/GitHub/metabolism-PQ")
temp<-20#C
pres<-1#atm
data<-read.csv("burris_pr_change.csv")
data$PQ.delta<-data$PQ-1.2
data$O2.pp.kpa<-data$O2.pp*101
data$O2.mgL<-data$O2.pp.kpa*0.36610

###Transform oxy concentration to percent saturation
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


# Set up NO3 effects-to be calculated for each panel ----------------------

##Create NO3 proportion of DIN. Proportion from 0 to 1
DIN<-seq(from=0, to=1, length.out=400)


# Estimate effects for panel A where C:N=5 --------------------------------

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


##Plot heatmap of change in PQ
p1<-levelplot(Z ~ X*Y, data=act.data,scales=list(x=list(at=seq(40,180,20), cex=1.3),tck = c(1,0), y = list(at=seq(0,1,.2),cex=1.3)),
               xlab=list(label=expression(paste("Maximum Diel DO (",O[2]," ","%"," ",saturation,")")), cex=1.4),#list(draw = FALSE),
               ylab=list(expression(paste("Proportion of DIN present as ",NO[3],)), cex=1.4),
          main=list(label="C:N=5", cex=1.5),
          col.regions =colorRampPalette(rev(c("#332288", "#5f48cb","#785EF0", "#DC267F", "#FE6100","#FFB000")))(100),
          #col.regions =colorRampPalette(rev(c("#4575B4", "#91BFDB", "#E0F3F8", "#FFFFBF","#FEE090" , "#FC8D59","#D73027" )))(100),
          at=seq(from=min(act.data$Z), to=max(act.data$Z),  length.out=11),
          aspect="fill",
          panel = function(...){
            panel.levelplot(...)
            panel.abline(coef = coef(lm(c(0.0,0.53)~c(85,102))), lwd=4,lty = "dotted")},
          colorkey=list(labels=list(cex=1.3,at=seq(0.8,-1.2,-.2))))
   


# Estimate effects for panel B where C:N=20 -------------------------------


##Create NO3 proportion of DIN. Proportion from 0 to 1
DIN<-seq(from=0, to=1, length.out=400)

#Calculate PQ change due to NO3
alg_CN<-20 #must define algal C:N ratio  
y<-DIN*2/alg_CN

##Create O2 and PQ gradient (based on Burris 1981)
#O2<-seq(from=90, to=150, length.out=400)
#PQ<-seq(from=1.2, to=0.5, length.out=400)

#Calculate the change in PQ due to photorespiration assume PQc is 1.1 and PQn is 0 (???)



##Expand to all combinations
act.data<-expand.grid(X=O2,Y= DIN)
z<-expand.grid(PQ, y)

#add two individual effects together
act.data$Z<-z$Var1+z$Var2



##Plot heatmap of change in PQ
p2<-levelplot(Z ~ X*Y, data=act.data,scales=list(x=list(at=seq(40,180,20), cex=1.3),tck = c(1,0), y = list(at=seq(0,1,.2),cex=1.3)),
               xlab=list(label=expression(paste("Maximum Diel DO (",O[2]," ","%"," ",saturation,")")), cex=1.4),#list(draw = FALSE),
               ylab=list(expression(paste("Proportion ",NO[3]," of DIN")), cex=1.4),
                    main=list(label="C:N=20", cex=1.5), 
               col.regions =colorRampPalette(rev(c( "#5f48cb","#785EF0", "#DC267F", "#FE6100","#FFB000")))(100),
               #col.regions = colorRampPalette(rev(c( "#91BFDB","#E0F3F8", "#FFFFBF","#FEE090" , "#FC8D59","#D73027" )))(100),
              at=seq(from=min(act.data$Z), to=max(act.data$Z),  length.out=8),
          aspect="fill",
              panel = function(...){
                panel.levelplot(...)
                panel.abline(coef = coef(lm(c(0,1)~c(85,95))), lwd=4,lty = "dotted")},
          colorkey=list(labels=list(cex=1.3,at=seq(0.8,-1.2,-.2))))     
      


trellis.device(windows, height=6, width=7)

grid.arrange(p1, p2)
print(p1, split=c(2,1,1,1),panel.height=list(x=2, units="in", y=20), more=TRUE)
print(p2, split=c(1,1,2,1),panel.height=list(x=2.05, units="in"))

plots<-list(p1,p2)
grobs <- lapply(plots, as_grob)
plot_widths <- lapply(grobs, function(x) {grobs[[1]][["list"]][["plot"]][["y.limits"]]})
aligned_widths <- align_margin(plot_widths, "first", greedy=TRUE)
aligned_widths <- align_margin(aligned_widths, "last", greedy=TRUE)


for (i in seq_along(plots)) {
  grobs[[i]]$heights <- aligned_widths[[i]]
}

plot_grid(plotlist = grobs, ncol = 1)


