# This scripts plots all the SIS index forecasts for the SubX models, the MME forecast and the
# MME spread (using all SubX models' ensemble mean)
#
# M. Alvarez - 2020
#_____________________________________________________________________________________________
rm(list=ls())
graphics.off()

setwd("/home/alvarez/")

# Call libraries to be used
library("ncdf4")
library("metR")
library('pracma')
library('lubridate')
library('reshape2')
library('data.table')
library("s2dverification") # No se si lo voy a usar al final
library("RColorBrewer")
library("maps")
library("ggplot2")
library("gridExtra")
library("grid")
library("verification")


# Slightly changed function for reliability plot from verification package

# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
# ** Copyright UCAR (c) 1992 - 2004 
# ** University Corporation for Atmospheric Research(UCAR) 
# ** National Center for Atmospheric Research(NCAR) 
# ** Research Applications Program(RAP) 
# ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA 
# ** 2004/1/7 11:29:42 
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
reliability.plot.2<- function(x.w,prob.w, titl = NULL, legend.names = NULL, ...){
  ## this function is similar to a attribute plot but somewhat simplified.
  nbins = 10
  bk = seq(0,1,1/nbins) #breaks
  
  obar = mean(x.w,na.rm=TRUE)
  
  h = hist(prob.w,breaks=bk,plot=F)$counts
  g = hist(prob.w[x.w==1],breaks=bk,plot=F)$counts # Forecast probabilities for the event when the event (SIS>umb) is observed
  
  obar.i = g/h 
  # E.g. (h) 46 times the event was forecast to ocurr with a 90-100% probability.
  #      (g) 43 times the event was forecast to ocurr with a 90-100% probability when IT ACTUALLY OCURRED  
  # observed frequency: 43/46 = 0.9348
  # forecast probability: 90-100% (considered 95% for plotting)
  
  x = seq((1/nbins)/2,1,1/nbins)

  # methods(  
  old.par <- par(no.readonly = TRUE) # all par settings which
  # could be changed.
  on.exit(par(old.par))
  
  obar.i<- as.matrix(obar.i)
  if(is.null(legend.names)) legend.names<- paste("Model", seq(1,dim(obar.i)[2])) 
  
  prob.y<- as.matrix(h)
  
  plot(x, obar.i[,1],  col = 2, lwd = 2, type = "n",
       xlim = c(0,1), ylim = c(0,1),
       xlab =  expression( paste("Forecast probability, ", y[i] ) ),
       ylab = expression( paste("Observed relative frequency, ", bar(o)[1] ))
  )
  if(is.null(titl)){title("Reliability Plot")}else{
    title(titl)
  }
  
  m<- dim(obar.i)[2]
  for(i in 1:m){
    points(x, obar.i[,i], type = "b", col = 1+i, lty = i, lwd = 2)
  }
  abline(0,1)
  abline(h=obar,lty="dashed")
  
  ## rank histogram plot in lower corner.
  
  pp<- par("plt")
  
  # plot lower box plot.
  
  par("plt" = c(pp[2] - 0.2 , pp[2],  pp[3], pp[3]+ 0.2) )
  par(new = TRUE)
  barplot(prob.y[,1]/sum(prob.y), axes = FALSE, axisnames = FALSE,ylim=c(0,0.6))
  axis(4)
  box() 
  
  invisible()  
}# close function






#---------------------------------------------------------------------------------------
#  Main Program  
#---------------------------------------------------------------------------------------
load("/home/alvarez/SubX_scripts/SubX_SIS_index_forecast_MME.rda")
dt.fcstSIS.MME = dt.fcstSIS
setnames(dt.fcstSIS.MME,old="SIS.index.fcst",new="SIS.fcst.MME")

rm(dt.fcstSIS)

load("/home/alvarez/SubX_scripts/SubX_SIS_index_forecast_ECCC.rda")
dt.fcstSIS.ECCC = dt.fcstSIS
setnames(dt.fcstSIS.ECCC,old="SIS.index.fcst",new="SIS.fcst.ECCC")
rm(dt.fcstSIS)
load("/home/alvarez/SubX_scripts/SubX_SIS_index_forecast_EMC.rda")
dt.fcstSIS.EMC = dt.fcstSIS
setnames(dt.fcstSIS.EMC,old="SIS.index.fcst",new="SIS.fcst.EMC")
rm(dt.fcstSIS)
load("/home/alvarez/SubX_scripts/SubX_SIS_index_forecast_ESRL.rda")
dt.fcstSIS.ESRL = dt.fcstSIS
setnames(dt.fcstSIS.ESRL,old="SIS.index.fcst",new="SIS.fcst.ESRL")
rm(dt.fcstSIS)
load("/home/alvarez/SubX_scripts/SubX_SIS_index_forecast_NRL.rda")
dt.fcstSIS.NRL = dt.fcstSIS
setnames(dt.fcstSIS.NRL,old="SIS.index.fcst",new="SIS.fcst.NRL")
rm(dt.fcstSIS)

load("/home/alvarez/SubX_scripts/SubX_SIS_index_forecast_GMAO.rda")
dt.fcstSIS.GMAO = dt.fcstSIS
setnames(dt.fcstSIS.GMAO,old="SIS.index.fcst.run1",new="SIS.fcst.GMAO.run1")
setnames(dt.fcstSIS.GMAO,old="SIS.index.fcst.run2",new="SIS.fcst.GMAO.run2")
rm(dt.fcstSIS)
load("/home/alvarez/SubX_scripts/SubX_SIS_index_forecast_RSMAS.rda")
dt.fcstSIS.RSMAS = dt.fcstSIS
setnames(dt.fcstSIS.RSMAS,old="SIS.index.fcst.run1",new="SIS.fcst.RSMAS.run1")
setnames(dt.fcstSIS.RSMAS,old="SIS.index.fcst.run2",new="SIS.fcst.RSMAS.run2")
rm(dt.fcstSIS)

# Cargo el RTSIS para la verificación
load("~/SubX_scripts/RTSIS_OA9899_1415.rda")


dt.fcstSIS=merge(dt.fcstSIS.ECCC,dt.fcstSIS.EMC,by=c("targetdate","startdate","L"),all=TRUE)
dt.fcstSIS=merge(dt.fcstSIS,dt.fcstSIS.ESRL,by=c("targetdate","startdate","L"),all=TRUE)
dt.fcstSIS=merge(dt.fcstSIS,dt.fcstSIS.NRL,by=c("targetdate","startdate","L"),all=TRUE)
dt.fcstSIS=merge(dt.fcstSIS,dt.fcstSIS.GMAO,by=c("targetdate","startdate","L"),all=TRUE)
dt.fcstSIS=merge(dt.fcstSIS,dt.fcstSIS.RSMAS,by=c("targetdate","startdate","L"),all=TRUE)
# dt.fcstSIS=merge(dt.fcstSIS,dt.fcstSIS.MME,by=c("targetdate","startdate","L"),all=TRUE)


# data should have in the first column
# roc.plot(x,pred): x is a binary observation (0,1), pred are forecast probabilities of event
#compute ROC area using the verification package
# roc.area(data[,2],data[,1])

bin.obs.SIS = dt.RTSIS
umb=0.75
bin.obs.SIS = bin.obs.SIS[,bin.obs := SIS.index>=umb]
bin.obs.SIS$bin.obs.num=0
bin.obs.SIS[bin.obs==TRUE]$bin.obs.num = 1

# Determine probability of event
mme.members=c("SIS.fcst.ECCC","SIS.fcst.EMC","SIS.fcst.ESRL","SIS.fcst.NRL","SIS.fcst.GMAO.run1","SIS.fcst.GMAO.run2","SIS.fcst.RSMAS.run1","SIS.fcst.RSMAS.run2")

# dt.fcstSIS[,prob := max(.SD, na.rm=TRUE), by = c("targetdate","startdate","L"), .SDcols = mme.members]

dt.fcstSIS.bin.neg = copy(dt.fcstSIS)
# Insert operation to break copy reference and not to overwrite
dt.fcstSIS.bin.neg$SIS.fcst.GMAO.run2[1] <- NA_real_  # new operation

dt.fcstSIS.bin = copy(dt.fcstSIS)

dt.fcstSIS.bin = dt.fcstSIS[SIS.fcst.ECCC>=umb, SIS.fcst.ECCC:=1]
dt.fcstSIS.bin = dt.fcstSIS[SIS.fcst.ECCC<umb, SIS.fcst.ECCC:=0]

dt.fcstSIS.bin = dt.fcstSIS.bin[SIS.fcst.EMC>=umb, SIS.fcst.EMC:=1]
dt.fcstSIS.bin = dt.fcstSIS.bin[SIS.fcst.EMC<umb, SIS.fcst.EMC:=0]

dt.fcstSIS.bin = dt.fcstSIS.bin[SIS.fcst.ESRL>=umb, SIS.fcst.ESRL:=1]
dt.fcstSIS.bin = dt.fcstSIS.bin[SIS.fcst.ESRL<umb, SIS.fcst.ESRL:=0]

dt.fcstSIS.bin = dt.fcstSIS.bin[SIS.fcst.NRL>=umb, SIS.fcst.NRL:=1]
dt.fcstSIS.bin = dt.fcstSIS.bin[SIS.fcst.NRL<umb, SIS.fcst.NRL:=0]

dt.fcstSIS.bin = dt.fcstSIS.bin[SIS.fcst.GMAO.run1>=umb, SIS.fcst.GMAO.run1:=1]
dt.fcstSIS.bin = dt.fcstSIS.bin[SIS.fcst.GMAO.run1<umb, SIS.fcst.GMAO.run1:=0]

dt.fcstSIS.bin = dt.fcstSIS.bin[SIS.fcst.GMAO.run2>=umb, SIS.fcst.GMAO.run2:=1]
dt.fcstSIS.bin = dt.fcstSIS.bin[SIS.fcst.GMAO.run2<umb, SIS.fcst.GMAO.run2:=0]

dt.fcstSIS.bin = dt.fcstSIS.bin[SIS.fcst.RSMAS.run1>=umb, SIS.fcst.RSMAS.run1:=1]
dt.fcstSIS.bin = dt.fcstSIS.bin[SIS.fcst.RSMAS.run1<umb, SIS.fcst.RSMAS.run1:=0]

dt.fcstSIS.bin = dt.fcstSIS.bin[SIS.fcst.RSMAS.run2>=umb, SIS.fcst.RSMAS.run2:=1]
dt.fcstSIS.bin = dt.fcstSIS.bin[SIS.fcst.RSMAS.run2<umb, SIS.fcst.RSMAS.run2:=0]

# Calculo la prob de ocurrencia

dt.fcstSIS.bin[,prob := rowMeans(.SD, na.rm=TRUE), by = c("targetdate","startdate","L"), .SDcols = mme.members]

bin.obs.SIS$SIS.index=NULL
bin.obs.SIS$bin.obs=NULL

dt.fcstSIS.bin = merge(dt.fcstSIS.bin,bin.obs.SIS,by="targetdate",all=TRUE)

# Compute reliability diagrams
# Decide number of bins
nbins = 10
bk = seq(0,1,1/nbins) #breaks

# Hacemos los roc.plot por leads
x.w1 = dt.fcstSIS.bin[L>=1 & L<=7,]$bin.obs.num
prob.w1 = dt.fcstSIS.bin[L>=1 & L<=7,]$prob

x.w2 = dt.fcstSIS.bin[L>=8 & L<=14,]$bin.obs.num
prob.w2 = dt.fcstSIS.bin[L>=8 & L<=14,]$prob

x.w3 = dt.fcstSIS.bin[L>=15 & L<=21,]$bin.obs.num
prob.w3 = dt.fcstSIS.bin[L>=15 & L<=21,]$prob

x.w4 = dt.fcstSIS.bin[L>=22 & L<=28,]$bin.obs.num
prob.w4 = dt.fcstSIS.bin[L>=22 & L<=28,]$prob

# # Me queda
# # remover los NA de obs
# 
# h = hist(prob.w1,breaks=bk,plot=F)$counts
# g = hist(prob.w1[x.w1==1],breaks=bk,plot=F)$counts # Forecast probabilities for the event when the event (SIS>umb) is observed
# 
# obari = g/h 
# # E.g. (h) 46 times the event was forecast to ocurr with a 90-100% probability.
# #      (g) 43 times the event was forecast to ocurr with a 90-100% probability when IT ACTUALLY OCURRED  
# # observed frequency: 43/46 = 0.9348
# # forecast probability: 90-100% (considered 95% for plotting)
# 
# yi = seq((1/nbins)/2,1,1/nbins)
# obar = mean(x.w1,na.rm=TRUE)



par(pty='s',las=1)
reliability.plot.2(x.w1,prob.w1,titl=paste0("Week 1. SIS>",as.character(umb)),legend.names="")
par(pty='s',las=1)
reliability.plot.2(x.w2,prob.w2,titl=paste0("Week 2. SIS>",as.character(umb)),legend.names="")
par(pty='s',las=1)
reliability.plot.2(x.w3,prob.w3,titl=paste0("Week 3. SIS>",as.character(umb)),legend.names="")
par(pty='s',las=1)
reliability.plot.2(x.w4,prob.w4,titl=paste0("Week 4. SIS>",as.character(umb)),legend.names="")


