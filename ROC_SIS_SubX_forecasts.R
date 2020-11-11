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


# Hacemos los roc.plot por leads
x.w1 = dt.fcstSIS.bin[L>=1 & L<=7,]$bin.obs.num
prob.w1 = dt.fcstSIS.bin[L>=1 & L<=7,]$prob

x.w2 = dt.fcstSIS.bin[L>=8 & L<=14,]$bin.obs.num
prob.w2 = dt.fcstSIS.bin[L>=8 & L<=14,]$prob

x.w3 = dt.fcstSIS.bin[L>=15 & L<=21,]$bin.obs.num
prob.w3 = dt.fcstSIS.bin[L>=15 & L<=21,]$prob

x.w4 = dt.fcstSIS.bin[L>=22 & L<=28,]$bin.obs.num
prob.w4 = dt.fcstSIS.bin[L>=22 & L<=28,]$prob

#----------------------------------------------------****
# Ahora el evento negativo
umb=-0.75

dt.fcstSIS.bin.neg = dt.fcstSIS.bin.neg[SIS.fcst.ECCC>umb, SIS.fcst.ECCC:=0]
dt.fcstSIS.bin.neg = dt.fcstSIS.bin.neg[SIS.fcst.ECCC<=umb, SIS.fcst.ECCC:=1]

dt.fcstSIS.bin.neg = dt.fcstSIS.bin.neg[SIS.fcst.EMC>umb, SIS.fcst.EMC:=0]
dt.fcstSIS.bin.neg = dt.fcstSIS.bin.neg[SIS.fcst.EMC<=umb, SIS.fcst.EMC:=1]

dt.fcstSIS.bin.neg = dt.fcstSIS.bin.neg[SIS.fcst.ESRL>umb, SIS.fcst.ESRL:=0]
dt.fcstSIS.bin.neg = dt.fcstSIS.bin.neg[SIS.fcst.ESRL<=umb, SIS.fcst.ESRL:=1]

dt.fcstSIS.bin.neg = dt.fcstSIS.bin.neg[SIS.fcst.NRL>umb, SIS.fcst.NRL:=0]
dt.fcstSIS.bin.neg = dt.fcstSIS.bin.neg[SIS.fcst.NRL<=umb, SIS.fcst.NRL:=1]

dt.fcstSIS.bin.neg = dt.fcstSIS.bin.neg[SIS.fcst.GMAO.run1>umb, SIS.fcst.GMAO.run1:=0]
dt.fcstSIS.bin.neg = dt.fcstSIS.bin.neg[SIS.fcst.GMAO.run1<=umb, SIS.fcst.GMAO.run1:=1]

dt.fcstSIS.bin.neg = dt.fcstSIS.bin.neg[SIS.fcst.GMAO.run2>umb, SIS.fcst.GMAO.run2:=0]
dt.fcstSIS.bin.neg = dt.fcstSIS.bin.neg[SIS.fcst.GMAO.run2<=umb, SIS.fcst.GMAO.run2:=1]

dt.fcstSIS.bin.neg = dt.fcstSIS.bin.neg[SIS.fcst.RSMAS.run1>umb, SIS.fcst.RSMAS.run1:=0]
dt.fcstSIS.bin.neg = dt.fcstSIS.bin.neg[SIS.fcst.RSMAS.run1<=umb, SIS.fcst.RSMAS.run1:=1]

dt.fcstSIS.bin.neg = dt.fcstSIS.bin.neg[SIS.fcst.RSMAS.run2>umb, SIS.fcst.RSMAS.run2:=0]
dt.fcstSIS.bin.neg = dt.fcstSIS.bin.neg[SIS.fcst.RSMAS.run2<=umb, SIS.fcst.RSMAS.run2:=1]

# Calculo la prob de ocurrencia

dt.fcstSIS.bin.neg[,prob := rowMeans(.SD, na.rm=TRUE), by = c("targetdate","startdate","L"), .SDcols = mme.members]

bin.obs.SIS = dt.RTSIS
bin.obs.SIS = bin.obs.SIS[,bin.obs := SIS.index<=umb]
bin.obs.SIS$bin.obs.num=0
bin.obs.SIS[bin.obs==TRUE]$bin.obs.num = 1

dt.fcstSIS.bin.neg = merge(dt.fcstSIS.bin.neg,bin.obs.SIS,by="targetdate",all=TRUE)

# Hacemos los roc.plot por leads
x.w1.n = dt.fcstSIS.bin.neg[L>=1 & L<=7,]$bin.obs.num
prob.w1.n = dt.fcstSIS.bin.neg[L>=1 & L<=7,]$prob

x.w2.n = dt.fcstSIS.bin.neg[L>=8 & L<=14,]$bin.obs.num
prob.w2.n = dt.fcstSIS.bin.neg[L>=8 & L<=14,]$prob

x.w3.n = dt.fcstSIS.bin.neg[L>=15 & L<=21,]$bin.obs.num
prob.w3.n = dt.fcstSIS.bin.neg[L>=15 & L<=21,]$prob

x.w4.n = dt.fcstSIS.bin.neg[L>=22 & L<=28,]$bin.obs.num
prob.w4.n = dt.fcstSIS.bin.neg[L>=22 & L<=28,]$prob


# par(mfrow=c(2,4))

par(pty='s',las=1)
roc.plot(x.w1,prob.w1,main=paste0("SIS>0.75 d1-7. ROCA= ",as.character(round(roc.area(x.w1,prob.w1)$A,2))))
par(pty='s',las=1)
roc.plot(x.w2,prob.w2,main=paste0("SIS>0.75 d8-14. ROCA= ",as.character(round(roc.area(x.w2,prob.w2)$A,2))))
par(pty='s',las=1)
roc.plot(x.w2,prob.w3,main=paste0("SIS>0.75 d15-21. ROCA= ",as.character(round(roc.area(x.w3,prob.w3)$A,2))))


roc.plot(x.w2,prob.w4,main=paste0("d) SIS>0.75 d22-28. ROCA= ",as.character(round(roc.area(x.w4,prob.w4)$A,2))))



roc.plot(x.w1.n,prob.w1.n,main=paste0("e) SIS<-0.75 d1-7. ROCA= ",as.character(round(roc.area(x.w1.n,prob.w1.n)$A,2))))
roc.plot(x.w2.n,prob.w2.n,main=paste0("f) SIS<-0.75 d8-14. ROCA= ",as.character(round(roc.area(x.w2.n,prob.w2.n)$A,2))))
roc.plot(x.w2.n,prob.w3.n,main=paste0("g) SIS<-0.75 d15-21. ROCA= ",as.character(round(roc.area(x.w3.n,prob.w3.n)$A,2))))
roc.plot(x.w2.n,prob.w4.n,main=paste0("h) SIS<-0.75 d22-28. ROCA= ",as.character(round(roc.area(x.w4.n,prob.w4.n)$A,2))))




# Este grid.arrange no funciona pq no son ggplots tengo que buscar el alternativo. Esto lo habia hecho en el trabajo con Caio pero 
# no puedo acceder a Kyle
# fig <- grid.arrange(g1,g2,g3,g4, ncol = 1,top = textGrob(paste0("ROC curves for daily SIS index SubX EM forecasts "),gp=gpar(fontsize=13,font=3)))
# ggsave(filename=paste0("/home/alvarez/SubX_processed_Rdata/SIS_index_forecast_ROC.png"),plot=fig,width = 11, height = 14)

